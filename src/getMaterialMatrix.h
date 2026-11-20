/** @file getMaterialMatrix.h

    @brief Provides a simple class to obtain a material matrix pointer

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsMaterialMatrixComposite.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{


/**
 * @brief Plane-stress stiffness matrix of an orthotropic ply in its material axes (x,y)
 *
 * Returns the 3x3 matrix G relating (sigma_xx, sigma_yy, sigma_xy) to
 * (eps_xx, eps_yy, gamma_xy):
 * \f[
 *   G = \frac{1}{1-\nu_{xy}\nu_{yx}}
 *       \begin{bmatrix} E_{xx} & \nu_{yx}E_{xx} & 0 \\
 *                       \nu_{xy}E_{yy} & E_{yy} & 0 \\
 *                       0 & 0 & (1-\nu_{xy}\nu_{yx})G_{xy} \end{bmatrix}.
 * \f]
 * G is symmetric only if the Poisson ratios satisfy the reciprocity relation
 * \f$\nu_{yx}E_{xx} = \nu_{xy}E_{yy}\f$ (Maxwell–Betti). Datasheet values are
 * usually rounded, so the relation is enforced up to a relative mismatch of 1%:
 * within it, both off-diagonal entries are set to their mean and G is exactly
 * symmetric; beyond it, the input is inconsistent (typically \f$\nu_{xy}=\nu_{yx}\f$
 * with \f$E_{xx}\neq E_{yy}\f$) and the function throws.
 *
 * @param Exx   Young's modulus in the x (fibre) direction
 * @param Eyy   Young's modulus in the y (transverse) direction
 * @param Gxy   In-plane shear modulus
 * @param nuxy  Poisson ratio: contraction in y under stress in x
 * @param nuyx  Poisson ratio: contraction in x under stress in y
 * @throws std::runtime_error if \f$|\nu_{yx}E_{xx}-\nu_{xy}E_{yy}| > 10^{-2}\max(|\nu_{yx}E_{xx}|,|\nu_{xy}E_{yy}|)\f$
 */
template<class T>
gsMatrix<T> gsCompositeMatrix(  const T Exx,
                              const T Eyy,
                              const T Gxy,
                              const T nuxy,
                              const T nuyx
                            )
{
    const T a = nuyx*Exx;
    const T b = nuxy*Eyy;
    GISMO_ENSURE(math::abs(a-b) <= (T)1e-2 * math::max(math::abs(a),math::abs(b)),
            "No symmetry in material properties for ply! (nu12*E2!=nu21*E1 within 1%):\n"<<
            "\tnu12 = "<<nuxy<<"\t E22 = "<<Eyy<<"\t nu12*E22 = "<<b<<"\n"
          <<"\tnu21 = "<<nuyx<<"\t E11 = "<<Exx<<"\t nu21*E11 = "<<a);

    gsMatrix<T> G(3,3);
    G.setZero();

    G(0,0) = Exx / (1-nuxy*nuyx);
    G(1,1) = Eyy / (1-nuxy*nuyx);
    G(2,2) = Gxy;
    G(0,1) = G(1,0) = (T)0.5 * (a+b) / (1-nuxy*nuyx);
    G(2,0) = G(0,2) = G(2,1) = G(1,2) = 0.0;
    return G;
}

template<class T>
gsMatrix<T> gsCompositeMatrixRaw(   const T G11,
                                    const T G22,
                                    const T G33,
                                    const T G12,
                                    const T G13,
                                    const T G23,
                                    const T G21,
                                    const T G31,
                                    const T G32 )
{
    gsMatrix<T> G(3,3);
    G.setZero();

    G(0,0) = G11;
    G(1,1) = G22;
    G(2,2) = G33;
    G(1,0) = G12;
    G(0,1) = G21;
    G(0,2) = G13;
    G(2,0) = G31;
    G(1,2) = G23;
    G(2,1) = G32;
    return G;
}

template<class T>
gsMatrix<T> gsCompositeMatrixRaw(   const T G11,
                                    const T G22,
                                    const T G33,
                                    const T G12,
                                    const T G13,
                                    const T G23 )
{
    gsMatrix<T> G(3,3);
    G.setZero();

    G(0,0) = G11;
    G(1,1) = G22;
    G(2,2) = G33;
    G(1,0) = G(0,1) = G12;
    G(0,2) = G(2,0) = G13;
    G(1,2) = G(2,1) = G23;
    return G;
}

/**
 * @brief      Gets a material matrix based on \a options
 *
 * @param[in]  mp          The undeformed geometry
 * @param[in]  thickness   The thickness
 * @param[in]  parameters  The parameters
 * @param[in]  rho         The density
 * @param[in]  options     The option list
 *
 * @tparam     d           The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T           Real type
 *
 * @return     The material matrix.
 *
 * @ingroup    KLShell
 *
 */
template<short_t d, class T>
typename gsMaterialMatrixBase<T>::uPtr
getMaterialMatrix(
                    const gsMultiPatch<T>               & mp,
                    const gsFunctionSet<T>                 & thickness,
                    const std::vector<gsFunctionSet<T> *>  & parameters,
                    const gsFunctionSet<T>                 & rho,
                    const gsOptionList                  & options
                    )
{
    index_t material        = options.askInt("Material",0); // SvK by default
    bool compressibility    = options.askSwitch("Compressibility",false);
    index_t implementation  = options.askInt("Implementation",1);

    enum Material       mat     = static_cast<enum Material>(material);
    enum Implementation impl    = static_cast<enum Implementation>(implementation);

    if      (mat==Material::SvK)
    {
        if (impl==Implementation::Composite)
            GISMO_ERROR("Construct composite material models using the constructor of gsMaterialMatrixComposite directly.");
        else
            return memory::make_unique(new gsMaterialMatrixLinear<d,T>(mp,thickness,*parameters.at(0),*parameters.at(1),rho));
    }
    else
    {
        if (impl==Implementation::Composite)
            GISMO_ERROR("Hyperelastic composites not available");
        //--------------------------------------NH Incompressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters,rho));
        }
        //--------------------------------------NH Compressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        //
        //--------------------------------------NH_ext Incompressible--------------------------------------------------
        else if ((mat==Material::NH_ext) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Extended NH with not available");
        }
        //---------------------------------------NH_ext Compressible-------------------------------------------------
        else if ((mat==Material::NH_ext) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH_ext, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH_ext, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH_ext, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        //
        //--------------------------------------MR Incompressible--------------------------------------------------
        else if ((mat==Material::MR) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters,rho));
        }
        //---------------------------------------MR Compressible-------------------------------------------------
        else if      ((mat==Material::MR) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        //
        //--------------------------------------OG Incompressible--------------------------------------------------
        else if ((mat==Material::OG) && (impl==Implementation::Analytical) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Ogden with analytical implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Generalized) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Ogden with generalized implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters,rho));
        }
        //---------------------------------------OG Compressible-------------------------------------------------
        else if      ((mat==Material::OG) && (impl==Implementation::Analytical) && (compressibility))
        {
            GISMO_ERROR("Compressible Ogden with analytical implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Generalized) && (compressibility))
        {
            GISMO_ERROR("Compressible Ogden with generalized implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters,rho));
        }
        else
            return NULL;
    }
}

/**
 * @brief      Gets a material matrix based on \a options
 *
 * @param[in]  mp          The undeformed geometry
 * @param[in]  thickness   The thickness
 * @param[in]  parameters  The parameters
 * @param[in]  options     The option list
 *
 * @tparam     d           The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T           Real type
 *
 * @return     The material matrix.
 *
 * @ingroup    KLShell
 *
 */
template<short_t d, class T>
typename gsMaterialMatrixBase<T>::uPtr
getMaterialMatrix(
                                const gsMultiPatch<T>               & mp,
                                const gsFunctionSet<T>                 & thickness,
                                const std::vector<gsFunctionSet<T> *>  & parameters,
                                const gsOptionList                  & options
                                )
{
    index_t material        = options.askInt("Material",0); // SvK by default
    bool compressibility    = options.askSwitch("Compressibility",false);
    index_t implementation  = options.askInt("Implementation",1);

    enum Material       mat     = static_cast<enum Material>(material);
    enum Implementation impl    = static_cast<enum Implementation>(implementation);

    if      (mat==Material::SvK)
    {
        if (impl==Implementation::Composite)
                GISMO_ERROR("Construct composite material models using the constructor of gsMaterialMatrixComposite directly.");
        else
                return memory::make_unique(new gsMaterialMatrixLinear<d,T>(mp,thickness,*parameters.at(0),*parameters.at(1)));
    }
    else
    {
        if (impl==Implementation::Composite)
            GISMO_ERROR("Hyperelastic composites not available");
        //--------------------------------------NH Incompressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters));
        }
        //--------------------------------------NH Compressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        //
        //--------------------------------------NH_ext Incompressible--------------------------------------------------
        else if ((mat==Material::NH_ext) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Extended NH with not available");
        }
        //---------------------------------------NH_ext Compressible-------------------------------------------------
        else if ((mat==Material::NH_ext) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH_ext, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH_ext, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::NH_ext, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        //
        //--------------------------------------MR Incompressible--------------------------------------------------
        else if ((mat==Material::MR) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters));
        }
        //---------------------------------------MR Compressible-------------------------------------------------
        else if      ((mat==Material::MR) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        //
        //--------------------------------------OG Incompressible--------------------------------------------------
        else if ((mat==Material::OG) && (impl==Implementation::Analytical) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Ogden with analytical implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Generalized) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Ogden with generalized implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,false>(mp,thickness,parameters));
        }
        //---------------------------------------OG Compressible-------------------------------------------------
        else if      ((mat==Material::OG) && (impl==Implementation::Analytical) && (compressibility))
        {
            GISMO_ERROR("Compressible Ogden with analytical implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Generalized) && (compressibility))
        {
            GISMO_ERROR("Compressible Ogden with generalized implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr index_t id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return memory::make_unique(new gsMaterialMatrixNonlinear<d,real_t,id,true>(mp,thickness,parameters));
        }
        else
            return NULL;
    }
}



} // namespace

