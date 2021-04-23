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

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixComposite.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{


template<class T>
gsMatrix<T> gsCompositeMatrix(  const T Exx,
                              const T Eyy,
                              const T Gxy,
                              const T nuxy,
                              const T nuyx
                            )
{
    GISMO_ASSERT(nuyx*Exx == nuxy*Eyy, "No symmetry in material properties for ply! (nu12*E2!=nu21*E1):\n"<<
            "\tnu12 = "<<nuxy<<"\t E22 = "<<Eyy<<"\t nu12*E22 = "<<nuxy*Eyy<<"\n"
          <<"\tnu21 = "<<nuyx<<"\t E11 = "<<Exx<<"\t nu21*E11 = "<<nuyx*Exx);

    gsMatrix<T> G(3,3);
    G.setZero();

    G(0,0) = Exx / (1-nuxy*nuyx);
    G(1,1) = Eyy / (1-nuxy*nuyx);
    G(2,2) = Gxy;
    G(0,1) = nuyx*Exx / (1-nuxy*nuyx);
    G(1,0) = nuxy*Eyy / (1-nuxy*nuyx);
    G(2,0) = G(0,2) = G(2,1) = G(1,2) = 0.0;
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
 * @ingroup    MaterialMatrix
 *
 */
template<short_t d, class T>
gsMaterialMatrixBase<T> * getMaterialMatrix(
                                const gsMultiPatch<T>               & mp,
                                const gsFunction<T>                 & thickness,
                                const std::vector<gsFunction<T> *>  & parameters,
                                const gsFunction<T>                 & rho,
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
            return new gsMaterialMatrixLinear<d,T>(mp,thickness,parameters,rho);
    }
    else
    {
        if (impl==Implementation::Composite)
            GISMO_ERROR("Hyperelastic composites not available");
        //--------------------------------------NH Incompressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters,rho);
        }
        //--------------------------------------NH Compressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
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
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
        }
        //
        //--------------------------------------MR Incompressible--------------------------------------------------
        else if ((mat==Material::MR) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters,rho);
        }
        //---------------------------------------MR Compressible-------------------------------------------------
        else if      ((mat==Material::MR) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
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
            constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters,rho);
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
            constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters,rho);
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
 * @ingroup    MaterialMatrix
 *
 */
template<short_t d, class T>
gsMaterialMatrixBase<T> * getMaterialMatrix(
                                const gsMultiPatch<T>               & mp,
                                const gsFunction<T>                 & thickness,
                                const std::vector<gsFunction<T> *>  & parameters,
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
                return new gsMaterialMatrixLinear<d,T>(mp,thickness,parameters);
    }
    else
    {
        if (impl==Implementation::Composite)
            GISMO_ERROR("Hyperelastic composites not available");
        //--------------------------------------NH Incompressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters);
        }
        //--------------------------------------NH Compressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
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
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        //
        //--------------------------------------MR Incompressible--------------------------------------------------
        else if ((mat==Material::MR) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters);
        }
        //---------------------------------------MR Compressible-------------------------------------------------
        else if      ((mat==Material::MR) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
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
            constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,thickness,parameters);
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
            constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,thickness,parameters);
        }
        else
            return NULL;
    }
}



} // namespace

