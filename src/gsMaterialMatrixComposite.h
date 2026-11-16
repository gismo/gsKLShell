/** @file gsMaterialMatrixComposite.h

    @brief Provides a material matrix for laminates

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
#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

/**
 * @brief      This class defines a linear material laminate (classical laminate theory)
 *
 * The laminate consists of plies stacked along the unit normal n = a1 x a2 / |a1 x a2|
 * (for dim == 2, a planar mid-surface, n is the out-of-plane unit vector). Ply 0 is
 * the ply at the most negative z along n, the last ply is at the most positive z.
 *
 * Per ply, the following functions are given and evaluated at physical points:
 *  - thickness: scalar;
 *  - G: the 3x3 Voigt [11,22,12] stiffness in the fibre axes with engineering shear, passed as a
 *    9x1 column. For example, the output of gsCompositeMatrix(Exx,Eyy,Gxy,nuxy,nuyx) reshaped
 *    with G.resize(9,1);
 *  - alpha: the fibre angle in radians, measured in the tangent plane from e1 = a1/|a1| towards
 *    e2 = a^2/|a^2|;
 *  - rho: the density (optional).
 *
 * The matrices A, B and D are returned through-thickness integrated (MatIntegration::Integrated)
 * in contravariant components, with A = sum_i t_i C_i, B = sum_i t_i z_i C_i and
 * D = sum_i (t_i z_i^2 + t_i^3/12) C_i, where t_i is the thickness of ply i, C_i its rotated
 * stiffness in contravariant components and z_i the signed height of its mid-surface above the
 * laminate mid-plane. The vectors are N = A eps + B kappa and M = B^T eps + D kappa,
 * with eps the membrane strain and kappa = b_ori - b_def the curvature change.
 *
 * @tparam     dim   The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T     Real type
 *
 * @ingroup    KLShell
 *
 */
template <  short_t dim,
            class T
         >
class gsMaterialMatrixComposite : public gsMaterialMatrixBaseDim<dim,T>
{
public:

    GISMO_CLONE_FUNCTION(gsMaterialMatrixComposite)

    using Base = gsMaterialMatrixBaseDim<dim,T>;

    typedef typename Base::function_ptr function_ptr;

    /// Default empty constructor
    gsMaterialMatrixComposite() { }

    /** @brief Constructor of the laminate with density.

        @param[in] mp         The undeformed geometry
        @param[in] thickness  Per-ply scalar thickness, ply 0 at the most negative z
        @param[in] G          Per-ply 3x3 Voigt [11,22,12] stiffness in fibre axes (engineering
                              shear), as a 9x1 column
        @param[in] alpha      Per-ply fibre angle in radians, from e1 = a1/|a1| towards e2 = a^2/|a^2|
        @param[in] rho        Per-ply density
    */
    gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const std::vector< gsFunctionSet<T> *>              & thickness,
                            const std::vector< gsFunctionSet<T> *>              & G,
                            const std::vector< gsFunctionSet<T> *>              & alpha,
                            const std::vector< gsFunctionSet<T> *>              & rho           );

    /** @brief Constructor of the laminate without density.

        @param[in] mp         The undeformed geometry
        @param[in] thickness  Per-ply scalar thickness, ply 0 at the most negative z
        @param[in] G          Per-ply 3x3 Voigt [11,22,12] stiffness in fibre axes (engineering
                              shear), as a 9x1 column
        @param[in] alpha      Per-ply fibre angle in radians, from e1 = a1/|a1| towards e2 = a^2/|a^2|
    */
    gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const std::vector< gsFunctionSet<T> *>              & thickness,
                            const std::vector< gsFunctionSet<T> *>              & G,
                            const std::vector< gsFunctionSet<T> *>              & alpha         );

    enum MatIntegration isMatIntegrated() const override {return MatIntegration::Integrated; }
    enum MatIntegration isVecIntegrated() const override {return MatIntegration::Integrated; }

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() override {return m_options;}
    void setOptions(gsOptionList opt) override {m_options.update(opt,gsOptionList::addIfUnknown); }

    // template COM
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const override;
    // template COM
    void pstretch_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const override
    {GISMO_NO_IMPLEMENTATION;}
    void pstretchDir_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const override
    {GISMO_NO_IMPLEMENTATION;}
    void pstress_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const override
    {GISMO_NO_IMPLEMENTATION;}
    void pstressDir_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const override
    {GISMO_NO_IMPLEMENTATION;}

    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const override;

    /// See \ref gsMaterialMatrixBase for details
    void parameters_into(const index_t /*patch*/, const gsMatrix<T> & /*u*/, gsMatrix<T>& /*result*/) const override
    {GISMO_NO_IMPLEMENTATION;}

    /// See \ref gsMaterialMatrixBase for details
    void transform_into(const index_t /*patch*/, const gsMatrix<T> & /*u*/, gsMatrix<T>& /*result*/) const override
    {GISMO_NO_IMPLEMENTATION;}

    /// See \ref gsMaterialMatrixBase for details
    void covtransform_into(const index_t /*patch*/, const gsMatrix<T> & /*u*/, gsMatrix<T>& /*result*/) const
    {GISMO_NO_IMPLEMENTATION;}

    /// See \ref gsMaterialMatrixBase for details
    void pstressTransform_into(const index_t /*patch*/, const gsMatrix<T> & /*u*/, gsMatrix<T>& /*result*/) const
    {GISMO_NO_IMPLEMENTATION;}

    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_matrix;
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_vector;

    std::ostream &print(std::ostream &os) const override;

    bool initialized() const override
    {
        return m_patches!=nullptr;
    }

public:
    /// Shared pointer for gsMaterialMatrixComposite
    typedef memory::shared_ptr< gsMaterialMatrixComposite > Ptr;

    /// Unique pointer for gsMaterialMatrixComposite
    typedef memory::unique_ptr< gsMaterialMatrixComposite > uPtr;

protected:
    void _initialize(const index_t nLayers);
    void _defaultOptions();

    gsMatrix<T> _transformationMatrix(const gsMatrix<T> & phi, const gsMatrix<T> & u) const;

protected:

    // template MAT
    void _computePoints(const index_t patch, const gsMatrix<T> & u, bool basis = true) const;

protected:

    // constructor
    using Base::m_patches;
    using Base::m_defpatches;

    // Composite
    index_t m_nLayers;

    // Geometric data
    using Base::m_data;


    gsOptionList m_options;

    std::vector< function_ptr > m_Ts;
    std::vector< function_ptr > m_Gs;
    std::vector< function_ptr > m_As;
    std::vector< function_ptr > m_Rs;
    mutable util::gsThreaded<std::vector< gsMatrix<T> >> m_Gcontainer, m_Tcontainer, m_Acontainer, m_Rcontainer;
};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixComposite.hpp)
#endif
