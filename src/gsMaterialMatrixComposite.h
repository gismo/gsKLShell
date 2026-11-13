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
 * @brief      This class defines a linear material laminate
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

    /** @brief Constructor of the assembler object.

        \param[in] ...
        \param[in] ...

    */
    /// Default empty constructor
    gsMaterialMatrixComposite() { }

    gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const std::vector< gsFunctionSet<T> *>              & thickness,
                            const std::vector< gsFunctionSet<T> *>              & G,
                            const std::vector< gsFunctionSet<T> *>              & alpha,
                            const std::vector< gsFunctionSet<T> *>              & rho           );

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

    /// Thickness-integrated strain energy per unit reference (mid-surface) area.
    ///
    /// For each in-plane point k,
    /// \f[ W_k = \tfrac12 \varepsilon^T A \varepsilon + \varepsilon^T B \kappa + \tfrac12 \kappa^T D \kappa . \f]
    /// This is an exception to the normalisation of the Base class (energy per unit
    /// reference volume at fibre z): the law reports isMatIntegrated() = MatIntegration::Integrated,
    /// i.e. it works with laminate resultants rather than fibre-point stresses.
    ///
    /// A, B and D are the matrices returned by eval3D_matrix (MaterialOutput::MatrixA, MatrixB, MatrixD).
    /// \f$\varepsilon\f$ is the membrane strain and \f$\kappa\f$ the curvature change, in engineering
    /// Voigt form [E11, E22, E12+E21] and [K11, K22, K12+K21], formed as in eval3D_vector
    /// for VectorN and VectorM.
    ///
    /// W vanishes when the deformed and the undeformed geometry coincide.
    /// For B = 0 and symmetric A and D, dW/d\f$\varepsilon\f$ and dW/d\f$\kappa\f$ equal the
    /// VectorN and VectorM resultants of eval3D_vector.
    ///
    /// \param patch the patch index
    /// \param u     in-plane parametric points (one per column)
    /// \param z     only z.rows() is used, it sets the output width; the values are ignored
    /// \param out   ignored
    /// \returns a 1 x (z.rows()*u.cols()) matrix; column j*u.cols()+k holds W of point k,
    ///          identical for every fibre row j
    gsMatrix<T> eval3D_energy(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_energy;

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
    // NOTE: it could be that these matrices should be transposed!!
    gsMatrix<T> _cart2cov(const gsVector<T> & a1, const gsVector<T> & a2, const gsVector<T> & e1, const gsVector<T> & e2) const;
    gsMatrix<T> _con2cart(const gsVector<T> & ac1, const gsVector<T> & ac2, const gsVector<T> & e1, const gsVector<T> & e2) const;

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
