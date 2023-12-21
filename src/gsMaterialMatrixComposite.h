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

    enum MatIntegration isMatIntegrated() const {return MatIntegration::Integrated; }
    enum MatIntegration isVecIntegrated() const {return MatIntegration::Integrated; }

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    // template COM
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // template COM
    void pstretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}
    void pstretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}
    void pstress_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}
    void pstressDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void parameters_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    /// See \ref gsMaterialMatrixBase for details
    void transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    /// See \ref gsMaterialMatrixBase for details
    void covtransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    /// See \ref gsMaterialMatrixBase for details
    void pstressTransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

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
