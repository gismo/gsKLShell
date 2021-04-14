/** @file gsMaterialMatrixBase.h

    @brief Provides a base class for material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsCore/gsFuncData.h>

namespace gismo
{

/**
 * @brief      This class defines the base class for material matrices
 *
 * @tparam     T     Real type
 *
 * @ingroup    MaterialMatrix
 */
template <short_t dim, class T>
class gsMaterialMatrixBaseDim : public gsMaterialMatrixBase<T>
{
    using Base = gsMaterialMatrixBase<T>;

public:

    gsMaterialMatrixBaseDim() {};


    gsMaterialMatrixBaseDim(const gsFunctionSet<T> & mp)
    :
    m_patches(&mp)
    { };

    gsMaterialMatrixBaseDim(const gsFunctionSet<T> & mp,
                            const gsFunctionSet<T> & mp_def)
    :
    m_patches(&mp)
    {
        this->setDeformed(mp_def);
    };


    /// Destructor
    virtual ~gsMaterialMatrixBaseDim() {};

public:

    /// Computes metric quantities on the deformed geometry
    void _computeMetricDeformed(const index_t patch, const gsMatrix<T> & u) const;

    /// Computes metric quantities on the undeformed geometry
    void _computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u) const;

    /// Gets metric quantities on the deformed and undeformed geometries
    void _getMetric(index_t k, T z) const;

    /// Gets metric quantities on the deformed geometry
    void _getMetricDeformed(index_t k, T z) const;

    /// Gets metric quantities on the undeformed geometry
    void _getMetricUndeformed(index_t k, T z) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computeStretch(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalStretch(const gsMatrix<T> & C ) const;

private:
    /// Implementation of \ref _computeMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricDeformed_impl(const index_t patch, const gsMatrix<T> & u) const;

    /// Implementation of \ref _computeMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricDeformed_impl(const index_t patch, const gsMatrix<T> & u) const;

    /// Implementation of \ref _getMetric for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricUndeformed_impl(const index_t patch, const gsMatrix<T> & u) const;

    /// Implementation of \ref _getMetric for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricUndeformed_impl(const index_t patch, const gsMatrix<T> & u) const;

    /// Implementation of \ref _getMetricDeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetric_impl(index_t k, T z) const;

    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetric_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricDeformed_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricDeformed_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricUndeformed_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricUndeformed_impl(index_t k, T z) const;


protected:
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;

    // Geometric data point
    mutable gsMapData<T> m_map, m_map_def;

    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T,dim,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def;
    mutable gsMatrix<T,3,2> m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T,3,3> m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T,3,3> m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;

    mutable gsMatrix<T> m_stretches, m_stretchvec;

    mutable T           m_J0_sq, m_J_sq;


};

} // namespace

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixBaseDim.hpp)
#endif