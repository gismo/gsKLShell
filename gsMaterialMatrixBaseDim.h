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
#include <gsUtils/gsThreaded.h>

namespace gismo
{

template <short_t dim, class T>
class gsMaterialMatrixBaseDimData;

/**
 * @brief      This class defines the base class for material matrices
 *
 * @tparam     T     Real type
 *
 * @ingroup    KLShell
 */
template <short_t dim, class T>
class gsMaterialMatrixBaseDim : public gsMaterialMatrixBase<T>
{
    using Base = gsMaterialMatrixBase<T>;

public:

    gsMaterialMatrixBaseDim() : m_patches(nullptr), m_thickness(nullptr), m_density(nullptr)
    {
        membersSetZero();
    }

    gsMaterialMatrixBaseDim(const gsFunctionSet<T> & mp)
    :
    m_patches(&mp), m_thickness(nullptr), m_density(nullptr)
    {
        GISMO_ASSERT(mp.targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        membersSetZero();
    }

    gsMaterialMatrixBaseDim(    const gsFunctionSet<T> * mp,
                                const gsFunctionSet<T> * mp_def,
                                const gsFunction<T> * thickness,
                                const gsFunction<T> * Density)
    :
    m_patches(mp),
    m_thickness(thickness),
    m_density(Density)
    {
        membersSetZero();
        this->setDeformed(mp_def);
    }

    /// Destructor
    virtual ~gsMaterialMatrixBaseDim() {}

public:
    /// See \ref gsMaterialMatrixBase for details
    virtual void    density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void  thickness_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void parameters_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void deformation_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    virtual void setDensity(const gsFunction<T> & Density)
    {
        m_density = const_cast<gsFunction<T> *>(&Density);
    }
    /// Gets the Density
    virtual gsFunction<T> * getDensity() {return const_cast<gsFunction<T> *>(m_density);}

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    inline virtual void setParameters(const std::vector<gsFunction<T>*> &pars)
    {
        m_pars = pars;
    }
    /**
     * @brief      Gets the number of parameters
     *
     */
    inline virtual index_t numParameters() { return m_pars.size(); }
    /// See \ref gsMaterialMatrixBase for details
    inline virtual void resetParameters()
    {
        m_pars.clear();
        m_pars.resize(0);
    }

    inline virtual bool initialized() const override
    {
        return m_thickness != nullptr && m_patches!= nullptr;
    }

public:

    void _computePoints(const index_t patch, const gsMatrix<T> & u) const;

    /// Computes metric quantities on the deformed geometry
    void _computeMetricDeformed(const index_t patch, const gsMatrix<T> & u) const;

    /// Computes metric quantities on the undeformed geometry
    void _computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u) const;

    /// Gets metric quantities on the deformed and undeformed geometries
    void _getMetric(const index_t k, const T z) const;

    /// Gets metric quantities on the deformed geometry
    void _getMetricDeformed(const index_t k, const T z) const;

    /// Gets metric quantities on the undeformed geometry
    void _getMetricUndeformed(const index_t k, const T z) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalStretch(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computeStretch(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    gsMatrix<T> _transformation(const gsMatrix<T> & basis1, const gsMatrix<T> & basis2 ) const;

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
    typename std::enable_if<_dim==2, void>::type _getMetric_impl(const index_t k, const T z) const;

    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetric_impl(const index_t k, const T z) const;

    /// Implementation of \ref _getMetricDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricDeformed_impl(const index_t k, const T z) const;

    /// Implementation of \ref _getMetricDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricDeformed_impl(const index_t k, const T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricUndeformed_impl(const index_t k, const T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricUndeformed_impl(const index_t k, const T z) const;

protected:

    void membersSetZero()
    {
        m_data.mine().membersSetZero();
    }


    const gsFunctionSet<T> * m_patches;
    using Base::m_defpatches;

    std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_thickness;
    const gsFunction<T> * m_density;

    // Geometric data point
    mutable util::gsThreaded< gsMaterialMatrixBaseDimData<dim,T> > m_data;
public:
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen //must be present whenever the class contains fixed size matrices

};

template<short_t dim, class T>
class gsMaterialMatrixBaseDimData
{
public:

    void membersSetZero()
    {
        m_Acov_ori.setZero(); m_Acon_ori.setZero(); m_Acov_def.setZero(); m_Acon_def.setZero(); m_Bcov_ori.setZero(); m_Bcon_ori.setZero(); m_Bcov_def.setZero(); m_Bcon_def.setZero();
        m_acov_ori.setZero(); m_acon_ori.setZero(); m_acov_def.setZero(); m_acon_def.setZero();
        m_ncov_ori.setZero(); m_ncov_def.setZero();
        m_Gcov_ori.setZero(); m_Gcon_ori.setZero(); m_Gcov_def.setZero(); m_Gcon_def.setZero(); m_Gcov_ori_L.setZero(); m_Gcov_def_L.setZero();
        m_gcov_ori.setZero(); m_gcov_def.setZero();m_gcon_ori.setZero(); m_gcon_def.setZero();
        m_Acov_ori_mat.setZero(); m_Acon_ori_mat.setZero(); m_Acov_def_mat.setZero(); m_Acon_def_mat.setZero(); m_Bcov_ori_mat.setZero(); m_Bcov_def_mat.setZero();
        m_acov_ori_mat.setZero(); m_acon_ori_mat.setZero(); m_acov_def_mat.setZero(); m_acon_def_mat.setZero(); m_ncov_ori_mat.setZero(); m_ncov_def_mat.setZero(); m_normal_ori_mat.setZero(); m_normal_def_mat.setZero();

    }
    // Material parameters and kinematics
    mutable gsMatrix<T> m_parmat;
    mutable gsVector<T> m_parvals;
    mutable gsMatrix<T> m_Tmat,m_rhomat;

    mutable gsMatrix<T> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def;
    mutable gsMatrix<T> m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T> m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T> m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat, m_normal_ori_mat, m_normal_def_mat;

    mutable gsMatrix<T> m_stretches, m_stretchvec;

    mutable T           m_J0_sq, m_J_sq;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen //must be present whenever the class contains fixed size matrices
};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsMaterialMatrixBaseDim
   */
  void pybind11_init_gsMaterialMatrixBaseDim2(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixBaseDim3(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixBaseDim.hpp)
#endif
