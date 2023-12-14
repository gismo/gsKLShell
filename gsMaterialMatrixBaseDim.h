/** @file gsMaterialMatrixBaseDim.h

    @brief Base class with dimension in template; used for metric computations

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

public:

    typedef T Scalar_t;

    using Base = gsMaterialMatrixBase<T>;

    typedef typename Base::function_ptr function_ptr;

    // enum {Linear=0}; // If true (1), this property entails that S = C *˙E

    gsMaterialMatrixBaseDim() 
    : 
    Base(nullptr,nullptr,nullptr,nullptr)
    {
        this->defaultOptions();
        membersSetZero();
    }

    gsMaterialMatrixBaseDim(const gsFunctionSet<T> * mp)
    :
    Base(mp,nullptr,nullptr,nullptr)
    {
        GISMO_ASSERT(mp->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        this->defaultOptions();
        membersSetZero();
    }

    gsMaterialMatrixBaseDim(const gsFunctionSet<T> * mp,
                            const gsFunctionSet<T> * mp_def)
    :
    Base(mp,mp_def,nullptr,nullptr)
    {
        GISMO_ASSERT(mp->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        GISMO_ASSERT(mp_def->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        this->defaultOptions();
        membersSetZero();
    }

    gsMaterialMatrixBaseDim(const gsFunctionSet<T> * mp,
                            const gsFunctionSet<T> * thickness,
                            const gsFunctionSet<T> * Density)
    :
    Base(mp,nullptr,thickness,Density)
    {
        GISMO_ASSERT(mp->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        this->defaultOptions();
        membersSetZero();
    }

    gsMaterialMatrixBaseDim(const gsFunctionSet<T> * mp,
                            const gsFunctionSet<T> * mp_def,
                            const gsFunctionSet<T> * thickness,
                            const gsFunctionSet<T> * Density)
    :
    Base(mp,mp_def,thickness,Density)
    {
        GISMO_ASSERT(mp->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        GISMO_ASSERT(mp_def->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        this->defaultOptions();
        membersSetZero();
    }

    /// Destructor
    virtual ~gsMaterialMatrixBaseDim() {}

public:
    /// See \ref gsMaterialMatrixBase for details
    virtual void defaultOptions() override;

    /// See \ref gsMaterialMatrixBase for details
    virtual void    density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void  thickness_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void parameters_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_spec2cov(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_spec2con(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_cov2cart(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_con2cart(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_deformation(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const;
    
    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_pstretch(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_pstretchDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual gsMatrix<T> eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    virtual bool initialized() const override
    {
        return ((m_thickness != nullptr)); // && (m_patches!= nullptr) //HV: could be defined without geometry, then call "setUndeformed"
    }

public:

    void _computePoints(const index_t patch, const gsMatrix<T> & u) const;

    /// Returns the covariant a tensor on the deformed geometry
    gsMatrix<T> _getAcov_def(index_t k, T z) const;

    /// Returns the contravariant a tensor on the deformed geometry
    gsMatrix<T> _getAcon_def(index_t k, T z) const;

    /// Returns the covariant b tensor on the deformed geometry
    gsMatrix<T> _getBcov_def(index_t k, T z) const;

    /// Returns the covariant n tensor on the deformed geometry
    gsMatrix<T> _getncov_def(index_t k, T z) const;

    /// Returns the covariant metric tensor on the deformed geometry
    gsMatrix<T> _getGcov_def(index_t k, T z) const;

    /// Returns the contravariant metric tensor on the deformed geometry
    gsMatrix<T> _getGcon_def(index_t k, T z) const;

    /// Returns the covariant basis vector a on the deformed geometry
    gsMatrix<T> _getacov_def(index_t k, T z) const;

    /// Returns the contravariant basis vector a on the deformed geometry
    gsMatrix<T> _getacon_def(index_t k, T z) const;

    /// Returns the covariant basis vector g on the deformed geometry
    gsMatrix<T> _getgcov_def(index_t k, T z) const;

    /// Returns the contravariant basis vector g on the deformed geometry
    gsMatrix<T> _getgcon_def(index_t k, T z) const;

    /// Returns the covariant a tensor on the original geometry
    gsMatrix<T> _getAcov_ori(index_t k, T z) const;

    /// Returns the contravariant a tensor on the original geometry
    gsMatrix<T> _getAcon_ori(index_t k, T z) const;

    /// Returns the covariant b tensor on the original geometry
    gsMatrix<T> _getBcov_ori(index_t k, T z) const;

    /// Returns the covariant n tensor on the original geometry
    gsMatrix<T> _getncov_ori(index_t k, T z) const;

    /// Returns the covariant metric tensor on the original geometry
    gsMatrix<T> _getGcov_ori(index_t k, T z) const;

    /// Returns the contravariant metric tensor on the original geometry
    gsMatrix<T> _getGcon_ori(index_t k, T z) const;

    /// Returns the covariant basis vector a on the original geometry
    gsMatrix<T> _getacov_ori(index_t k, T z) const;

    /// Returns the contravariant basis vector a on the original geometry
    gsMatrix<T> _getacon_ori(index_t k, T z) const;

    /// Returns the covariant basis vector g on the original geometry
    gsMatrix<T> _getgcov_ori(index_t k, T z) const;

    /// Returns the contravariant basis vector g on the original geometry
    gsMatrix<T> _getgcon_ori(index_t k, T z) const;

    /// Computes metric quantities on the deformed geometry
    void _computeMetricDeformed(const index_t patch, const gsMatrix<T> & u) const;

    /// Computes metric quantities on the undeformed geometry
    void _computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u) const;

    /// Gets metric quantities on the deformed and undeformed geometries
    void _getMetric(const index_t k, const T z) const;

    /// Gets metric quantities on the deformed and undeformed geometries
    void _getMetric(index_t k, T z, const gsMatrix<T> & C) const;

    /// Gets metric quantities on the deformed geometry
    void _getMetricDeformed(const index_t k, const T z) const;

    void _getMetricDeformed(const gsMatrix<T> & C) const;

    /// Gets metric quantities on the undeformed geometry
    void _getMetricUndeformed(const index_t k, const T z) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalStretch(const gsMatrix<T> & C, const gsMatrix<T> & gcon_ori ) const;

    /// Computes the principal strain given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalPStrain(const gsMatrix<T> & C ) const;

    /// Computes the principal stress given stress tensor S, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalPStress(const gsMatrix<T> & S ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computeStretch(const gsMatrix<T> & C, const gsMatrix<T> & gcon_ori ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computePStrain(const gsMatrix<T> & C ) const;

    /// Computes the principal stresses of a given stress tensor S, into class members m_pstress and m_pstressvec
    void _computePStress(const gsMatrix<T> & C ) const;


    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    gsMatrix<T> _transformation(const gsMatrix<T> & basis1, const gsMatrix<T> & basis2 ) const;

    virtual void setUndeformed(const gsFunctionSet<T> * undeformed) override
    {
        GISMO_ASSERT(undeformed->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        function_ptr f_ptr = memory::make_shared_not_owned(undeformed);
        Base::setUndeformed(f_ptr);
    }

    virtual void setDeformed(const gsFunctionSet<T> * deformed) override
    {
        GISMO_ASSERT(deformed->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        function_ptr f_ptr = memory::make_shared_not_owned(deformed);
        Base::setDeformed(f_ptr);
    }

    virtual void setUndeformed(const function_ptr undeformed) override
    {
        GISMO_ASSERT(undeformed->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        Base::setUndeformed(undeformed);
    }

    virtual void setDeformed(const function_ptr deformed) override
    {
        GISMO_ASSERT(deformed->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
        Base::setDeformed(deformed);
    }

private:

    template <short_t _type>
    typename std::enable_if<_type==0, index_t>::type _tensionField(const gsVector<T> &, const gsVector<T> & Ep) const;

    template <short_t _type>
    typename std::enable_if<_type==1, index_t>::type _tensionField(const gsVector<T> &, const gsVector<T> & Ep) const;

    template <short_t _type>
    typename std::enable_if<_type==2, index_t>::type _tensionField(const gsVector<T> &, const gsVector<T> & Ep) const;

    template <short_t _type>
    typename std::enable_if<_type!=0 &&
                            _type!=1 &&
                            _type!=2, index_t>::type _tensionField(const gsVector<T> &, const gsVector<T> & Ep) const
    {GISMO_NO_IMPLEMENTATION;}

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

    // -------------------------------------------------------------------------------------------------------------
    // Separate getters
    // -------------------------------------------------------------------------------------------------------------

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getBcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getBcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getncov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getncov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getGcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getGcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getGcon_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getGcon_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getgcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getgcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getBcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getBcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getncov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getncov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getGcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getGcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getGcon_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getGcon_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T>>::type _getgcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T>>::type _getgcov_ori_impl(index_t k, T z) const;


protected:

    void membersSetZero()
    {
        m_data.mine().membersSetZero();
    }

    using Base::m_patches;
    using Base::m_defpatches;

    using Base::m_options;

    using Base::m_pars;
    using Base::m_thickness;
    using Base::m_density;

    // Geometric data point
    mutable util::gsThreaded< gsMaterialMatrixBaseDimData<dim,T> > m_data;
public:

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //must be present whenever the class contains fixed size matrices
#   undef Eigen
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

    mutable gsMatrix<T> m_stretches, m_stretchvec, m_pstress, m_pstressvec, m_pstrain, m_pstrainvec;

    mutable T           m_J0_sq, m_J_sq;

    mutable gsVector<T> m_thetas, m_gammas;

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
