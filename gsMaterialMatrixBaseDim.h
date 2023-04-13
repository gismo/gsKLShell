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

    // enum {Linear=0}; // If true (1), this property entails that S = C *˙E

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
        GISMO_ASSERT(mp->targetDim()==dim,"Geometric dimension and the template dimension are not the same!");
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
    virtual void covtransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    virtual void contransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

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
    inline virtual index_t numParameters() const override { return m_pars.size(); }

    /// See \ref gsMaterialMatrixBase for details
    inline virtual void resetParameters()
    {
        m_pars.clear();
        m_pars.resize(0);
    }


public:

    void _computePoints(const index_t patch, const gsMatrix<T> & u, bool basis = true) const;

    /// Sets the Tension Field Theory parameter theta for the evaluation points to be used (note; can this be done better?)
    void setTheta(gsVector<T> & thetas) { m_thetas = thetas; }

    /// Sets the Tension Field Theory parameter gamma for the evaluation points to be used (note; can this be done better?)
    void setGamma(gsVector<T> & gammas) { m_gammas = gammas; }

    /// Returns the covariant a tensor on the deformed geometry
    gsMatrix<T,2,2> _getAcov_def(index_t k, T z) const;

    /// Returns the contravariant a tensor on the deformed geometry
    gsMatrix<T,2,2> _getAcon_def(index_t k, T z) const;

    /// Returns the covariant b tensor on the deformed geometry
    gsMatrix<T,2,2> _getBcov_def(index_t k, T z) const;

    /// Returns the covariant n tensor on the deformed geometry
    gsMatrix<T,3,2> _getncov_def(index_t k, T z) const;

    /// Returns the covariant metric tensor on the deformed geometry
    gsMatrix<T,3,3> _getGcov_def(index_t k, T z) const;

    /// Returns the contravariant metric tensor on the deformed geometry
    gsMatrix<T,3,3> _getGcon_def(index_t k, T z) const;

    /// Returns the covariant basis vector a on the deformed geometry
    gsMatrix<T,3,2> _getacov_def(index_t k, T z) const;

    /// Returns the contravariant basis vector a on the deformed geometry
    gsMatrix<T,3,2> _getacon_def(index_t k, T z) const;

    /// Returns the covariant basis vector g on the deformed geometry
    gsMatrix<T,3,3> _getgcov_def(index_t k, T z) const;

    /// Returns the contravariant basis vector g on the deformed geometry
    gsMatrix<T,3,3> _getgcon_def(index_t k, T z) const;

    /// Returns the covariant a tensor on the original geometry
    gsMatrix<T,2,2> _getAcov_ori(index_t k, T z) const;

    /// Returns the contravariant a tensor on the original geometry
    gsMatrix<T,2,2> _getAcon_ori(index_t k, T z) const;

    /// Returns the covariant b tensor on the original geometry
    gsMatrix<T,2,2> _getBcov_ori(index_t k, T z) const;

    /// Returns the covariant n tensor on the original geometry
    gsMatrix<T,3,2> _getncov_ori(index_t k, T z) const;

    /// Returns the covariant metric tensor on the original geometry
    gsMatrix<T,3,3> _getGcov_ori(index_t k, T z) const;

    /// Returns the contravariant metric tensor on the original geometry
    gsMatrix<T,3,3> _getGcon_ori(index_t k, T z) const;

    /// Returns the covariant basis vector a on the original geometry
    gsMatrix<T,3,2> _getacov_ori(index_t k, T z) const;

    /// Returns the contravariant basis vector a on the original geometry
    gsMatrix<T,3,2> _getacon_ori(index_t k, T z) const;

    /// Returns the covariant basis vector g on the original geometry
    gsMatrix<T,3,3> _getgcov_ori(index_t k, T z) const;

    /// Returns the contravariant basis vector g on the original geometry
    gsMatrix<T,3,3> _getgcon_ori(index_t k, T z) const;

    /// Computes metric quantities on the deformed geometry
    void _computeMetricDeformed(const index_t patch, const gsMatrix<T> & u, bool basis = true) const;

    /// Computes metric quantities on the undeformed geometry
    void _computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u, bool basis = true) const;

    /// Gets metric quantities on the deformed and undeformed geometries
    void _getMetric(index_t k, T z, bool basis = true) const;

    /// Gets metric quantities on the deformed and undeformed geometries
    void _getMetric(index_t k, T z, const gsMatrix<T> & C, bool basis = true) const;

    /// Gets metric quantities on the deformed geometry
    void _getMetricDeformed(index_t k, T z, bool basis = true) const;

    void _getMetricDeformed(const gsMatrix<T> & C) const;

    /// Gets metric quantities on the undeformed geometry
    void _getMetricUndeformed(index_t k, T z, bool basis = true) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalStretch(const gsMatrix<T> & C, const gsMatrix<T> & gcon_ori ) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalPStress(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalPStrain(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computeStretch(const gsMatrix<T> & C, const gsMatrix<T> & gcon_ori ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computePStress(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computePStrain(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    gsMatrix<T> _transformation(const gsMatrix<T> & basis1, const gsMatrix<T> & basis2 ) const;

    // void setDeformed(const gsFunctionSet<T> & deformed) {this->setDeformed(deformed);}

    bool hasDeformed() const { return m_defpatches->nPieces()!=0; }

private:
    /// Implementation of \ref _computeMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricDeformed_impl(const index_t patch, bool basis) const;

    /// Implementation of \ref _computeMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricDeformed_impl(const index_t patch, bool basis) const;

    /// Implementation of \ref _getMetric for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricUndeformed_impl(const index_t patch, bool basis) const;

    /// Implementation of \ref _getMetric for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricUndeformed_impl(const index_t patch, bool basis) const;

    /// Implementation of \ref _getMetricDeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetric_impl(index_t k, T z, bool basis) const;

    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetric_impl(index_t k, T z, bool basis) const;

    /// Implementation of \ref _getMetricDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricDeformed_impl(index_t k, T z, bool basis) const;

    /// Implementation of \ref _getMetricDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricDeformed_impl(index_t k, T z, bool basis) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricUndeformed_impl(index_t k, T z, bool basis) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricUndeformed_impl(index_t k, T z, bool basis) const;

    // -------------------------------------------------------------------------------------------------------------
    // Separate getters
    // -------------------------------------------------------------------------------------------------------------

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,2,2>>::type _getBcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,2,2>>::type _getBcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,2>>::type _getncov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,2>>::type _getncov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type _getGcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type _getGcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type _getGcon_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type _getGcon_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type _getgcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type _getgcov_def_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,2,2>>::type _getBcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,2,2>>::type _getBcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,2>>::type _getncov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,2>>::type _getncov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type _getGcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type _getGcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type _getGcon_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type _getGcon_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type _getgcov_ori_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type _getgcov_ori_impl(index_t k, T z) const;


protected:

    void membersSetZero()
    {
        m_Acov_ori.setZero(); m_Acon_ori.setZero(); m_Acov_def.setZero(); m_Acon_def.setZero(); m_Bcov_ori.setZero(); m_Bcon_ori.setZero(); m_Bcov_def.setZero(); m_Bcon_def.setZero();
        m_acov_ori.setZero(); m_acon_ori.setZero(); m_acov_def.setZero(); m_acon_def.setZero();
        m_ncov_ori.setZero(); m_ncov_def.setZero();
        m_Gcov_ori.setZero(); m_Gcon_ori.setZero(); m_Gcov_def.setZero(); m_Gcon_def.setZero(); m_Gcov_ori_L.setZero(); m_Gcov_def_L.setZero();
        m_gcov_ori.setZero(); m_gcov_def.setZero();m_gcon_ori.setZero(); m_gcon_def.setZero();
        m_Acov_ori_mat.setZero(); m_Acon_ori_mat.setZero(); m_Acov_def_mat.setZero(); m_Acon_def_mat.setZero(); m_Bcov_ori_mat.setZero(); m_Bcov_def_mat.setZero();
        m_acov_ori_mat.setZero(); m_acon_ori_mat.setZero(); m_acov_def_mat.setZero(); m_acon_def_mat.setZero(); m_ncov_ori_mat.setZero(); m_ncov_def_mat.setZero();
    }

    mutable gsVector<T> m_thetas, m_gammas;

    const gsFunctionSet<T> * m_patches;
    using Base::m_defpatches;

    std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_thickness;
    const gsFunction<T> * m_density;
    // Material parameters and kinematics
    mutable gsMatrix<T> m_parmat;
    mutable gsVector<T> m_parvals;
    mutable gsMatrix<T> m_Tmat,m_rhomat;

    // Geometric data point
    mutable util::gsThreaded<gsMapData<T> > m_map_ori, m_map_def;

    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T,dim,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def;
    mutable gsMatrix<T,3,2> m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T,3,3> m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T,3,3> m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;

    mutable gsMatrix<T> m_stretches, m_stretchvec, m_pstress, m_pstressvec, m_pstrain, m_pstrainvec;

    mutable T           m_J0_sq, m_J_sq;
// public:
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW //must be present whenever the class contains fixed size matrices

};

#ifdef GISMO_BUILD_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsMaterialMatrixBaseDim
   */
  void pybind11_init_gsMaterialMatrixBaseDim2(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixBaseDim3(pybind11::module &m);

#endif // GISMO_BUILD_PYBIND11

} // namespace

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixBaseDim.hpp)
#endif
