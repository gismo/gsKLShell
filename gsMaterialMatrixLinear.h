/** @file gsMaterialMatrixLinear.h

    @brief Provides material matrices for the thin shell class.

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
#include <gsIO/gsOptionList.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{


/**
 * @brief      This class defines a linear material
 *
 * @tparam     dim   The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T     Real type
 *
 * @ingroup    MaterialMatrix
 *
 */
template <  short_t dim,
            class T
         >
class gsMaterialMatrixLinear : public gsMaterialMatrixBase<T>
{
public:
    /**
     * @brief      Constructor without deformed multipatch and density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    /**
     * @brief      Constructor without density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  mp_def     Deformed geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    /**
     * @brief      Full constructor
     *
     * @param[in]  mp         Original geometry
     * @param[in]  mp_def     Deformed geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars,
                        const gsFunction<T> & Density);

    /// Destructor
    gsMaterialMatrixLinear() { }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isMatIntegrated() const {return MatIntegration::Constant; }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isVecIntegrated() const {return MatIntegration::Constant; }

    /// See \ref gsMaterialMatrixBase for details
    gsOptionList & options() {return m_options;}

    /// See \ref gsMaterialMatrixBase for details
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt

    /// See \ref gsMaterialMatrixBase for details
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const
    { GISMO_NO_IMPLEMENTATION; }

    /// See \ref gsMaterialMatrixBase for details
    void setParameters(const std::vector<gsFunction<T>*> &pars)
    {
        m_pars = pars;
        m_numPars = m_pars.size();
    }

    /// See \ref gsMaterialMatrixBase for details
    void info() const;

public:
    /// Shared pointer for gsMaterialMatrixLinear
    typedef memory::shared_ptr< gsMaterialMatrixLinear > Ptr;

    /// Unique pointer for gsMaterialMatrixLinear
    typedef memory::unique_ptr< gsMaterialMatrixLinear > uPtr;

protected:
    /**
     * @brief      Initializes the object.
     *
     * Initializes options, flags and defines the number of parameters
     *
     */
    void _initialize();

    /**
     * @brief      Sets default options
     */
    void _defaultOptions();

protected:
    /**
     * @brief      Computes the linear material matrix entry with indices \a i \a j \a k \a l
     *
     * The entry is computed by \f$ \mathcal{C}^{ijkl} = \frac{2\lambda\mu}{\lambda+2\mu}a^{ij}a^{kl} + \mu (a^{ik}a^{jl} + a^{il}a^{jk})\f$ where \f$\lambda\f$ and \f$\mu\f$
     * are the Lamé parameters and \f$a^{ij}=\mathbf{a}^i\cdot \mathbf{a}^j$ with \f$\mathbf{a}^i\f$ is the contravariant vector \f$ i \f$
     *
     * @param[in]  i,j,k,l  Indices
     *
     * @return     Cijkl
     */
    T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;

    /**
     * @brief      Computes the linear force/moment entry with indices \a i \a j \a k \a l at height z
     *
     * Computes the thickness-integrated stress tensor, i.e. the normal force (0th thickness-moment) or the bending moment (1st thickness-moment).
     * Sij is computed as \f$ \mathcal{C}^{ijkl} : \mathbf{E}_{ij} \f$ where \f$\mathbf{E}_{ij} = a_{ij} z b_{ij}\f$ with \f$ a_{ij}\f$ the in-plane metric and \f$b_{ij}\f$ the curvature.
     * According to this definition, MaterialOutput out==VectorN would return the 0th moment, hence the \f$a_{ij}\f$-part would be relevant and for VectorM the \f$b_{ij}$ part is relevant.
     * The term could also be integrated (out==Generalized).
     *
     * @param[in]  i,j   Indices
     * @param[in]  z     Through-thickness coordinate
     * @param[in]  out   Output specification
     *
     * @return     Sij
     */
    T _Sij    (const index_t i, const index_t j, const T z, enum MaterialOutput out) const;

    /**
     * @brief      Computes the map, the metric quantities and the parameters on
     *             specified points.
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane point coordinates
     */
    void _computePoints(const index_t patch, const gsMatrix<T> & u) const;

    /// Computes metric quantities on the deformed geometry
    void _computeMetricDeformed() const;

    /// Computes metric quantities on the undeformed geometry
    void _computeMetricUndeformed() const;

    /// Gets metric quantities on the deformed and undeformed geometries
    void _getMetric(index_t k, T z) const;

    /// Gets metric quantities on the deformed geometry
    void _getMetricUndeformed(index_t k, T z) const;

    /// Gets metric quantities on the undeformed geometry
    void _getMetricDeformed(index_t k, T z) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computeStretch(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalStretch(const gsMatrix<T> & C ) const;

private:

    /// Implementation of \ref _computeMetricDeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricDeformed_impl() const;

    /// Implementation of \ref _computeMetricDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricDeformed_impl() const;

    /// Implementation of \ref _computeMetricUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricUndeformed_impl() const;

    /// Implementation of \ref _computeMetricUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricUndeformed_impl() const;

    /// Implementation of \ref _getMetric for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetric_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetric for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetric_impl(index_t k, T z) const;

    /// Implementation of \ref _getMetricDeformed for planar geometries
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
    // general
    index_t m_numPars; // how many parameters for the material model?

    // constructor
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;
    const gsFunction<T> * m_thickness;
    std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_density;

    // Linear material matrix
    mutable gsMapData<T> m_map, m_map_def;
    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T,dim,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def;
    mutable gsMatrix<T,3,2> m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;
    mutable gsMatrix<T> m_Emat,m_Nmat,m_Tmat,m_rhomat;
    mutable real_t m_lambda, m_mu, m_Cconstant;

    // Compressible material matrix
    mutable gsMatrix<T>                 m_deriv2_def, m_deriv2_ori;
    mutable gsMatrix<T,3,3>             m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T,3,3>             m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;


    mutable gsMatrix<T>                 m_parmat;
    mutable gsVector<T>                 m_parvals;
    mutable T                           m_J0, m_J0_sq, m_J, m_J_sq, m_Tval;
    mutable gsMatrix<T> m_stretches, m_stretchvec;

    gsOptionList m_options;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixLinear.hpp)
#endif
