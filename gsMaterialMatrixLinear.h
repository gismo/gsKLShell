/** @file gsMaterialMatrixLinear.h

    @brief Provides linear material matrices

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
#include <gsKLShell/gsMaterialMatrixBaseDim.h>
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
class gsMaterialMatrixLinear : public gsMaterialMatrixBaseDim<dim,T>
{
public:
    using Base = gsMaterialMatrixBaseDim<dim,T>;

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
     * @brief      Full constructor
     *
     * @param[in]  mp         Original geometry
     * @param[in]  mp_def     Deformed geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
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
    void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    /// See \ref gsMaterialMatrixBase for details
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

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
    void _computePoints(const index_t patch, const gsMatrix<T> & u, bool basis = true) const;

    /// Computes the stretch given deformation tensor C, into class members m_stretches and m_stretchDirs
    void _computePStress(const gsMatrix<T> & C ) const;

    /// Computes the stretch given deformation tensor C, into a pair
    std::pair<gsVector<T>,gsMatrix<T>> _evalPStress(const gsMatrix<T> & C ) const;

protected:
    // general
    index_t m_numPars; // how many parameters for the material model?

    // constructor
    using Base::m_patches;
    using Base::m_defpatches;
    const gsFunction<T> * m_thickness;
    std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_density;

    mutable gsMatrix<T> m_Emat,m_Nmat,m_Tmat,m_rhomat;
    mutable real_t m_lambda, m_mu, m_Cconstant;

    mutable gsMatrix<T>                 m_parmat;
    mutable gsVector<T>                 m_parvals;

    mutable gsMatrix<T> m_pstress, m_pstressvec;

    // Geometric data point
    using Base::m_map;
    using Base::m_map_def;

    using Base::m_Acov_ori;
    using Base::m_Acon_ori;
    using Base::m_Acov_def;
    using Base::m_Acon_def;
    using Base::m_Bcov_ori;
    using Base::m_Bcon_ori;
    using Base::m_Bcov_def;
    using Base::m_Bcon_def;
    using Base::m_acov_ori;
    using Base::m_acon_ori;
    using Base::m_acov_def;
    using Base::m_acon_def;
    using Base::m_ncov_ori;
    using Base::m_ncov_def;
    using Base::m_Gcov_ori;
    using Base::m_Gcon_ori;
    using Base::m_Gcov_def;
    using Base::m_Gcon_def;
    using Base::m_Gcov_ori_L;
    using Base::m_Gcov_def_L;
    using Base::m_gcov_ori;
    using Base::m_gcov_def;
    using Base::m_gcon_ori;
    using Base::m_gcon_def;
    using Base::m_Acov_ori_mat;
    using Base::m_Acon_ori_mat;
    using Base::m_Acov_def_mat;
    using Base::m_Acon_def_mat;
    using Base::m_Bcov_ori_mat;
    using Base::m_Bcov_def_mat;
    using Base::m_acov_ori_mat;
    using Base::m_acon_ori_mat;
    using Base::m_acov_def_mat;
    using Base::m_acon_def_mat;
    using Base::m_ncov_ori_mat;
    using Base::m_ncov_def_mat;

    using Base::m_stretches;
    using Base::m_stretchvec;

    using Base::m_J0_sq;
    using Base::m_J_sq;



    gsOptionList m_options;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixLinear.hpp)
#endif
