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


template <  short_t dim,
            class T
         >
class gsMaterialMatrixLinear : public gsMaterialMatrixBase<T>
{
public:
    /** @brief Constructor of the assembler object.

        \param[in] ...
        \param[in] ...

    */
    /// Default empty constructor
    gsMaterialMatrixLinear() { }

    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars,
                        const gsFunction<T> & Density);

    inline enum MatIntegration isMatIntegrated() const {return MatIntegration::Constant; }
    inline enum MatIntegration isVecIntegrated() const {return MatIntegration::Constant; }

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt

    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const
    { GISMO_NO_IMPLEMENTATION; }

    void setParameters(const std::vector<gsFunction<T>*> &pars)
    {
        m_pars = pars;
        m_numPars = m_pars.size();
    }

    void info() const;

public:
    /// Shared pointer for gsMaterialMatrixLinear
    typedef memory::shared_ptr< gsMaterialMatrixLinear > Ptr;

    /// Unique pointer for gsMaterialMatrixLinear
    typedef memory::unique_ptr< gsMaterialMatrixLinear > uPtr;

protected:
    void _initialize();
    void _defaultOptions();

protected:

    T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    T _Sij    (const index_t i, const index_t j, const T z, enum MaterialOutput out) const;

    void _computePoints(const index_t patch, const gsMatrix<T> & u, bool deformed=true) const;
    void _computeMetricDeformed() const;
    void _computeMetricUndeformed() const;
    void _getMetric(index_t k, T z) const;
    void _getMetricDeformed(index_t k, T z) const;
    void _getMetricUndeformed(index_t k, T z) const;

    private:
        template<enum Material _mat>
        typename std::enable_if<_mat==Material::OG, void>::type _computePoints_impl(const gsMatrix<T> & u, bool deformed) const;
        template<enum Material _mat>
        typename std::enable_if<_mat!=Material::OG, void>::type _computePoints_impl(const gsMatrix<T> & u, bool deformed) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type _computeMetricDeformed_impl() const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type _computeMetricDeformed_impl() const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type _computeMetricUndeformed_impl() const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type _computeMetricUndeformed_impl() const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type _getMetric_impl(index_t k, T z) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type _getMetric_impl(index_t k, T z) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type _getMetricDeformed_impl(index_t k, T z) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type _getMetricDeformed_impl(index_t k, T z) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type _getMetricUndeformed_impl(index_t k, T z) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type _getMetricUndeformed_impl(index_t k, T z) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type computeStretch_impl(const gsMatrix<T> & C ) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type computeStretch_impl(const gsMatrix<T> & C ) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, std::pair<gsVector<T>,gsMatrix<T>>>::type evalStretch_impl(const gsMatrix<T> & C ) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, std::pair<gsVector<T>,gsMatrix<T>>>::type evalStretch_impl(const gsMatrix<T> & C ) const;

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
    mutable gsMatrix<T,3,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def, m_ncov_ori, m_ncov_def;
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

    gsOptionList m_options;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixLinear.hpp)
#endif
