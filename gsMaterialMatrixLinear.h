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

/** @brief Assembles system matrices for thin shell linear and nonlinear elasticity problems.

    \tparam T coefficient type

    \ingroup gsThinShell
*/

/*
Desired template parameters
- compressible                      bool    COM
- material model                    int     MAT
- output type                       int     OUT / VEC
- thickness integration method      int     INT
- geometric dimension               int     DIM
*/

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

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt



    /// Shared pointer for gsMaterialMatrixLinear
    typedef memory::shared_ptr< gsMaterialMatrixLinear > Ptr;

    /// Unique pointer for gsMaterialMatrixLinear
    typedef memory::unique_ptr< gsMaterialMatrixLinear > uPtr;

    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    short_t domainDim() const;

    // template OUT
    short_t targetDim() const;

    // template COM
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // template COM
    // void eval_into_vector(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // // template COM
    // void eval_into_matrix(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // // template COM
    // void eval_into_pstress(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // // template COM
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    // template COM
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    enum Material material() const {return Material::SvK;};

    void setMoment(const index_t moment) const {m_moment = moment;};

public:
    // template<short_t num=0>
    // void makeMatrix()                   { m_outputType=2; m_output = num;}

    void setParameters(const std::vector<gsFunction<T>*> &pars)
    {
        m_pars = pars;
        m_numPars = m_pars.size();
    }

    void info() const;

    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u) const;

    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u) const;

    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
    { GISMO_NO_IMPLEMENTATION; }
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u) const
    { GISMO_NO_IMPLEMENTATION; }

protected:
    void initialize();
    void defaultOptions();


protected:

    // template MAT
    T Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    // template MAT
    T Sij    (const index_t i, const index_t j, const T z) const;

    // template MAT
    void computePoints(const index_t patch, const gsMatrix<T> & u, bool deformed=true) const;
    // template DIM
    void computeMetricDeformed() const;
    // template DIM
    void computeMetricUndeformed() const;
    // template DIM
    void getMetric(index_t k, T z) const;
    // template DIM
    void getMetricDeformed(index_t k, T z) const;
    // template DIM
    void getMetricUndeformed(index_t k, T z) const;
    // template DIM
    void computeStretch(const gsMatrix<T> & C ) const;

    // template DIM
    std::pair<gsVector<T>,gsMatrix<T>> evalStretch(const gsMatrix<T> & C ) const;

    private:
        template<enum Material _mat>
        typename std::enable_if<_mat==Material::OG, void>::type computePoints_impl(const gsMatrix<T> & u, bool deformed) const;
        template<enum Material _mat>
        typename std::enable_if<_mat!=Material::OG, void>::type computePoints_impl(const gsMatrix<T> & u, bool deformed) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type computeMetricDeformed_impl() const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type computeMetricDeformed_impl() const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type computeMetricUndeformed_impl() const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type computeMetricUndeformed_impl() const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type getMetric_impl(index_t k, T z) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type getMetric_impl(index_t k, T z) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type getMetricDeformed_impl(index_t k, T z) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type getMetricDeformed_impl(index_t k, T z) const;

        template<short_t _dim>
        typename std::enable_if<_dim==2, void>::type getMetricUndeformed_impl(index_t k, T z) const;
        template<short_t _dim>
        typename std::enable_if<_dim==3, void>::type getMetricUndeformed_impl(index_t k, T z) const;

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
    mutable index_t m_moment;
    mutable gsMatrix<T> m_result;

    // constructor
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;
    const gsFunction<T> * m_thickness;
    mutable std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_density;

    // Linear material matrix
    mutable gsMapData<T> m_map, m_map_def;
    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T,3,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def, m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;
    mutable gsMatrix<T> m_Emat,m_Nmat,m_Tmat,m_rhomat;
    mutable real_t m_lambda, m_mu, m_Cconstant;


    // Composite material matrix
    const std::vector<std::pair<T,T>>   m_YoungsModuli;
    const std::vector<T>                m_ShearModuli;
    const std::vector<std::pair<T,T>>   m_PoissonRatios;
    const std::vector<T>                m_thickValues;
    const std::vector<T>                m_phis;
    const std::vector<T>                m_densities;
    mutable gsVector<T>                 m_e1, m_e2, m_ac1, m_ac2, m_a1, m_a2;
    mutable gsMatrix<T>                 m_covBasis, m_covMetric, m_conBasis, m_conMetric;
    mutable gsMatrix<T>                 m_stretches, m_stretchvec;
    mutable T                           m_E1, m_E2, m_G12, m_nu12, m_nu21, m_t, m_t_tot,
                                        m_t_temp, m_z, m_z_mid, m_phi, m_rho;
    mutable gsMatrix<T>                 m_Dmat, m_Transform;

    // Compressible material matrix
    mutable gsMatrix<T>                 m_deriv2_def, m_deriv2_ori;
    mutable gsMatrix<T,3,3>             m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T,3,3>             m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;
    mutable gsMatrix<T>                 m_parmat;
    mutable gsVector<T>                 m_parvals;
    mutable T                           m_J0, m_J0_sq, m_J, m_J_sq, m_Tval;
    // integrateZ
    mutable gsMatrix<T> m_points2D, m_points3D, m_evalPoints;
    // mutable gsMatrix<T> m_quNodes;
    // mutable gsVector<T> m_quWeights;
    mutable gsGaussRule<T> m_gauss;
    mutable index_t m_numGauss;
    mutable T m_tHalf;


    mutable gsOptionList m_options;

    mutable index_t m_material;
    mutable bool m_compressible;
    mutable int m_compFun;
    mutable int m_outputType, m_output;

private:
    static int delta(const int a, const int b)
    {
        return (a==b) ? 1 : 0;
    }

    static int idelta(const int a, const int b)
    {
        return (a!=b) ? 1 : 0;
    }

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixLinear.hpp)
#endif
