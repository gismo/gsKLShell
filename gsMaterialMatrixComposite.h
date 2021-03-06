/** @file gsMaterialMatrixComposite.h

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
class gsMaterialMatrixComposite : public gsMaterialMatrixBase<T>
{
public:
    /** @brief Constructor of the assembler object.

        \param[in] ...
        \param[in] ...

    */
    /// Default empty constructor
        gsMaterialMatrixComposite() { }

        gsMaterialMatrixComposite(
                                const gsFunctionSet<T>                              & mp,
                                const gsFunction<T>                                 & thickness,
                                const gsFunction<T>                                 & E11,
                                const gsFunction<T>                                 & E22,
                                const gsFunction<T>                                 & G12,
                                const gsFunction<T>                                 & nu12,
                                const gsFunction<T>                                 & nu21,
                                const gsFunction<T>                                 & phi           );

        gsMaterialMatrixComposite(
                                const gsFunctionSet<T>                              & mp,
                                const gsFunctionSet<T>                              & mp_def,
                                const gsFunction<T>                                 & thickness,
                                const gsFunction<T>                                 & E11,
                                const gsFunction<T>                                 & E22,
                                const gsFunction<T>                                 & G12,
                                const gsFunction<T>                                 & nu12,
                                const gsFunction<T>                                 & nu21,
                                const gsFunction<T>                                 & phi           );

        gsMaterialMatrixComposite(
                                const gsFunctionSet<T>                              & mp,
                                const gsFunction<T>                                 & thickness,
                                const gsFunction<T>                                 & E11,
                                const gsFunction<T>                                 & E22,
                                const gsFunction<T>                                 & G12,
                                const gsFunction<T>                                 & nu12,
                                const gsFunction<T>                                 & nu21,
                                const gsFunction<T>                                 & phi,
                                const gsFunction<T>                                 & rho           );

        gsMaterialMatrixComposite(
                                const gsFunctionSet<T>                              & mp,
                                const gsFunctionSet<T>                              & mp_def,
                                const gsFunction<T>                                 & thickness,
                                const gsFunction<T>                                 & E11,
                                const gsFunction<T>                                 & E22,
                                const gsFunction<T>                                 & G12,
                                const gsFunction<T>                                 & nu12,
                                const gsFunction<T>                                 & nu21,
                                const gsFunction<T>                                 & phi,
                                const gsFunction<T>                                 & rho           );

        gsMaterialMatrixComposite(
                                const gsFunctionSet<T>                              & mp,
                                const gsFunction<T>                                 & thickness,
                                const std::vector<gsFunction<T>*>                   & pars          );

        gsMaterialMatrixComposite(
                                const gsFunctionSet<T>                              & mp,
                                const gsFunctionSet<T>                              & mp_def,
                                const gsFunction<T>                                 & thickness,
                                const std::vector<gsFunction<T>*>                   & pars          );

        gsMaterialMatrixComposite(
                                const gsFunctionSet<T>                              & mp,
                                const gsFunctionSet<T>                              & mp_def,
                                const gsFunction<T>                                 & thickness,
                                const std::vector<gsFunction<T>*>                   & pars,
                                const gsFunction<T>                                 & rho           );

        // // ---------------------------------------------------------------------------------------------

        // gsMaterialMatrixComposite(
        //                         const gsFunctionSet<T>                              & mp,
        //                         const gsVector<T>                                   & E11,
        //                         const gsVector<T>                                   & E22,
        //                         const gsVector<T>                                   & G12,
        //                         const gsVector<T>                                   & nu12,
        //                         const gsVector<T>                                   & nu21,
        //                         const gsVector<T>                                   & thickness,
        //                         const gsVector<T>                                   & phi           );

        // gsMaterialMatrixComposite(
        //                         const gsFunctionSet<T>                              & mp,
        //                         const gsFunctionSet<T>                              & mp_def,
        //                         const gsVector<T>                                   & E11,
        //                         const gsVector<T>                                   & E22,
        //                         const gsVector<T>                                   & G12,
        //                         const gsVector<T>                                   & nu12,
        //                         const gsVector<T>                                   & nu21,
        //                         const gsVector<T>                                   & thickness,
        //                         const gsVector<T>                                   & phi           );

        // gsMaterialMatrixComposite(
        //                         const gsFunctionSet<T>                              & mp,
        //                         const gsVector<T>                                   & E11,
        //                         const gsVector<T>                                   & E22,
        //                         const gsVector<T>                                   & G12,
        //                         const gsVector<T>                                   & nu12,
        //                         const gsVector<T>                                   & nu21,
        //                         const gsVector<T>                                   & thickness,
        //                         const gsVector<T>                                   & phi,
        //                         const gsVector<T>                                   & rho           );

        // gsMaterialMatrixComposite(
        //                         const gsFunctionSet<T>                              & mp,
        //                         const gsFunctionSet<T>                              & mp_def,
        //                         const gsVector<T>                                   & E11,
        //                         const gsVector<T>                                   & E22,
        //                         const gsVector<T>                                   & G12,
        //                         const gsVector<T>                                   & nu12,
        //                         const gsVector<T>                                   & nu21,
        //                         const gsVector<T>                                   & thickness,
        //                         const gsVector<T>                                   & phi,
        //                         const gsVector<T>                                   & rho           );

    enum MatIntegration isMatIntegrated() const {return MatIntegration::Integrated; }
    enum MatIntegration isVecIntegrated() const {return MatIntegration::Integrated; }

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt

    // template COM
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // template COM
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
    { GISMO_NO_IMPLEMENTATION; }
    void info() const;

public:
    /// Shared pointer for gsMaterialMatrixComposite
    typedef memory::shared_ptr< gsMaterialMatrixComposite > Ptr;

    /// Unique pointer for gsMaterialMatrixComposite
    typedef memory::unique_ptr< gsMaterialMatrixComposite > uPtr;

protected:
    void initialize();
    void initializeParameters();
    void defaultOptions();

    gsMatrix<T> computeMatrix(const T E11, const T E22, const T G12, const T nu12, const T nu21, const T phi) const;
    gsMatrix<T> cart2cov(const gsVector<T> a1, const gsVector<T> a2, const gsVector<T> e1, const gsVector<T> e2) const;
    gsMatrix<T> con2cart(const gsVector<T> ac1, const gsVector<T> ac2, const gsVector<T> e1, const gsVector<T> e2) const;


protected:

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
    const gsFunction<T> * m_E11;
    const gsFunction<T> * m_E22;
    const gsFunction<T> * m_G12;
    const gsFunction<T> * m_nu12;
    const gsFunction<T> * m_nu21;
    const gsFunction<T> * m_phi;
    const gsFunction<T> * m_rho;
    mutable std::vector<gsFunction<T>* > m_pars;


    // Linear material matrix
    mutable gsMapData<T> m_map, m_map_def;

    // Compute points
    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T,3,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def, m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;

    // Composite
    index_t m_nLayers;
    mutable gsVector<T>                 m_e1, m_e2, m_ac1, m_ac2, m_a1, m_a2;
    mutable gsMatrix<T>                 m_covBasis, m_covMetric, m_conBasis, m_conMetric;
    mutable gsMatrix<T>                 m_stretches, m_stretchvec;
    mutable T                           m_t, m_t_tot,m_t_temp, m_z, m_z_mid;
    mutable gsMatrix<T>                 m_Dmat, m_Transform;
    mutable gsMatrix<T>                 m_Tmat, m_E1mat, m_E2mat, m_G12mat, m_nu12mat, m_nu21mat, m_phiMat, m_rhoMat;



    mutable gsMatrix<T>                 m_deriv2_def, m_deriv2_ori;
    mutable gsMatrix<T,3,3>             m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T,3,3>             m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;

    mutable T                           m_J0, m_J0_sq, m_J, m_J_sq, m_Tval;

    mutable gsOptionList m_options;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixComposite.hpp)
#endif
