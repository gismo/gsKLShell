/** @file gsMaterialMatrix.h

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

namespace gismo
{

/** @brief Assembles system matrices for thin shell linear and nonlinear elasticity problems.

    \tparam T coefficient type

    \ingroup gsThinShell
*/
template <class T>
class gsMaterialMatrix : public gismo::gsFunction<T>
{
public:
    /** @brief Constructor of the assembler object.

        \param[in] ...
        \param[in] ...

    */
    // Without deformed geometry
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const gsFunction<T> & YoungsModulus,
                        const gsFunction<T> & PoissonRatio);
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const gsFunction<T> & YoungsModulus,
                        const gsFunction<T> & PoissonRatio,
                        const gsFunction<T> & Density);

    // With deformed geometry
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const gsFunction<T> & YoungsModulus,
                        const gsFunction<T> & PoissonRatio);
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const gsFunction<T> & YoungsModulus,
                        const gsFunction<T> & PoissonRatio,
                        const gsFunction<T> & Density);

    // Laminates without deformed
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi);
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi,
                        const std::vector<T>                density);

    // Laminates with deformed
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const gsFunctionSet<T>            & mp_def,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi);
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const gsFunctionSet<T>            & mp_def,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi,
                        const std::vector<T>                density);


    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}



    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // CORRECT???
    // ~gsMaterialMatrix() { delete m_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrix)

    short_t domainDim() const;

    short_t targetDim() const;

    mutable gsMaterialMatrix<T> * m_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t p) const
    {
        // delete m_piece;
        // m_piece = new gsMaterialMatrix(m_patches->piece(k), *m_thickness, *m_YoungsModulus, *m_PoissonRatio);
        m_piece = new gsMaterialMatrix(*this);
        m_piece->setPatch(p);
        return *m_piece;
    }

    ~gsMaterialMatrix() { delete m_piece; }

    void setPatch(index_t p) {m_pIndex = p; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

public:
    void setMoment(int m)   { m_moment = m; }
    void makeMatrix()     { m_output=2; }
    void makeVector()     { m_output=1; }
    void makeDensity()    { m_output=0; }

protected:
    void initialize();
    void defaultOptions();
    void getOptions() const;

    gsMatrix<T> eval3D(const index_t i, const gsMatrix<T>& z) const;
    gsMatrix<T> eval3D(const index_t i) const;
    gsMatrix<T> eval_Compressible(const index_t i, const gsMatrix<T>& z) const;
    gsMatrix<T> eval_Incompressible(const index_t i, const gsMatrix<T>& z) const;
    gsMatrix<T> eval_Composite(const gsMatrix<T>& u) const;
    gsMatrix<T> eval3D_Incompressible(const gsMatrix<T>& u) const;
    gsMatrix<T> eval3D_Compressible(const gsMatrix<T>& u) const;

    T Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    T Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T Sij    (const index_t i, const index_t j) const;
    T Sij    (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    T dPsi   (const index_t i, const index_t j) const;
    T dPsi   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    T d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    // T Sij  (const index_t i, const index_t j, gsMatrix<T> & cinv) const;

    // void eval_Linear(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // void eval_Incompressible(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // void eval_Compressible(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    gsMatrix<T> integrateZ(const gsMatrix<T>& u) const;
    gsMatrix<T> multiplyZ (const gsMatrix<T>& u) const;

    void computeMetricDeformed(const gsMatrix<T>& u) const;
    void computeMetricUndeformed(const gsMatrix<T>& u) const;
    void getMetric(index_t k, T z) const;
    void getMetricDeformed(index_t k, T z) const;
    void getMetricUndeformed(index_t k, T z) const;

    void computeStretch(const gsMatrix<T> & C ) const;
    void computePoints(const gsMatrix<T> & u, bool deformed=true) const;

protected:
    // general
    index_t m_pIndex;
    index_t m_numParameters; // how many parameters for the material model?
    int m_moment;
    mutable gsMatrix<T> m_result;

    // constructor
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;
    const gsFunction<T> * m_thickness;
    const gsFunction<T> * m_par2;
    const gsFunction<T> * m_par1;
    const gsFunction<T> * m_density;


    // Linear material matrix
    mutable gsMapData<T> m_map, m_map_def;
    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
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
    mutable gsMatrix<T,3,3>             m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def;
    mutable gsMatrix<T>                 m_par1mat, m_par2mat;
    mutable T                           m_par1val, m_par2val, m_J0, m_J;
    // integrateZ
    mutable gsMatrix<T> m_points2D, m_points3D, m_evalPoints;
    // mutable gsMatrix<T> m_quNodes;
    // mutable gsVector<T> m_quWeights;
    mutable gsGaussRule<T> m_gauss;
    mutable index_t m_numGauss;
    mutable T m_tHalf;


    mutable gsOptionList m_options;

    /// @brief Specifies the material law to use
    struct material_law
    {
        enum type
        {
            SvK_Isotropic = 0,          /// Psi = ........ S = 2*mu*E + lambda*tr(E)*I
            SvK_Orthotropic = 1,        /// Psi = ........ S = 2*mu*E + lambda*tr(E)*I
            NHK = 2,       /// Psi = ........ S = lambda*ln(J)*C^-1 + mu*(I-C^-1)
            MR = 3        /// Psi = ........ S = lambda*ln(J)*C^-1 + mu*(I-C^-1)
        };
    };
    /// @brief Specifies (in)compressibility
    struct compressibility
    {
        enum type
        {
            incompressible = 0,
            compressible = 1
        };
    };

    mutable index_t m_material;
    mutable bool m_compressible;
    mutable int m_output;



};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrix.hpp)
#endif


// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
// template <class T>
// class gsMaterialMatrix : public gismo::gsFunction<T>
// {
//   // Computes the material matrix for different material models
//   //
// protected:
//     const gsFunctionSet<T> * _mp;
//     const gsFunction<T> * _YoungsModulus;
//     const gsFunction<T> * _PoissonRatio;
//     mutable gsMapData<T> _tmp;
//     mutable gsMatrix<real_t,3,3> F0;
//     mutable gsMatrix<T> Emat,Nmat;
//     mutable real_t lambda, mu, E, nu, C_constant;

// public:
//     /// Shared pointer for gsMaterialMatrix
//     typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

//     /// Unique pointer for gsMaterialMatrix
//     typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//     gsMaterialMatrix(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
//                    const gsFunction<T> & PoissonRatio) :
//     _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio), _mm_piece(nullptr)
//     {
//         _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
//     }

//     ~gsMaterialMatrix() { delete _mm_piece; }

//     GISMO_CLONE_FUNCTION(gsMaterialMatrix)

//     short_t domainDim() const {return 2;}

//     short_t targetDim() const {return 9;}

//     mutable gsMaterialMatrix<T> * _mm_piece; // todo: improve the way pieces are accessed

//     const gsFunction<T> & piece(const index_t k) const
//     {
//         delete _mm_piece;
//         _mm_piece = new gsMaterialMatrix(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
//         return *_mm_piece;
//     }

//     //class .. matMatrix_z
//     // should contain eval_into(thickness variable)

//     // Input is parametric coordinates of the surface \a mp
//     void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
//     {
//         // NOTE 1: if the input \a u is considered to be in physical coordinates
//         // then we first need to invert the points to parameter space
//         // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
//         // otherwise we just use the input paramteric points
//         _tmp.points = u;

//         static_cast<const gsFunction<T>&>(_mp->piece(0)).computeMap(_tmp); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

//         // NOTE 2: in the case that parametric value is needed it suffices
//         // to evaluate Youngs modulus and Poisson's ratio at
//         // \a u instead of _tmp.values[0].
//         _YoungsModulus->eval_into(_tmp.values[0], Emat);
//         _PoissonRatio->eval_into(_tmp.values[0], Nmat);

//         result.resize( targetDim() , u.cols() );
//         for( index_t i=0; i< u.cols(); ++i )
//         {
//             gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

//             F0.leftCols(2) = _tmp.jacobian(i);
//             F0.col(2)      = _tmp.normal(i).normalized();
//             F0 = F0.inverse();
//             F0 = F0 * F0.transpose(); //3x3

//             // Evaluate material properties on the quadrature point
//             E = Emat(0,i);
//             nu = Nmat(0,i);
//             lambda = E * nu / ( (1. + nu)*(1.-2.*nu)) ;
//             mu     = E / (2.*(1. + nu)) ;

//             C_constant = 2*lambda*mu/(lambda+2*mu);

//             C(0,0) = C_constant*F0(0,0)*F0(0,0) + 1*mu*(2*F0(0,0)*F0(0,0));
//             C(1,1) = C_constant*F0(1,1)*F0(1,1) + 1*mu*(2*F0(1,1)*F0(1,1));
//             C(2,2) = C_constant*F0(0,1)*F0(0,1) + 1*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
//             C(1,0) =
//             C(0,1) = C_constant*F0(0,0)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(0,1));
//             C(2,0) =
//             C(0,2) = C_constant*F0(0,0)*F0(0,1) + 1*mu*(2*F0(0,0)*F0(0,1));
//             C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(1,1));

//             //gsDebugVar(C);
//         }
//     }

//     // piece(k) --> for patch k

// }; //! [Include namespace]

