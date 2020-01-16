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
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const gsFunction<T> & YoungsModulus,
                        const gsFunction<T> & PoissonRatio);

    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const gsFunction<T> & YoungsModulus,
                        const gsFunction<T> & PoissonRatio);

    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi);


    /// @brief Returns the list of default options for assembly
    gsOptionList defaultOptions();

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
    void makeMatrix()     { m_matrix = true; }
    void makeVector()     { m_matrix = false; }

protected:
    void eval3D_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    gsMatrix<T> eval3D_Linear(const gsMatrix<T>& u) const;
    gsMatrix<T> eval_Composite(const gsMatrix<T>& u) const;
    gsMatrix<T> eval3D_Incompressible(const gsMatrix<T>& u) const;
    gsMatrix<T> eval3D_Compressible(const gsMatrix<T>& u) const;

    T Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    T Sij  (const index_t i, const index_t j) const;
    gsMatrix<T> S() const;

    // void eval_Linear(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // void eval_Incompressible(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // void eval_Compressible(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    gsMatrix<T> integrateZ(const gsMatrix<T>& u) const;
    gsMatrix<T> multiplyZ (const gsMatrix<T>& u) const;

    void computeMetric(index_t k, bool computedeformed=false, bool computefull=false, T z=0.0) const;
    void computePoints(const gsMatrix<T> & u, bool deformed=true) const;

protected:
    // general
    index_t m_model;
    index_t m_pIndex;
    index_t m_numParameters; // how many parameters for the material model?
    int m_moment;
    bool m_matrix;
    mutable gsMatrix<T> m_result;

    // constructor
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;
    const gsFunction<T> * m_thickness;
    const gsFunction<T> * m_par2;
    const gsFunction<T> * m_par1;


    // Linear material matrix
    mutable gsMapData<T> m_map, m_map_def;
    mutable gsMatrix<T> m_metricA, m_metricA_def;
    mutable gsMatrix<T> m_Emat,m_Nmat,m_Tmat;
    mutable real_t m_lambda, m_mu, m_Cconstant;

    // Composite material matrix
    const std::vector<std::pair<T,T>>   m_YoungsModuli;
    const std::vector<T>                m_ShearModuli;
    const std::vector<std::pair<T,T>>   m_PoissonRatios;
    const std::vector<T>                m_thickValues;
    const std::vector<T>                m_phis;
    mutable gsVector<T>                 m_normal, m_e1, m_e2, m_ac1, m_ac2, m_a1, m_a2;
    mutable gsMatrix<T>                 m_covBasis, m_covMetric, m_conBasis, m_conMetric;
    mutable T                           m_E1, m_E2, m_G12, m_nu12, m_nu21, m_t, m_t_tot,
                                        m_t_temp, m_z, m_z_mid, m_phi;
    mutable gsMatrix<T>                 m_Dmat, m_Transform;

    // Compressible material matrix
    mutable gsMatrix<T>                 m_deriv2_def, m_deriv2;
    mutable gsMatrix<T>                 m_metricB_def, m_metricB, m_metricG, m_metricG_def;
    mutable gsVector<T>                 m_normal_def;
    mutable gsMatrix<T>                 m_par1mat, m_par2mat;
    mutable T                           m_par1val, m_par2val, m_J0;

    // integrateZ
    mutable gsMatrix<T> m_points, m_evalPoints;
    mutable gsMatrix<T> m_quNodes;
    mutable gsVector<T> m_quWeights;
    mutable gsGaussRule<T> m_gauss;
    mutable index_t m_numGauss;
    mutable T m_tHalf;


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

