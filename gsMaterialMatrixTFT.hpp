/** @file gsMaterialMatrixTFT.hpp

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

#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>

#ifdef GISMO_WITH_IPOPT
#include <gsIpopt/gsOptProblem.h>
#endif

namespace gismo
{

#ifdef GISMO_WITH_IPOPT

template <typename T>
class gsOptProblemTFT : public gsOptProblem<T>
{
public:

    gsOptProblemTFT(const gsMatrix<T> & C, const gsMatrix<T> & e)
    :
    m_C(C),
    m_e(e)
    {
        m_numDesignVars  = 1;
        m_numConstraints = 1;
        m_numConJacNonZero = 1;

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        // theta has a lower bound of -1 and an upper bound of 1
    const T inf = std::numeric_limits<T>::infinity();
        m_desLowerBounds[0] = -0.5*M_PI;
        m_desUpperBounds[0] =  0.5*M_PI;

        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);

        // we have one equality constraint, so we set the bounds on
        // this constraint to be equal (and zero).
        m_conLowerBounds[0] = 0;
        m_conUpperBounds[0] = inf;

        // m_conLowerBounds[0] = -1;
        // m_conUpperBounds[0] = 1;

        // we initialize x in bounds, in the upper right quadrant
        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign(0,0) = 0.0;

        //
        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);
    }


public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        T theta = u(0,0);
        T n1 = math::cos(theta);
        T n2 = math::sin(theta);
        T m1 = -math::sin(theta);
        T m2 = math::cos(theta);

        gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<T,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<T,1,1> gamma = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );

        gsMatrix<T,1,1> result = n2_vec.transpose() * m_C * m_e + gamma * n2_vec.transpose() * m_C * n1_vec;
        GISMO_ASSERT(result.rows()==1 && result.cols() ==1,"f is not scalar!");

        return result(0,0) * result(0,0);
    }

    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        T theta = u(0,0);
        result.resize(m_numConstraints,1);
        // return the value of the constraints: g(x)
        real_t n1 = math::cos(theta);
        real_t n2 = math::sin(theta);
        real_t m1 = -math::sin(theta);
        real_t m2 = math::cos(theta);

        gsVector<real_t,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<real_t,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<real_t,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<real_t,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<real_t,1,1> gamma = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );
        result(0) = gamma(0,0);
    }

    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // at the moment only a full finite difference matrix is returned.

        gsVector<T> uu = u;
        gsVector<T> e1( this->m_numConstraints );
        gsVector<T> e2( this->m_numConstraints );
        gsAsVector<T> ee1( e1.data() , e1.rows() );
        gsAsVector<T> ee2( e2.data() , e2.rows() );

        index_t lastDesginVar = -1;

        // TODO: Replace by a better value or use AD...
        const T h = T(0.00001);

        for( index_t i = 0 ; i < this->m_numConJacNonZero; ++i )
        {
            index_t row = this->m_conJacRows[i];  // constrains
            index_t col = this->m_conJacCols[i];  // designVariables

            if( lastDesginVar != col )
            {
                gsAsConstVector<T> uuMap( uu.data() , uu.rows() );

                uu(col) -= h/2.;
                evalCon_into( uuMap, ee1);
                uu(col) += h;
                evalCon_into( uuMap, ee2 );
                uu(col) = u(col);

                lastDesginVar = col;
            }

            result(i) = (0.5*e1(row) - 0.5*e2(row)) / h;

        }

    }

private:

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;

    const gsMatrix<T> m_C;
    const gsMatrix<T> m_e;
};
#endif

// template <short_t dim, class T >
// gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
//                                         const gsFunctionSet<T> & mp,
//                                         const gsFunction<T> & thickness,
//                                         const std::vector<gsFunction<T>*> &pars
//                                         )
//                                         :
//                                         Base(mp),
//                                         m_thickness(&thickness),
//                                         m_pars(pars)
// {
//     _initialize();
// }

// template <short_t dim, class T >
// gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
//                                     const gsFunctionSet<T> & mp,
//                                     const gsFunction<T> & thickness,
//                                     const std::vector<gsFunction<T>*> &pars,
//                                     const gsFunction<T> & density
//                                     )
//                                     :
//                                     Base(mp),
//                                     m_thickness(&thickness),
//                                     m_pars(pars),
//                                     m_density(&density)
// {
//     _initialize();
// }

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixTFT<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    gsMatrix<T> Cs = m_materialMat->eval3D_matrix(patch,u,z,out);

    gsMatrix<T> Es;
    if (out == MaterialOutput::MatrixA)
    {
        gsDebugVar("MatrixA");
        Es = m_materialMat->eval3D_strain(patch,u,z,MaterialOutput::VectorN); // returns flexural strain (although type is VectorN)
    }
    else if (out == MaterialOutput::MatrixD)
    {
        gsDebugVar("MatrixD");
        Es = m_materialMat->eval3D_strain(patch,u,z,MaterialOutput::VectorM); // returns flexural strain (although type is VectorM)
    }
    else
    {
        gsDebugVar("MatrixB/MatrixC");
        Es = gsMatrix<T>::Zero(3,Cs.cols()); // FOR LINEAR MATERIALS
    }

    gsDebugVar(Es);

    gsMatrix<T> TF = m_materialMat->eval3D_tensionfield(patch,u,z,MaterialOutput::TensionField);

    gsMatrix<T> result = Cs;

#ifdef GISMO_WITH_IPOPT

    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            if (TF(0,j*u.cols() + k) == 1) // taut
            {
                // do nothing
            }
            else if (TF(0,j*u.cols() + k) == -1) // slack
            {
                // Set to zero
                result.col(j*u.cols()).setZero();
            }
            else if (TF(0,j*u.cols() + k) == 0) // wrinkled
            {
                // Make modified form
                gsAsMatrix<T, Dynamic, Dynamic> res = result.reshapeCol(j*u.cols()+k,3,3);

                gsMatrix<T> C = Cs.reshapeCol(j*u.cols()+k,3,3);
                gsMatrix<T> E = Es.col(j*u.cols()+k);

                gsOptProblemTFT<real_t> opt(C,E);
                opt.solve();

                T theta = opt.currentDesign()(0,0);

                T n1 = math::cos(theta);
                T n2 = math::sin(theta);
                T m1 = -math::sin(theta);
                T m2 = math::cos(theta);
                gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
                gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
                gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;
                gsMatrix<T,1,1> denum = n1_vec.transpose() * C * n1_vec;
                gsMatrix<T,3,3> C_I = C - 1 / (  n1_vec.transpose() * C * n1_vec ) * C * ( n1_vec * n1_vec.transpose() ) * C;
                gsMatrix<T,1,1> gamma = - ( n1_vec.transpose() * C * E ) / ( n1_vec.transpose() * C * n1_vec );
                gsMatrix<T,1,1> tmp2 = (n1_vec.transpose() * C * n2_vec);

                GISMO_ASSERT(tmp2.rows()==1 && tmp2.cols()==1,"Must be scalar");

                gsMatrix<T,1,1> df = n4_vec.transpose() * C * (E + gamma(0,0) * n1_vec)
                                    + 2 * gamma * ( n2_vec.transpose() * C * n2_vec
                                    - math::pow(tmp2(0,0),2) / (n1_vec.transpose() * C * n1_vec) );

                gsMatrix<T,3,1> b = n2_vec - ( (n1_vec.transpose() * C * n2_vec)(0,0) / ( n1_vec.transpose() * C * n1_vec )(0,0)) * n1_vec;
                gsMatrix<T,3,3> C_II = C_I + 2 * gamma(0,0) / df(0,0) * (C * b * b.transpose() * C);
                res = C_II;
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,j*u.cols()));
        }
    }

#else
    GISMO_ERROR("IpOpt not enabled");
#endif

    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixTFT<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    // WARNING: Check which values of out correspond to the output we want
    gsMatrix<T> Cs;
    if (out == MaterialOutput::VectorN)
        Cs = this->eval3D_matrix(patch,u,z,MaterialOutput::MatrixA);
    else if (out == MaterialOutput::VectorM)
        Cs = this->eval3D_matrix(patch,u,z,MaterialOutput::MatrixD);
    else
        GISMO_ERROR("Output type not understood");

    gsMatrix<T> Es = m_materialMat->eval3D_strain(patch,u,z,out);

    gsMatrix<T> result = Es;

#ifdef GISMO_WITH_IPOPT

    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            result.col(j*u.cols()+k) = Cs.reshapeCol(j*u.cols()+k,3,3) * Es.col(j*u.cols()+k);
        }
    }

#else
    GISMO_ERROR("IpOpt not enabled");
#endif

    return result;
}

} // end namespace
