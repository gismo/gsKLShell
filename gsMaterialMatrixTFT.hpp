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

using namespace gismo;

template <typename T>
class objective : public gsFunction<T>
{
public:
    objective(const gsMatrix<T> & C, const gsMatrix<T> & S)
    :
    m_C(C),
    m_S(S)
    {
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsVector<T,3> n1_vec, n2_vec;
        T theta, gamma;
        T n1, n2, m1, m2;

        for (index_t k = 0; k!=u.cols(); k++)
        {
            theta = u(0,k);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            gamma = - ( n1_vec.transpose() * m_S ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();

            T res = (n2_vec.transpose() * m_S).value() + gamma * (n2_vec.transpose() * m_C * n1_vec).value();

            result(0,k) = res;
        }
    }

    short_t domainDim() const
    {
        return 1;
    }

    short_t targetDim() const
    {
        return 1;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsVector<T,3> n1_vec, n2_vec, n3_vec, n4_vec;
        T theta, gamma, dgammadT;
        T n1, n2, m1, m2;

        for (index_t k = 0; k!=u.cols(); k++)
        {
            theta = u(0,k);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
            n4_vec<<m1*m1, m2*m2, 2*m1*m2;

            gamma = - ( n1_vec.transpose() * m_S ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();
            T tmp = ( n2_vec.transpose() * m_C * n1_vec ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();
            dgammadT = - 2 * gamma * tmp;
            // gsDebugVar(dgammadT);

            T res = (n3_vec.transpose() * m_S).value() + dgammadT * (n2_vec.transpose() * m_C * n1_vec).value()
                  + gamma * (n3_vec.transpose() * m_C * n1_vec).value() + 2*gamma * (n2_vec.transpose() * m_C * n2_vec).value();
            // T res = (n2_vec.transpose() * m_S).value() + gamma * (n2_vec.transpose() * m_C * n1_vec).value();

            result(0,k) = res;
        }
    }

private:
    const gsMatrix<T> m_C;
    const gsMatrix<T> m_S;
};

template <typename T>
class gsTFTMat
{
public:

    gsTFTMat()
    { }


    gsTFTMat(const gsMatrix<T> & C, const gsMatrix<T> & e)
    :
    m_C(C),
    m_e(e)
    { }

    void compute(T t)
    {
        this->compute(t,m_C,m_e);
    }

    void compute(T t, const gsMatrix<T> & C, const gsMatrix<T> & e )
    {
        theta = t;
        T n1 = math::cos(theta);
        T n2 = math::sin(theta);
        T m1 = -math::sin(theta);
        T m2 = math::cos(theta);
        gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<T,1,1> denum = n1_vec.transpose() * C * n1_vec;

        C_I = C - 1 / (  n1_vec.transpose() * C * n1_vec ) * C * ( n1_vec * n1_vec.transpose() ) * C;

        T tmp = ( n2_vec.transpose() * C * n1_vec ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

        gamma = - ( n1_vec.transpose() * C * e ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

        // gsDebugVar((n1_vec.transpose() * C * n1_vec).value() + gamma * (n1_vec.transpose() * C * n1_vec).value());
        // gsDebugVar((n2_vec.transpose() * C * n1_vec).value() + gamma * (n2_vec.transpose() * C * n1_vec).value());

        gsMatrix<T> I(3,3); I.setIdentity();

        dgammadE = - ( n1_vec.transpose() * C * I ) / ( n1_vec.transpose() * C * n1_vec ).value();
        dgammadE.transposeInPlace();

        dgammadT = - 2 * gamma * tmp;

        dfdE = ( n2_vec.transpose() - tmp * n1_vec.transpose() ) * C * I;
        dfdE.transposeInPlace();

        dfdT = (n4_vec.transpose() * C * ( e + gamma * n1_vec )).value() +
                2 * gamma * ( n2_vec.transpose() * C * n2_vec - math::pow(( n1_vec.transpose() * C * n2_vec ).value(),2) / (n1_vec.transpose() * C * n1_vec).value() );

        dTdE = - dfdE / dfdT;

        gsMatrix<T,1,1> tmp2 = (n1_vec.transpose() * C * n2_vec);

        T df =  (n4_vec.transpose() * C * (e + gamma * n1_vec)).value()
                + 2 * gamma * ( (n2_vec.transpose() * C * n2_vec).value()
                - math::pow(tmp2(0,0),2) / (n1_vec.transpose() * C * n1_vec)
                )
                                ;

        GISMO_ASSERT(tmp2.rows()==1 && tmp2.cols()==1,"Must be scalar");

        gsMatrix<T,3,1> b = n2_vec - tmp * n1_vec;

        gsMatrix<T> C_I2 = C + C * (n1_vec * dgammadE.transpose());

        C_II = C_I + 2 * gamma / df * (C * b * b.transpose() * C);
        Ew = gamma * n1_vec;
        E  = e + Ew;
        Sp = C * (e + Ew);

        dEw = n1_vec.transpose() * dgammadE
                + dgammadT * n1_vec.transpose() * dTdE
                + 2 * gamma * n2_vec.transpose() * dTdE;
    }

public:
    mutable gsMatrix<T> C_I, C_II;
    mutable T gamma, theta;
    mutable gsMatrix<T> Ew, E, Sp, dEw;
    mutable gsMatrix<T> dgammadE, dfdE, dTdE;
    mutable T dgammadT, dfdT;

private:
    const gsMatrix<T> m_C;
    const gsVector<T> m_e;
};


namespace gismo
{

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // GISMO_ASSERT(out==MaterialOutput::MatrixA,"Tension Field Theory only works for membrane models, hence only outputs the A matrix");
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    return this->_eval3D_matrix_impl<true>(patch,u,z);
    // return this->_eval3D_matrix_impl<linear>(patch,u,z);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // GISMO_ASSERT(out==MaterialOutput::VectorN,"Tension Field Theory only works for membrane models, hence only outputs the N vector");

    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    // WARNING: Check which values of out correspond to the output we want
    return this->_eval3D_vector_impl<true>(patch,u,z);
    // return this->_eval3D_vector_impl<linear>(patch,u,z);
}

/// Provides an implementation of eval3D_matrix for \a gsMaterialMatrixLinear
template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if<_linear, gsMatrix<T> >::type
gsMaterialMatrixTFT<dim,T,linear>::_eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    gsMatrix<T> result = m_materialMat->eval3D_matrix(patch,u,z,out);
    return result;
    gsMatrix<T> TF = m_materialMat->eval3D_tensionfield(patch,u,z,MaterialOutput::TensionField);

    index_t colIdx;
    T gamma, theta;
    gsTFTMat<T> tftData;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j*u.cols()+k;
            if (TF(0,colIdx) == 1) // taut
            {
                // do nothing
            }
            else if (TF(0,colIdx) == -1) // slack
            {
                // Set to zero
                result.col(colIdx).setZero();
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                gsMatrix<T> N = m_materialMat->eval3D_vector(patch,u.col(k),z(j,k),MaterialOutput::VectorN);
                gsMatrix<T> thetas = eval_theta(C,N);

                gsMatrix<T> dC = m_materialMat->eval3D_dmatrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                theta = thetas(0,0);
                result.col(colIdx) = this->_compute_C(theta,C.reshape(3,3),N.reshape(3,1),dC.reshape(9,3)).reshape(9,1);
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if<!_linear, gsMatrix<T> >::type
gsMaterialMatrixTFT<dim,T,linear>::_eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    GISMO_NO_IMPLEMENTATION;
}

/// Provides an implementation of eval3D_matrix for \a gsMaterialMatrixLinear
template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if<_linear, gsMatrix<T> >::type
gsMaterialMatrixTFT<dim,T,linear>::_eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    gsMatrix<T> result = m_materialMat->eval3D_vector(patch,u,z,out);
    return result;
    gsMatrix<T> TF = m_materialMat->eval3D_tensionfield(patch,u,z,MaterialOutput::TensionField);

    index_t colIdx;
    T gamma, theta;
    gsTFTMat<T> tftData;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j*u.cols()+k;
            if (TF(0,colIdx) == 1) // taut
            {
                // do nothing
            }
            else if (TF(0,colIdx) == -1) // slack
            {
                // Set to zero
                result.col(colIdx).setZero();
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                gsMatrix<T> N = m_materialMat->eval3D_vector(patch,u.col(k),z(j,k),MaterialOutput::VectorN);
                gsMatrix<T> thetas = eval_theta(C,N);
                theta = thetas(0,0);
                gsDebugVar(N);
                gsDebugVar(C);
                result.col(colIdx) = this->_compute_S(theta,C,N);
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    gsDebugVar(result);

    return result;
}

template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if< !_linear, gsMatrix<T> >::type
gsMaterialMatrixTFT<dim,T,linear>::_eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    gsDebugVar("Bye");
    // gsDebugVar(MaterialMat::Linear);
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval_theta(const gsMatrix<T> & Cs, const gsMatrix<T> & Ns)
{
    GISMO_ASSERT(Cs.cols()==Ns.cols(),"Number of C matrices and N vectors is different");
    gsMatrix<T> result(1,Cs.cols());
    gsVector<T> value(1),init(1);
    gsVector<T,3> n1_vec, n2_vec;
    T n1, n2, m1, m2;

    index_t iter;

    T theta;
    T gamma;

    struct
    {
        bool operator()(const gsVector<T> & a, const gsVector<T> & b) const
        {
            GISMO_ASSERT(a.size()==b.size(),"Sizes must agree!");
            return std::lexicographical_compare(  a.begin(), a.end(),b.begin(), b.end());
        };
    }
    lexcomp;

    T comp_tol = 1e-5;
    auto comp = [&comp_tol](const gsVector<T> & a, const gsVector<T> & b)
    {
        return (a-b).norm() < comp_tol;
    };

    gsVector<T> x = gsVector<T>::LinSpaced(10,-0.5 * M_PI,0.5 * M_PI);

    for (index_t k = 0; k!=Cs.cols(); k++)
    {
        init.setZero();
        value.setZero();

        gsMatrix<T> C = Cs.col(k);
                    C.resize(3,3);
        gsMatrix<T> N = Ns.col(k);
                    N.resize(3,1);

        std::vector<gsVector<T>> results; results.reserve(x.size());
        for (index_t t=0; t!=x.size(); t++)
        {
            init.at(0) = x.at(t);

            objective<T> obj(C,N);
            iter = obj.newtonRaphson(value,init,false,1e-12,1000);
            GISMO_ASSERT(iter!=-1,"Newton iterations did not converge");

            theta = init.at(0);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            gamma = - ( n1_vec.transpose() * N ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

            gsVector<> res(2);
            res.at(0) = std::fmod(theta,M_PI);
            res.at(1) = gamma;

            if (res(1) > 0)
            {
                results.push_back(res);
            }
        }

        std::sort(results.begin(), results.end(),lexcomp);

        auto last = std::unique(results.begin(), results.end(),comp);
        results.erase(last, results.end());

        GISMO_ASSERT(results.size()>=1,"No suitable theta found");

        theta = results.at(0)(0);

        result(0,k) = theta;
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & dC)
{
    GISMO_ASSERT(C .rows()==3,"Number of rows of C must be 3, but is "<<C.rows());
    GISMO_ASSERT(C .cols()==3,"Number of cols of C must be 3, but is "<<C.cols());
    GISMO_ASSERT(S .rows()==3,"Number of rows of S must be 3, but is "<<S.rows());
    GISMO_ASSERT(S .cols()==1,"Number of cols of S must be 1, but is "<<S.cols());
    GISMO_ASSERT(dC.cols()==3,"Number of cols of dC must be 3, but is "<<dC.cols());
    GISMO_ASSERT(dC.rows()==9,"Number of rows of dC must be 1, but is "<<dC.rows());

    gsMatrix<T> result = C;

    T n1 = math::cos(theta);
    T n2 = math::sin(theta);
    T m1 = -math::sin(theta);
    T m2 = math::cos(theta);
    gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
    gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
    gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

    // Derivative of C to C, dC/dC, times N1
    gsMatrix<T> dCN1(3,3);
    for (index_t d = 0; d!=dC.cols(); d++)
    {
        dCN1.col(d) = dC.reshapeCol(d,3,3) * n1_vec;
    }

    // Gamma
    T gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

    // dGamma/dE
    gsMatrix<T> I(3,3); I.setIdentity();
    gsMatrix<T> dgammadE = - ( n1_vec.transpose() * C * I) / ( n1_vec.transpose() * C * n1_vec ).value()
                           + ( n1_vec.transpose() * S ).value() * ( n1_vec.transpose() * dCN1 * I) / (std::pow( (n1_vec.transpose() * C * n1_vec).value(),2));

    dgammadE.transposeInPlace();

    // dGamma/dT
    T tmp = ( n2_vec.transpose() * C * n1_vec ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
    T dgammadT = - 2 * gamma * tmp;

    // dF/dE
    // gsMatrix<T> dfdE = ( n2_vec.transpose() - tmp * n1_vec.transpose() ) * C * I;
    gsMatrix<T> dfdE = n2_vec.transpose() * C * I + dgammadE.transpose() * (n2_vec.transpose() * C * n1_vec).value();
    dfdE.transposeInPlace();

    // dF/dT
    T dfdT = (n4_vec.transpose() * ( S + C * gamma * n1_vec )).value() +
            2 * gamma * ( n2_vec.transpose() * C * n2_vec - math::pow(( n1_vec.transpose() * C * n2_vec ).value(),2) / (n1_vec.transpose() * C * n1_vec).value() );

    // dT/dE
    gsMatrix<T> dTdE = - dfdE / dfdT;

    // dF
    gsMatrix<T,1,1> tmp2 = (n1_vec.transpose() * C * n2_vec);
    T df =  (n4_vec.transpose() * (S + C * gamma * n1_vec)).value()
            + 2 * gamma * ( (n2_vec.transpose() * C * n2_vec).value() - math::pow(tmp2(0,0),2)
                / (n1_vec.transpose() * C * n1_vec)
            );

    result += C * ( n1_vec * dgammadE.transpose() + dgammadT * n1_vec * dTdE.transpose() - 2*gamma * n2_vec * dTdE.transpose() );

    return result;
}


template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_S(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S)
{
    gsMatrix<T> result = S;

    T n1 = math::cos(theta);
    T n2 = math::sin(theta);
    gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;

    // Gamma
    T gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

    gsMatrix<T> Ew = gamma * n1_vec;
    result += C * Ew;
    return result;
}


// template <short_t dim, class T, bool linear >
// gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval_theta(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
// {
//     gsMatrix<T> result(1,u.cols());
//     gsMatrix<T> E, Evec;
//     gsVector<T> value(1),init(1);
//     gsVector<T,3> n1_vec, n2_vec;
//     T n1, n2, m1, m2;

//     index_t iter;

//     T theta;
//     T gamma;

//     struct
//     {
//         bool operator()(const gsVector<T> & a, const gsVector<T> & b) const
//         {
//             gsDebugVar(a);
//             gsDebugVar(b);
//             GISMO_ASSERT(a.size()==b.size(),"Sizes must agree!");
//             return std::lexicographical_compare(  a.begin(), a.end(),b.begin(), b.end());
//         };
//     }
//     lexcomp;

//     real_t comp_tol = 1e-5;
//     auto comp = [&comp_tol](const gsVector<T> & a, const gsVector<T> & b)
//     {
//         return std::pow((a(0)-b(0)),2) + std::pow((a(1)-b(1)),2) < comp_tol;
//     };

//     gsMatrix<T> Cs = m_materialMat->eval3D_matrix(patch,u,z,MaterialOutput::MatrixA);
//     gsMatrix<T> Ns = m_materialMat->eval3D_vector(patch,u,z,MaterialOutput::VectorN);

//     gsVector<T> x = gsVector<T>::LinSpaced(10,-0.5 * M_PI,0.5 * M_PI);

//     for (index_t k = 0; k!=u.cols(); k++)
//     {
//         init.setZero();
//         value.setZero();

//         gsAsMatrix<T> C = Cs.reshapeCol(k,3,3);
//         gsAsMatrix<T> N = Ns.reshapeCol(k,3,1);

//         std::vector<gsVector<T>> results; results.reserve(x.size());
//         for (index_t t=0; t!=x.size(); t++)
//         {
//             init.at(0) = x.at(t);

//             objective<T> obj(C,N);
//             iter = obj.newtonRaphson(value,init,false,1e-12,1000);
//             GISMO_ASSERT(iter!=-1,"Newton iterations did not converge");

//             theta = init.at(0);
//             n1 = math::cos(theta);
//             n2 = math::sin(theta);
//             m1 = -math::sin(theta);
//             m2 = math::cos(theta);

//             n1_vec<<n1*n1, n2*n2, 2*n1*n2;
//             n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

//             gamma = - ( n1_vec.transpose() * N ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

//             gsVector<> res(2);
//             res.at(0) = std::fmod(theta,M_PI);
//             res.at(1) = gamma;

//             if (res(1) > 0)
//             {
//                 results.push_back(res);
//             }
//         }

//         for (auto it = results.begin(); it!=results.end(); it++)
//             gsDebugVar(*it);

//         std::sort(results.begin(), results.end(),lexcomp);

//         auto last = std::unique(results.begin(), results.end(),comp);
//         results.erase(last, results.end());

//         GISMO_ASSERT(results.size()>=1,"No suitable theta found");

//         theta = results.at(0)(0);

//         result(0,k) = theta;
//     }

// }

} // end namespace
