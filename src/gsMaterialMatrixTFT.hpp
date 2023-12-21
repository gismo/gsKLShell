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

#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsCore/gsFunction.h>

using namespace gismo;

template <typename T>
T mod(const T x, const T y)
{
    return x - math::floor(x/y)*y;
}

template <typename T> class objective;

namespace gismo
{

template <short_t dim, class T, bool linear >
std::ostream & gsMaterialMatrixTFT<dim,T,linear>::print(std::ostream &os) const
{
    os  <<"---------------------------------------------------------------------\n"
        <<"---------------------Tension-Field Theory Material-------------------\n"
        <<"---------------------------------------------------------------------\n\n";

    m_materialMat->print(os);

    os  <<"---------------------------------------------------------------------\n\n";
    return os;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // GISMO_ASSERT(out==MaterialOutput::MatrixA,"Tension Field Theory only works for membrane models, hence only outputs the A matrix");
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    // return this->_eval3D_matrix_impl<true>(patch,u,z);
    return this->_eval3D_matrix_impl<linear>(patch,u,z);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return this->eval3D_stress(patch,u,z,out);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);
    gsMatrix<T> result = m_materialMat->eval3D_pstress(patch,u,z,out);
    index_t colIdx;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j*u.cols()+k;
            if (TF(0,colIdx) == 1 || TF(0,colIdx) == 0) // taut & wrinkled
            {
                // do nothing
            }
            else if (TF(0,colIdx) == -1) // slack
            {
                // Set to zero
                // result.col(colIdx).setZero();
                result.col(colIdx) *= m_options.getReal("SlackMultiplier");
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);
    gsMatrix<T> result = m_materialMat->eval3D_pstrain(patch,u,z);
    index_t colIdx;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j*u.cols()+k;
            if (TF(0,colIdx) == 1 || TF(0,colIdx) == 0) // taut & wrinkled
            {
                // do nothing
            }
            else if (TF(0,colIdx) == -1) // slack
            {
                // Set to zero
                // result.col(colIdx).setZero();
                result.col(colIdx) *= m_options.getReal("SlackMultiplier");
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    gsMatrix<T> result = m_materialMat->eval3D_strain(patch,u,z);
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);
    index_t colIdx;
    T theta;
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
                // result.col(colIdx).setZero();
                result.col(colIdx) *= m_options.getReal("SlackMultiplier");
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                gsMatrix<T> N = m_materialMat->eval3D_stress(patch,u.col(k),z(j,k),MaterialOutput::VectorN);
                gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k));
                gsMatrix<T> thetas = eval_theta(C,N,E);
                theta = thetas(0,0);
                result.col(colIdx) = this->_compute_E(theta,C.reshape(3,3),N,E);
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    gsMatrix<T> result = m_materialMat->eval3D_stress(patch,u,z,MaterialOutput::VectorN);
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);
    index_t colIdx;
    T theta;
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
                // result.col(colIdx).setZero();
                result.col(colIdx) *= m_options.getReal("SlackMultiplier");
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                gsAsMatrix<T,Dynamic,Dynamic> N = result.reshapeCol(colIdx,3,1);
                gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k));

                gsMatrix<T> thetas = eval_theta(C,N,E);
                theta = thetas(0,0);
                result.col(colIdx) = this->_compute_S(theta,C.reshape(3,3),N);
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    gsMatrix<T> result = m_materialMat->eval3D_CauchyStress(patch,u,z,MaterialOutput::CauchyStressN);
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);
    index_t colIdx;
    T theta;

    this->_computePoints(patch,u);
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
                // result.col(colIdx).setZero();
                result.col(colIdx) *= m_options.getReal("SlackMultiplier");
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k)); // on point i, on height z(0,j)
                gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                gsMatrix<T> N = m_materialMat->eval3D_stress(patch,u.col(k),z(j,k),MaterialOutput::VectorN);;
                gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k));
                gsMatrix<T> thetas = eval_theta(C,N,E);
                theta = thetas(0,0);
                gsMatrix<T> S = this->_compute_S(theta,C.reshape(3,3),N);
                T detF = math::sqrt(m_data.mine().m_J0_sq)*1.0;
                result.col(colIdx) = S / math::sqrt(detF);
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

// template <short_t dim, class T, bool linear >
// gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
// {
//     // Input: u in-plane points
//     //        z matrix with, per point, a column with z integration points
//     // Output: (n=u.cols(), m=z.rows())
//     //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

//     gsMatrix<T> result = m_materialMat->eval3D_CauchyPStress(patch,u,z,MaterialOutput::PCauchyStressN);
//     gsMatrix<T> TF = this->_compute_TF(patch,u,z);
//     index_t colIdx;
//     T theta;
//     std::pair<gsVector<T>,gsMatrix<T>> res;

//     this->_computePoints(patch,u);
//     for (index_t k=0; k!=u.cols(); k++)
//     {
//         for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
//         {
//             colIdx = j*u.cols()+k;
//             if (TF(0,colIdx) == 1) // taut
//             {
//                 // do nothing
//             }
//             else if (TF(0,colIdx) == -1) // slack
//             {
//                 // Set to zero
//                 // result.col(colIdx).setZero();
//                 result.col(colIdx) *= m_options.getReal("SlackMultiplier");
//             }
//             else if (TF(0,colIdx) == 0) // wrinkled
//             {
//                 this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k)); // on point i, on height z(0,j)
//                 gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
//                 gsMatrix<T> N = m_materialMat->eval3D_stress(patch,u.col(k),z(j,k),MaterialOutput::VectorN);;
//                 gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k),MaterialOutput::StrainN);
//                 gsMatrix<T> thetas = eval_theta(C,N,E);
//                 theta = thetas(0,0);
//                 gsMatrix<T> S = this->_compute_S(theta,C.reshape(3,3),N);
//                 res = this->_evalPStress(S);

//                 T detF = math::sqrt(m_data.mine().m_J0_sq)*1.0;
//                 result.col(colIdx) = res.first / math::sqrt(detF);
//             }
//             else
//                 GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
//         }
//     }
//     return result;
// }

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    return this->_compute_TF(patch,u,z);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_theta(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    gsMatrix<T> result(1,u.cols());
    result.setZero();
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);
    index_t colIdx;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j*u.cols()+k;
            if (TF(0,colIdx) == 1 || TF(0,colIdx) == -1) // taut
            {
                result(0,colIdx) = NAN;
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                gsMatrix<T> N = m_materialMat->eval3D_stress(patch,u.col(k),z(j,k),MaterialOutput::VectorN);;
                gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k));

                gsMatrix<T> thetas = eval_theta(C,N,E);
                result(0,colIdx) = thetas(0,0);
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval3D_gamma(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    gsMatrix<T> result(1,u.cols());
    result.setZero();
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);
    index_t colIdx;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j*u.cols()+k;
            if (TF(0,colIdx) == 1 || TF(0,colIdx) == -1) // taut
            {
                result(0,colIdx) = NAN;
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                gsMatrix<T> N = m_materialMat->eval3D_stress(patch,u.col(k),z(j,k),MaterialOutput::VectorN);;
                gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k));

                gsMatrix<T> thetas = eval_theta(C,N,E);
                result(0,colIdx) = _compute_gamma(thetas(0,0),C.reshape(3,3),N.reshape(3,1),E.reshape(3,1));
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }
    return result;
}

/// Provides an implementation of eval3D_matrix for \a gsMaterialMatrixLinear
template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if<_linear, gsMatrix<T> >::type
gsMaterialMatrixTFT<dim,T,linear>::_eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out) const
{
    gsMatrix<T> result = m_materialMat->eval3D_matrix(patch,u,z,MaterialOutput::MatrixA);
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);

    index_t colIdx;
    T theta;
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
                // result.col(colIdx).setZero();
                result.col(colIdx) *= m_options.getReal("SlackMultiplier");
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = result.col(colIdx);
                gsMatrix<T> N = m_materialMat->eval3D_stress(patch,u.col(k),z(j,k),MaterialOutput::VectorN);
                gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k));
                gsMatrix<T> thetas = eval_theta(C,N,E);

                theta = thetas(0,0);
                result.col(colIdx) = this->_compute_C(theta,C.reshape(3,3),N.reshape(3,1)).reshape(9,1);
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
    gsMatrix<T> result = m_materialMat->eval3D_matrix(patch,u,z,out);
    gsMatrix<T> TF = this->_compute_TF(patch,u,z);

    index_t colIdx;
    T theta;
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
                // result.col(colIdx).setZero();
                result.col(colIdx) *= m_options.getReal("SlackMultiplier");
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                gsMatrix<T> C = result.col(colIdx);
                gsMatrix<T> N = m_materialMat->eval3D_stress(patch,u.col(k),z(j,k),MaterialOutput::VectorN);
                gsMatrix<T> E = m_materialMat->eval3D_strain(patch,u.col(k),z(j,k));
                gsMatrix<T> thetas = eval_theta(C,N,E);

                gsMatrix<T> dC = m_materialMat->eval3D_dmatrix(patch,u.col(k),z(j,k),MaterialOutput::MatrixA);
                dC *= 2; // since eval3D_dmatrix gives dC/dC and dC/dE = 2*dC/dC
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
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::eval_theta(const gsMatrix<T> & Cs, const gsMatrix<T> & Ns, const gsMatrix<T> & Es) const
{
    GISMO_ASSERT(Cs.cols()==Ns.cols(),"Number of C matrices and N vectors is different");
    GISMO_ASSERT(Cs.cols()==Es.cols(),"Number of C matrices and E vectors is different");
    gsMatrix<T> result(1,Cs.cols());
    gsVector<T,3> n1_vec, n2_vec, n4_vec;

    T theta = 0;

    for (index_t k = 0; k!=Cs.cols(); k++)
    {
        gsMatrix<T> C = Cs.col(k);
                    C.resize(3,3);
        gsMatrix<T> N = Ns.col(k);
                    N.resize(3,1);
        gsMatrix<T> E = Es.col(k);
                    E.resize(3,1);

        gsVector<T> theta_interval = _theta_interval(C, N, E);
        objective<T> obj(C,N);

        gsVector<T> zeros(1); zeros<<0;
        gsVector<T> arg(1); //arg<<theta;
        T f;
        bool converged = false;
        if (theta_interval(0)!=0 && theta_interval(1)!=0)
        {
            arg<<theta_interval(0);
            converged = obj.findRootBrent(f,theta_interval(0),theta_interval(1),theta,1e-18);
            if (!converged)
            {
                converged = obj.findRootBisection(f,theta_interval(0),theta_interval(1),theta,1e-6);
                if (!converged)
                {
                    obj.newtonRaphson(zeros,arg,false,1e-12,100,1);
                    theta = arg(0);
                    f = obj.eval(theta);
                }
            }
        }

        if (!_check_theta_full(theta,C,N,E) || (theta_interval(0)==0 && theta_interval(1)==0))
        {
            // Take the whole interval and evaluate 50 points
            gsVector<T> theta_vec = gsVector<T>::LinSpaced(51,0,M_PI);
            gsVector<T> f_vec(theta_vec.size());
            // Compute f for each point
            T t;
            for (index_t i=0; i!=theta_vec.size(); i++)
            {
                t = theta_vec(i);
                f_vec(i) = obj.eval(t);
            }
            // For each interval, store the intervals where the sign of f changes
            std::vector<std::pair<T,T>> intervals; intervals.reserve(theta_vec.size());
            for (index_t i=0; i!=f_vec.size()-1; i++)
                if (f_vec(i)*f_vec(i+1) < 0)
                    intervals.push_back(std::make_pair(theta_vec(i),theta_vec(i+1)));

            // Compute all intervals
            std::vector<T> thetas(intervals.size());
            for (size_t i = 0; i!=intervals.size(); i++)
            {
                arg<<intervals[i].first;
                converged = obj.findRootBrent(f,intervals[i].first,intervals[i].second,theta,1e-18);
                if (!converged)
                {
                    converged = obj.findRootBisection(f,intervals[i].first,intervals[i].second,theta,1e-6);
                    if (!converged)
                    {
                        obj.newtonRaphson(zeros,arg,false,1e-12,100,1);
                        theta = arg(0);
                        f = obj.eval(theta);
                    }
            }
                thetas[i] = theta;
            }

            if ((f_vec(0)*f_vec(f_vec.size()-1) <= 0))
                thetas.push_back(0.0);

            std::vector<bool> full_check(thetas.size());
            for (size_t i = 0; i!=thetas.size(); i++)
                full_check[i] = _check_theta_full(thetas[i],C,N,E);

            index_t full_check_sum = std::accumulate(full_check.begin(),full_check.end(),0);
            if (full_check_sum!=0)
            {
                std::vector<bool>::iterator res = std::find_if(full_check.begin(), full_check.end(), [](bool i) { return i;});
                index_t i = std::distance(full_check.begin(),res);
                theta = thetas[i];
            }
            else
            {
                gsWarn<<"Everything failed??\n";
                // gsDebugVar(E);
                // gsDebugVar(N);
                // gsDebugVar(C);
                // gsVector<T> x = gsVector<T>::LinSpaced(1000,-M_PI,M_PI);
                // for (index_t XX=0; XX!=x.size(); XX++)
                // {
                //     T t = x(XX);
                //     f = obj.eval(t);

                //     n1 = math::cos(t);
                //     n2 = math::sin(t);
                //     m1 = -math::sin(t);
                //     m2 = math::cos(t);

                //     n1_vec<<n1*n1, n2*n2, 2*n1*n2;
                //     n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
                //     gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

                //     gamma = - ( n1_vec.transpose() * N ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

                //     gsInfo<<x(XX)<<","<<f<<","<<gamma<<","<<(n1_vec.transpose() * N).value()<<","<<(n4_vec.transpose() * Es.col(k)).value()<<"\n";
                // }
                theta = 0;
            }
        }
        result(0,k) = theta;
    }
    return result;
}

template <short_t dim, class T, bool linear >
T gsMaterialMatrixTFT<dim,T,linear>::_compute_gamma(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const
{
    T n1 = math::cos(theta);
    T n2 = math::sin(theta);
    T m1 = -math::sin(theta);
    T m2 = math::cos(theta);
    gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
    gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
    gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;
    
    return - ( n1_vec.transpose() * N ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
}

template <short_t dim, class T, bool linear >
bool gsMaterialMatrixTFT<dim,T,linear>::_check_theta_full(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const
{
    T n1 = math::cos(theta);
    T n2 = math::sin(theta);
    T m1 = -math::sin(theta);
    T m2 = math::cos(theta);
    gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
    gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;
    
    bool check = true;
    check &= ((n1_vec.transpose() * N).value() < 0); // Li et al eq. 63 / eq. 55
    check &= ((n4_vec.transpose() * E).value() > 0); // Li et al eq. 58
    return check;
}

template <short_t dim, class T, bool linear >
bool gsMaterialMatrixTFT<dim,T,linear>::_check_theta_gamma(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const
{
    T n1 = math::cos(theta);
    T n2 = math::sin(theta);
    gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
    return ((n1_vec.transpose() * N).value() < 0); // Li et al eq. 63 / eq. 55
}

template <short_t dim, class T, bool linear >
gsVector<T> gsMaterialMatrixTFT<dim,T,linear>::_theta_interval(const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const
{
    GISMO_ASSERT(C.rows()==C.cols(),"C must be a 3x3 square matrix");
    GISMO_ASSERT(C.rows()==3,"C must be a 3x3 square matrix");
    GISMO_ASSERT(N.rows()==3,"N must be a 3x1 vector");
    GISMO_ASSERT(E.rows()==3,"E must be a 3x1 vector");
    gsVector<T> result(2);
 
    objective<T> obj(C,N);

    // See Lu et al., Finite Element Analysis of Membrane Wrinkling, 2001, International Journal for numerical methods in engineering
    T E11 = E(0);
    T E22 = E(1);
    T E12 = 0.5 * E(2);
    T R_E = math::sqrt( math::pow((E11-E22)/2.,2) + E12*E12 );
    if (R_E==0) {R_E = 1; gsDebugVar("??");};
    T sin_E0 = E12/R_E;
    T cos_E0 = (E22-E11)/(2*R_E);
    std::complex<T> sin_Esqrt = E12*E12-E11*E22;
    std::complex<T> sin_E1_C = -math::sqrt(sin_Esqrt)/R_E; 
    T sin_E1 = sin_E1_C.real();
    T cos_E1 = -(E11+E22)/(2*R_E);
    std::complex<T> sin_E2_C = math::sqrt(sin_Esqrt)/R_E; 
    T sin_E2 = sin_E2_C.real();
    T cos_E2 = -(E11+E22)/(2*R_E);

    T theta_0E = math::atan2(sin_E0,cos_E0);
    T theta_1E = math::atan2(sin_E1,cos_E1);
    T theta_2E = math::atan2(sin_E2,cos_E2);

    T pi = 4*math::atan(1);
    T N11 = N(0);
    T N22 = N(1);
    T N12 = N(2);
    T R_N = math::sqrt( math::pow((N11-N22)/2.,2) + N12*N12 );
    if (R_N==0) R_N = 1;
    T sin_N0 = -N12/R_N;
    T cos_N0 = (N11-N22)/(R_N);
    std::complex<T> sin_Nsqrt = N12*N12-N11*N22;
    // if (sin_Nsqrt < 0) gsDebugVar("Oops");
    std::complex<T> sin_N1_C = math::sqrt(sin_Nsqrt)/R_N; 
    T sin_N1 = sin_N1_C.real();
    T cos_N1 = -(N11+N22)/(2*R_N);
    std::complex<T> sin_N2_C = -math::sqrt(sin_Nsqrt)/R_N; 
    T sin_N2 = sin_N2_C.real();
    T cos_N2 = -(N11+N22)/(2*R_N);

    T theta_0N = math::atan2(sin_N0,cos_N0);
    T theta_1N = math::atan2(sin_N1,cos_N1);
    T theta_2N = math::atan2(sin_N2,cos_N2);

    T theta_1 = mod((theta_1N - theta_0N + theta_0E + pi),(2*pi));
    T theta_2 = mod((theta_2N - theta_0N + theta_0E + pi),(2*pi));
    if (math::isnan(theta_1) || math::isnan(theta_2) || R_N==1 || R_E == 1)
    {
        gsDebugVar(E);
        gsDebugVar(C);
        gsDebugVar(N);


        gsDebugVar(R_E);

        gsDebugVar(sin_Esqrt);

        gsDebugVar(sin_E0);
        gsDebugVar(sin_E1);
        gsDebugVar(sin_E2);

        gsDebugVar(cos_E0);
        gsDebugVar(cos_E1);
        gsDebugVar(cos_E2);

        gsDebugVar(R_N);

        gsDebugVar(sin_Nsqrt);

        gsDebugVar(sin_N0);
        gsDebugVar(sin_N1);
        gsDebugVar(sin_N2);

        gsDebugVar(cos_N0);
        gsDebugVar(cos_N1);
        gsDebugVar(cos_N2);


        gsDebugVar(theta_1N);
        gsDebugVar(theta_0N);
        gsDebugVar(theta_0E);
    }

    T Q_E_min = theta_1E - theta_0E;
    T Q_E_max = theta_2E - theta_0E;
    T Q_S_min, Q_S_max;

    T min, max;
    if (theta_1 < theta_2)
    {
        Q_S_min = theta_1 - pi - theta_0E;
        Q_S_max = theta_2 - pi - theta_0E;

        // set difference:
        min = math::max(Q_E_min,Q_S_min);
        max = math::min(Q_E_max,Q_S_max);
        min /=2;
        max /=2;
        // gsDebugVar("One interval");
    }
    else if (theta_1 > theta_2)
    {
        // try first interval
        Q_S_min = theta_1 - pi - theta_0E;
        Q_S_max = pi - theta_0E;

        // set difference:
        min = math::max(Q_E_min,Q_S_min);
        max = math::min(Q_E_max,Q_S_max);
        min /=2;
        max /=2;
        if (obj.eval(min)*obj.eval(max)>0)
        {
            // second interval
            Q_S_min = - pi - theta_0E;
            Q_S_max = theta_2 - pi - theta_0E;
            // set difference:
            min = math::max(Q_E_min,Q_S_min);
            max = math::min(Q_E_max,Q_S_max);
            min /=2;
            max /=2;
        //     // if (obj.eval(min)*obj.eval(max)>0)
        //     // {
        //     //     Q_S_min = theta_1 - pi - theta_0E;
        //     //     Q_S_max = pi - theta_0E;
        //     //     gsDebugVar(obj.eval(math::max(Q_E_min,Q_S_min)/2));
        //     //     gsDebugVar(obj.eval(math::min(Q_E_max,Q_S_max)/2));
        //     //     gsDebugVar(obj.eval(math::max(Q_E_min,Q_S_min)/2)*obj.eval(math::min(Q_E_max,Q_S_max)/2));
        //     //     Q_S_min = - pi - theta_0E;
        //     //     Q_S_max = theta_2 - pi - theta_0E;
        //     //     gsDebugVar(obj.eval(math::max(Q_E_min,Q_S_min)/2));
        //     //     gsDebugVar(obj.eval(math::min(Q_E_max,Q_S_max)/2));
        //     //     gsDebugVar(obj.eval(math::max(Q_E_min,Q_S_min)/2)*obj.eval(math::min(Q_E_max,Q_S_max)/2));

        //     //     GISMO_ERROR("Root is not found??");
        //     // }
        }
        // gsDebugVar("Two intervals");
    }
    else
    {
        min = 0;
        max = 0;
        // gsWarn<<"Interval not specified\n";
        // GISMO_ERROR("Don't know what to do!");
    }
    result(0) = min;
    result(1) = max;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_TF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    if (m_options.getSwitch("Explicit"))
    {
        function_ptr tmp = m_materialMat->getDeformed();
        m_materialMat->setDeformed(m_defpatches0);
        gsMatrix<T> TF = m_materialMat->eval3D_tensionfield(patch,u,z,MaterialOutput::TensionField);
        m_materialMat->setDeformed(tmp);
        return TF;
    }
    else
        return m_materialMat->eval3D_tensionfield(patch,u,z,MaterialOutput::TensionField);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_TF(const index_t patch, const gsVector<T>& u, const T & z) const
{
    gsMatrix<T> zmat(1,1);
    zmat<<z;
    return _compute_TF(patch,u,zmat);
}


template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_E(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & E) const
{
    gsMatrix<T> result = E;

    T n1 = math::cos(theta);
    T n2 = math::sin(theta);
    gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;

    // Gamma
    T gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

    gsMatrix<T> Ew = gamma * n1_vec;
    result += Ew;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_S(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S) const
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

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S) const
{
    GISMO_ASSERT(C .rows()==3,"Number of rows of C must be 3, but is "<<C.rows());
    GISMO_ASSERT(C .cols()==3,"Number of cols of C must be 3, but is "<<C.cols());
    GISMO_ASSERT(S .rows()==3,"Number of rows of S must be 3, but is "<<S.rows());
    GISMO_ASSERT(S .cols()==1,"Number of cols of S must be 1, but is "<<S.cols());

    gsMatrix<T> result = C;

    T n1 = math::cos(theta);
    T n2 = math::sin(theta);
    T m1 = -math::sin(theta);
    T m2 = math::cos(theta);
    gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
    gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
    gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

    // Gamma
    T gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

    // dGamma/dE
    gsMatrix<T> I(3,3); I.setIdentity();
    gsMatrix<T> dgammadE = - ( n1_vec.transpose() * C) / ( n1_vec.transpose() * C * n1_vec ).value();
    dgammadE.transposeInPlace();

    // dGamma/dT
    T tmp = ( n2_vec.transpose() * C * n1_vec ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
    T dgammadT = - 2 * gamma * tmp;

    // dF/dE
    // gsMatrix<T> dfdE = ( n2_vec.transpose() - tmp * n1_vec.transpose() ) * C;
    gsMatrix<T> dfdE = n2_vec.transpose() * C - dgammadE.transpose() * (n2_vec.transpose() * C * n1_vec).value();
    dfdE.transposeInPlace();

    // dF/dT
    // NAKASHINO
    // T dfdT = (n4_vec.transpose() * ( S + C * gamma * n1_vec )).value() +
    //         2 * gamma * ( n2_vec.transpose() * C * n2_vec - math::pow(( n1_vec.transpose() * C * n2_vec ).value(),2) / ( n1_vec.transpose() * C * n1_vec ).value() );
    T dfdT = (n4_vec.transpose() * S).value() + dgammadT * (n2_vec.transpose() * C * n1_vec).value() +
                gamma * ( (n4_vec.transpose() * C * n1_vec).value() + 2*(n2_vec.transpose() * C * n2_vec).value());


    // dT/dE
    gsMatrix<T> dTdE = - dfdE / dfdT;

    result += C * ( n1_vec * dgammadE.transpose() + dgammadT * n1_vec * dTdE.transpose() + 2*gamma * n2_vec * dTdE.transpose() );
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFT<dim,T,linear>::_compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & dC) const
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
    gsMatrix<T> dgammadE = - ( n1_vec.transpose() * C + gamma * ( n1_vec.transpose() * dCN1) ) / ( n1_vec.transpose() * C * n1_vec ).value();
    // dgammadE.transposeInPlace();

    // dGamma/dT
    T tmp = ( n2_vec.transpose() * C * n1_vec ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
    T dgammadT = - 2 * gamma * tmp;

    // dF/dE
    // gsMatrix<T> dfdE = ( n2_vec.transpose() - tmp * n1_vec.transpose() ) * C;
    gsMatrix<T> dfdE = n2_vec.transpose() * C + dgammadE * (n2_vec.transpose() * C * n1_vec).value() + gamma * ( n2_vec.transpose() * dCN1);
    // dfdE.transposeInPlace();

    // dF/dT
    // NAKASHINO
    // T dfdT = (n4_vec.transpose() * ( S + C * gamma * n1_vec )).value() +
    //         2 * gamma * ( n2_vec.transpose() * C * n2_vec - math::pow(( n1_vec.transpose() * C * n2_vec ).value(),2) / ( n1_vec.transpose() * C * n1_vec ).value() );
    T dfdT = (n4_vec.transpose() * S).value() + dgammadT * (n2_vec.transpose() * C * n1_vec).value() +
                gamma * ( (n4_vec.transpose() * C * n1_vec).value() + 2*(n2_vec.transpose() * C * n2_vec).value());

    // dT/dE
    gsMatrix<T> dTdE = - dfdE / dfdT;

    // d(gamma*n1)dE
    gsMatrix<T> dgamman1dE = ( n1_vec * dgammadE + dgammadT * n1_vec * dTdE + 2*gamma * n2_vec * dTdE );
    result += dgamman1dE.transpose() * C;
    result += gamma * dCN1;
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
//     gsMatrix<T> Ns = m_materialMat->eval3D_stress(patch,u,z,MaterialOutput::VectorN);

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

// namespace internal
// {

// /// @brief get a Linear Material Matrix from XML data
// ///
// /// \ingroup KLShell
// template<short_t d, class T> // implemented for template linear=false
// class gsXml< gsMaterialMatrixTFT<d,T,false> >
// {
// private:
//     gsXml() { }
//     typedef gsMaterialMatrixTFT<d,T,false> Object;

// public:
//     GSXML_COMMON_FUNCTIONS(gsMaterialMatrixTFT<TMPLA3(d,T,false)>);
//     static std::string tag ()  { return "MaterialMatrix"; }
//     static std::string type () { return "Linear" +  to_string(d); }

//     GSXML_GET_POINTER(Object)

//     static void get_into(gsXmlNode * node, Object & obj)
//     {
//         gsMaterialMatrixBase<T> * base = gsXml< gsMaterialMatrixBase<T> >::get(node);
//         obj = Object(base);
//     }

//     static gsXmlNode * put (const Object & obj,
//                             gsXmlTree & data)
//     {
//         gsXml< gsMaterialMatrixBase<T> >::put(node,data);

//         return putMaterialMatrixToXml< Object >( obj,data );
//         // GISMO_NO_IMPLEMENTATION;
//         // return putGeometryToXml(obj,data);
//     }
// };

// }// namespace internal

} // end namespace

template <typename T>
class objective : public gsFunction<T>
{
public:
    objective(const gsMatrix<T> & C, const gsMatrix<T> & S)
    :
    m_C(C),
    m_S(S)
    {
        m_support.resize(1,2);
        m_support<<0,2*3.1415926535;
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

    T eval(const T & theta) const
    {
        gsMatrix<T> t(1,1);
        t<<theta;
        gsMatrix<T> res;
        this->eval_into(t,res);
        return res(0,0);
    }

    short_t domainDim() const
    {
        return 1;
    }

    short_t targetDim() const
    {
        return 1;
    }

    gsMatrix<T> support() const
    {
        return m_support;
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

    void setSupport(const gsMatrix<T> supp)
    {
        m_support = supp;
    }

    /**
     * @brief      Finds a root bisection.
     *
     * @param      f      The function value
     * @param[in]  A      Left end of the interval
     * @param[in]  B      Right end of the interval
     * @param      x      Final x coordinate of root
     * @param[in]  t      Tolerance
     * @param[in]  itmax  Maximum iterations
     *
     * @return     True when converged
     */
    bool findRootBisection(T & f, const T & A, const T & B, T & x, const T & t = 1e-12, const index_t & itmax = 1000)
    {
        T a = A;
        T b = B;
        T c = (b-a)/2;
        index_t it = 0;
        while (math::abs(b-a) > t && it++ < itmax)
        {
            if (this->eval(c)==0.0)
                break;
            else if (this->eval(c)*this->eval(a) < 0)
                b = c;
            else
                a = c;
        }
        x = c;
        if (math::abs(b-a) < t)
            return true;
        else
            return false;
    }

    /**
     * @brief      Finds a root using Brent's method.
     *
     * @param      f      The function value
     * @param[in]  A      Left end of the interval
     * @param[in]  B      Right end of the interval
     * @param      x      Final x coordinate of root
     * @param[in]  t      Tolerance
     * @param[in]  itmax  Maximum iterations
     *
     * @return     True when converged
     *
     * @return     { description_of_the_return_value }
     */
    bool findRootBrent(T & f, const T & a, const T & b, T & x, const T & t = 1e-12, const index_t & itmax = 1000)
    {
        /*
         * Implementation of Brent's method for finding the root \a c of a scalar equation \a *this within the interval [a,b]
         *
         * Adopted from: https://people.math.sc.edu/Burkardt/cpp_src/brent/brent.html
         * The brent.cpp file, function zero(  )
         * https://people.math.sc.edu/Burkardt/cpp_src/brent/brent.cpp
         */
        T c;
        T d;
        T e;
        T fa;
        T fb;
        T fc;
        T m;
        T macheps;
        T p;
        T q;
        T r;
        T s;
        T sa;
        T sb;
        T tol;
//
//  Make local copies of A and B.
//
        sa = a;
        sb = b;
        fa = this->eval( sa );
        fb = this->eval( sb );

        if (fa*fb >=0)
            return false;
        else if (fa > fb)
            std::swap(fa,fb);

        c = sa;
        fc = fa;
        e = sb - sa;
        d = e;

        macheps = std::numeric_limits<T>::epsilon();

        for ( index_t it = 0; it!=itmax; it++ )
        {
            if ( std::fabs ( fc ) < std::fabs ( fb ) )
            {
                sa = sb;
                sb = c;
                c = sa;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            tol = 2.0 * macheps * std::fabs ( sb ) + t;
            m = 0.5 * ( c - sb );

            if ( std::fabs ( m ) <= tol || fb == 0.0 )
            {
                x = sb;
                f = std::fabs ( m );
                return true;
                // break;
            }

            if ( std::fabs ( e ) < tol || std::fabs ( fa ) <= std::fabs ( fb ) )
            {
                e = m;
                d = e;
            }
            else
            {
                s = fb / fa;

                if ( sa == c )
                {
                    p = 2.0 * m * s;
                    q = 1.0 - s;
                }
                else
                {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
                    q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
                }

                if ( 0.0 < p )
                {
                    q = - q;
                }
                else
                {
                    p = - p;
                }

                s = e;
                e = d;

                if ( 2.0 * p < 3.0 * m * q - std::fabs ( tol * q ) &&
                    p < std::fabs ( 0.5 * s * q ) )
                {
                    d = p / q;
                }
                else
                {
                    e = m;
                    d = e;
                }
            }
            sa = sb;
            fa = fb;

            if ( tol < std::fabs ( d ) )
            {
                sb = sb + d;
            }
            else if ( 0.0 < m )
            {
                sb = sb + tol;
            }
            else
            {
                sb = sb - tol;
            }

            fb = this->eval( sb );

            if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
            {
                c = sa;
                fc = fa;
                e = sb - sa;
                d = e;
            }
        }


        GISMO_ERROR("Brent's method did not find a root");
        return false;
    }

private:
    const gsMatrix<T> m_C;
    const gsMatrix<T> m_S;
    mutable gsMatrix<T> m_support;
};
