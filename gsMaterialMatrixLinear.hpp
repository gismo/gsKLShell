/** @file gsMaterialMatrixLinear.hpp

    @brief Provides linear material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

/*
    To Do [updated 16-06-2020]:
    - Make beta (compressible materials) and material parameters universal for all integration points over the thickness. So get them out of the dPsi functions etc and move them into the integration loops as global variables.

*/



#pragma once

#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>

namespace gismo
{

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness
                                        )
                                        :
                                        Base(&mp,nullptr,&thickness,nullptr)
{

}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & YoungsModulus,
                                        const gsFunction<T> & PoissonsRatio
                                        )
                                        :
                                        Base(&mp,nullptr,&thickness,nullptr)
{
    m_pars.resize(2);
    m_pars[0] = const_cast<gsFunction<T> *>(&YoungsModulus);
    m_pars[1] = const_cast<gsFunction<T> *>(&PoissonsRatio);
    _initialize();
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & YoungsModulus,
                                        const gsFunction<T> & PoissonsRatio,
                                        const gsFunction<T> & Density
                                        )
                                        :
                                        Base(&mp,nullptr,&thickness,&Density)
{
    m_pars.resize(2);
    m_pars[0] = const_cast<gsFunction<T> *>(&YoungsModulus);
    m_pars[1] = const_cast<gsFunction<T> *>(&PoissonsRatio);
    _initialize();
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const std::vector<gsFunction<T>*> &pars
                                        )
                                        :
                                        Base(&mp,nullptr,&thickness,nullptr)
{
    GISMO_ASSERT(pars.size()==2,"Two material parameters should be assigned!");
    m_pars = pars;
    _initialize();
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                    const gsFunctionSet<T> & mp,
                                    const gsFunction<T> & thickness,
                                    const std::vector<gsFunction<T>*> &pars,
                                    const gsFunction<T> & density
                                    )
                                    :
                                    Base(&mp,nullptr,&thickness,&density)
{
    GISMO_ASSERT(pars.size()==2,"Two material parameters should be assigned!");
    m_pars = pars;
    _initialize();
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::info() const
{
    gsInfo  <<"---------------------------------------------------------------------\n"
            <<"---------------------Hyperelastic Material Info----------------------\n"
            <<"---------------------------------------------------------------------\n\n";

    gsInfo  <<"Material model: \t";
    gsInfo<<"Saint-Venant Kirchhoff";
    gsInfo<<"\n";

    gsInfo  <<"---------------------------------------------------------------------\n\n";

}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::_defaultOptions()
{
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::_initialize()
{
    // Set default options
    this->_defaultOptions();

    // set flags
    m_map.mine().flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    _computePoints(patch,u,false);
    gsMatrix<T> result(9, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
                // _getMetric(k, z(j, k), false); // on point i, on height z(0,j)
                _getMetric(k, z(j, k) * m_Tmat(0, k), false); // on point i, on height z(0,j)

                gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j*u.cols()+k,3,3);
                /*
                    C = C1111,  C1122,  C1112
                        symm,   C2222,  C2212
                        symm,   symm,   C1212
                */
                C(0,0)          = _Cijkl(0,0,0,0); // C1111
                C(1,1)          = _Cijkl(1,1,1,1); // C2222
                C(2,2)          = _Cijkl(0,1,0,1); // C1212
                C(1,0) = C(0,1) = _Cijkl(0,0,1,1); // C1122
                C(2,0) = C(0,2) = _Cijkl(0,0,0,1); // C1112
                C(2,1) = C(1,2) = _Cijkl(1,1,0,1); // C2212
        }
    }

    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    _computePoints(patch,u,false);
    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // _getMetric(k, z(j, k), false); // on point i, on height z(0,j)
            _getMetric(k, z(j, k) * m_Tmat(0, k), false); // on point i, on height z(0,j)

            result(0, j * u.cols() + k) = _Sij(0, 0, z(j, k), out);
            result(1, j * u.cols() + k) = _Sij(1, 1, z(j, k), out);
            result(2, j * u.cols() + k) = _Sij(0, 1, z(j, k), out);
        }
    }

    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_CauchyVector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    _computePoints(patch,u,true);
    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // _getMetric(k, z(j, k), false); // on point i, on height z(0,j)
            _getMetric(k, z(j, k) * m_Tmat(0, k), true); // on point i, on height z(0,j)

            T detF = math::sqrt(m_J0_sq)*1.0;
            result(0, j * u.cols() + k) = 1 / detF * _Sij(0, 0, z(j, k), out);
            result(1, j * u.cols() + k) = 1 / detF * _Sij(1, 1, z(j, k), out);
            result(2, j * u.cols() + k) = 1 / detF * _Sij(0, 1, z(j, k), out);
        }
    }

    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    _computePoints(patch,u,true);
    gsMatrix<T> result(2, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> S;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // _getMetric(k, z(j, k), true); // on point i, on height z(0,j)
            _getMetric(k, z(j, k) * m_Tmat(0, k), true); // on point i, on height z(0,j)

            S.setZero();
            S(0, 0) = _Sij(0, 0, 0, out);
            S(0, 1) = _Sij(0, 1, 0, out);
            S(1, 0) = _Sij(1, 0, 0, out);
            S(1, 1) = _Sij(1, 1, 0, out);
            S(2, 2) = 0;
            res = _evalPStress(S);
            result.col(j * u.cols() + k) = res.first;
        }
    }

    return result;

}

/*
    Available class members:
        - m_parvals
        - m_metric
        - m_metric_def
*/
template <short_t dim, class T >
T gsMaterialMatrixLinear<dim,T>::_Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ENSURE( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);


    // --------------------------
    // Saint Venant Kirchhoff
    // --------------------------
    T lambda, mu, Cconstant;

    mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    GISMO_ENSURE((1.-2.*m_parvals.at(1)) != 0, "Division by zero in construction of SvK material parameters! (1.-2.*nu) = "<<(1.-2.*m_parvals.at(1))<<"; m_parvals.at(1) = "<<m_parvals.at(1));
    lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1))) ;
    Cconstant = 2*lambda*mu/(lambda+2*mu);

    return Cconstant*m_Acon_ori(i,j)*m_Acon_ori(k,l) + mu*(m_Acon_ori(i,k)*m_Acon_ori(j,l) + m_Acon_ori(i,l)*m_Acon_ori(j,k));
}

template <short_t dim, class T >
T gsMaterialMatrixLinear<dim,T>::_Sij(const index_t i, const index_t j, const T z, enum MaterialOutput out) const
{
    gsMatrix<T> strain;
    GISMO_ENSURE( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
    if      (out == MaterialOutput::VectorN || out == MaterialOutput::PStressN) // To be used with multiplyZ_into
        strain = 0.5*(m_Acov_def - m_Acov_ori);
    else if (out == MaterialOutput::VectorM || out == MaterialOutput::PStressM) // To be used with multiplyZ_into
        strain = (m_Bcov_ori - m_Bcov_def);
    else if (out == MaterialOutput::Generic) // To be used with multiplyLinZ_into or integrateZ_into
        strain = 0.5*(m_Acov_def - m_Acov_ori) + z*(m_Bcov_ori - m_Bcov_def);
    else
        GISMO_ERROR("Output type is not VectorN, PstressN, VectorM, PstressM or Generic!");

    T result =  _Cijkl(i,j,0,0) * strain(0,0) + _Cijkl(i,j,0,1) * strain(0,1)
                + _Cijkl(i,j,1,0) * strain(1,0) + _Cijkl(i,j,1,1) * strain(1,1);

    return result;
}

//--------------------------------------------------------------------------------------------------------------------------------------


template <short_t dim, class T >
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixLinear<dim,T>::_evalPStress(const gsMatrix<T> & S) const
{
    gsVector<T> pstresses;
    gsMatrix<T> pstressvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    pstresses.resize(2,1);    pstresses.setZero();
    pstressvec.resize(3,3);   pstressvec.setZero();

    // Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

    // gsMatrix<T> B(3,3);
    // B.setZero();
    // for (index_t k = 0; k != 2; k++)
    //     for (index_t l = 0; l != 2; l++)
    //         B += S(k,l) * m_gcov_ori.col(k) * m_gcov_ori.col(l).transpose();

    // eigSolver.compute(B);

    // index_t zeroIdx = -1;
    // for (index_t k=0; k!=3; k++)
    //     zeroIdx = eigSolver.eigenvalues()[k] == 0 ? k : zeroIdx;

    // GISMO_ASSERT(zeroIdx!=-1,"No zero found?");

    // index_t count = 0;
    // pstressvec.col(2) = m_gcon_ori.col(2);
    // pstresses(2,0) = S(2,2);

    // for (index_t k=0; k!=3; k++)
    // {
    //     if (k==zeroIdx) continue;
    //     pstressvec.col(count) = eigSolver.eigenvectors().col(k);
    //     pstresses(count,0) = eigSolver.eigenvalues()(k,0);
    //     count++;
    // }

    // BUG: targetDim of pstresses is 2, not 3!

    result.first = pstresses;
    result.second = pstressvec;

    return result;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::_computePStress(const gsMatrix<T> & C) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalPStress(C);
    m_pstress = result.first;
    m_pstressvec = result.second;
}

} // end namespace
