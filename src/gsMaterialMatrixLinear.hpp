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

#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsIO/gsXml.h>
#include <gsKLShell/src/gsMaterialMatrixXml.hpp>

using namespace gismo;

namespace gismo
{

template <short_t dim, class T>
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear()
                                        :
                                        Base()
{
    _initialize();
}

template <short_t dim, class T>
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & thickness
                                        )
                                        :
                                        Base(&mp,&thickness,nullptr)
{
    _initialize();
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & thickness,
                                        const gsFunctionSet<T> & YoungsModulus,
                                        const gsFunctionSet<T> & PoissonsRatio
                                        )
                                        :
                                        gsMaterialMatrixLinear(&mp,&thickness,YoungsModulus,PoissonsRatio,nullptr)
{
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & thickness,
                                        const gsFunctionSet<T> & YoungsModulus,
                                        const gsFunctionSet<T> & PoissonsRatio,
                                        const gsFunctionSet<T> & Density
                                        )
                                        :
                                        gsMaterialMatrixLinear(&mp,&thickness,YoungsModulus,PoissonsRatio,&Density)
{
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> * mp,
                                        const gsFunctionSet<T> * thickness,
                                        const gsFunctionSet<T> & YoungsModulus,
                                        const gsFunctionSet<T> & PoissonsRatio,
                                        const gsFunctionSet<T> * Density
                                        )
                                        :
                                        Base(mp,thickness,Density)
{
    m_pars.resize(2);
    this->setYoungsModulus(YoungsModulus);
    this->setPoissonsRatio(PoissonsRatio);
    _initialize();
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                    const gsFunctionSet<T> & mp,
                                    const gsFunctionSet<T> & thickness,
                                    const std::vector<gsFunctionSet<T>*> &pars
                                    )
                                    :
                                    gsMaterialMatrixLinear(&mp,&thickness,pars,nullptr)
{
}

template <short_t dim, class T>
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                    const gsFunctionSet<T> & mp,
                                    const gsFunctionSet<T> & thickness,
                                    const std::vector<gsFunctionSet<T>*> &pars,
                                    const gsFunctionSet<T> & density
                                    )
                                    :
                                    gsMaterialMatrixLinear(&mp,&thickness,pars,&density)
{
}

template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                    const gsFunctionSet<T> & thickness,
                                    const std::vector<gsFunctionSet<T>*> &pars
                                    )
                                    :
                                    gsMaterialMatrixLinear(nullptr,&thickness,pars,nullptr)
{
}

template <short_t dim, class T>
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                    const gsFunctionSet<T> & thickness,
                                    const std::vector<gsFunctionSet<T>*> &pars,
                                    const gsFunctionSet<T> & density
                                    )
                                    :
                                    gsMaterialMatrixLinear(nullptr,&thickness,pars,&density)
{
}

template <short_t dim, class T>
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                    const gsFunctionSet<T> * mp,
                                    const gsFunctionSet<T> * thickness,
                                    const std::vector<gsFunctionSet<T>*> &pars,
                                    const gsFunctionSet<T> * density
                                    )
                                    :
                                    Base(mp,thickness,density)
{
    GISMO_ASSERT(pars.size()==2,"Two material parameters should be assigned!");
    this->setParameters(pars);
    _initialize();
}

template <short_t dim, class T>
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear( const gsMaterialMatrixLinear<dim,T> & other)
{
    *this =(other);
}

template <short_t dim, class T>
std::ostream & gsMaterialMatrixLinear<dim,T>::print(std::ostream &os) const
{
    os  <<"---------------------------------------------------------------------\n"
        <<"---------------------Elastic Material Info---------------------------\n"
        <<"---------------------------------------------------------------------\n\n";

    os  <<"Material model: \t";
    os  <<"Saint-Venant Kirchhoff";
    os  <<"\n";
    os  <<"---------------------------------------------------------------------\n\n";
    return os;
}

template <short_t dim, class T>
void gsMaterialMatrixLinear<dim,T>::defaultOptions()
{
    Base::defaultOptions();
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T>
void gsMaterialMatrixLinear<dim,T>::_initialize()
{
    // Set default options
    this->defaultOptions();
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->_computePoints(patch,u);
    gsMatrix<T> result(9, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
            m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
                this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k)); // on point i, on height z(0,j)

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

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->_computePoints(patch,u);

    gsMatrix<T> result(9, 1);
    result.setZero();

    // Evaluate material properties on the quadrature point
    for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
        m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,0);

    // Get metric
    this->_getMetricUndeformed(0, z * m_data.mine().m_Tmat(0, 0)); // on point i, on height z(0,j)
    this->_getMetricDeformed(Cmat);

    gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(0,3,3);
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
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_dmatrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return gsMatrix<T>::Zero(27,u.cols()*z.rows());
}

// template <short_t dim, class T>
// gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
// {
//     // gsInfo<<"TO DO: evaluate moments using thickness";
//     // Input: u in-plane points
//     //        z matrix with, per point, a column with z integration points
//     // Output: (n=u.cols(), m=z.cols())
//     //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
//     typedef std::function<gsMatrix<T> (gsMatrix<T> const &)> Stress_t;
//     Stress_t Sfun = [this]( const gsMatrix<T> & E )
//     {
//         gsMatrix<T> S(3,1);
//         S(0,0) = _Sij(0, 0, E);
//         S(1,0) = _Sij(1, 1, E);
//         S(2,0) = _Sij(0, 1, E);
//         return S;
//     };

//     this->_computePoints(patch,u);

//     gsMatrix<T> result(3, u.cols() * z.rows());
//     result.setZero();

//     gsMatrix<T> E;
//     for (index_t k=0; k!=u.cols(); k++)
//     {
//         // Evaluate material properties on the quadrature point
//         for (index_t v=0; v!=m_parmat.rows(); v++)
//             m_parvals.at(v) = m_parmat(v,k);

//         for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
//         {
//             this->_getMetric(k, z(j, k) * m_Tmat(0, k)); // on point i, on height z(0,j)

//             E = _E(z(j, k) * m_Tmat(0, k),out);

//             result(0, j * u.cols() + k) = ;
//             result(1, j * u.cols() + k) = _Sij(1, 1, E);
//             result(2, j * u.cols() + k) = _Sij(0, 1, E);
//         }
//     }

//     return result;
// }

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return eval3D_stress(patch,u,z,out);
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_CauchyVector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return eval3D_CauchyStress(patch,u,z,out);
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_vector_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->_computePoints(patch,u);

    gsMatrix<T> result(3, 1);
    result.setZero();

    gsMatrix<T> E;
    // Evaluate material properties on the quadrature point
    for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
        m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,0);

    // Get metric
    this->_getMetricUndeformed(0, z * m_data.mine().m_Tmat(0, 0)); // on point i, on height z(0,j)
    this->_getMetricDeformed(Cmat);

    // GISMO_ASSERT(   out == MaterialOutput::VectorN ||
    //                 out == MaterialOutput::PStrainN ||
    //                 out == MaterialOutput::PStressN ||
    //                 out == MaterialOutput::TensionField ||
    //                 out == MaterialOutput::StrainN,
    //                 "MaterialOutput is wrong!");
    GISMO_ASSERT(z==0.0,"This only works for z=0");

    m_data.mine().m_Acov_def = m_data.mine().m_Gcov_def.block(0,0,2,2);
    m_data.mine().m_Acon_def = m_data.mine().m_Gcon_def.block(0,0,2,2);
    m_data.mine().m_Bcov_def.setZero();
    m_data.mine().m_Bcon_def.setZero();

    E = _E(z * m_data.mine().m_Tmat(0, 0),out);

    result(0, 0) = _Sij(0, 0, E);
    result(1, 0) = _Sij(1, 1, E);
    result(2, 0) = _Sij(0, 1, E);
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    gsMatrix<T> Smat = eval3D_stress(patch,u,z,out);

    gsMatrix<T> result(2, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> S;
    index_t colIdx;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
            m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j * u.cols() + k;
            S.setZero();
            S(0,0) = Smat(0,colIdx);
            S(1,1) = Smat(1,colIdx);
            S(0,1) = S(1,0) = Smat(2,colIdx);
            S(2,2) = 0;
            res = this->_evalPStress(S);
            result.col(j * u.cols() + k) = res.first;
        }
    }

    return result;

}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_pstressDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    gsMatrix<T> Smat = eval3D_stress(patch,u,z,out);

    gsMatrix<T> result(9, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> S;
    index_t colIdx;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
            m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j * u.cols() + k;
            S.setZero();
            S(0,0) = Smat(0,colIdx);
            S(1,1) = Smat(1,colIdx);
            S(0,1) = S(1,0) = Smat(2,colIdx);
            S(2,2) = 0;
            res = this->_evalPStress(S);
            result.col(j * u.cols() + k) = res.second.reshape(9,1);
        }
    }

    return result;

}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    gsMatrix<T> Smat = eval3D_CauchyStress(patch,u,z,out);

    gsMatrix<T> result(2, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> S;
    index_t colIdx;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
            m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j * u.cols() + k;
            S.setZero();
            S(0,0) = Smat(0,colIdx);
            S(1,1) = Smat(1,colIdx);
            S(0,1) = S(1,0) = Smat(2,colIdx);
            S(2,2) = 0;
            res = this->_evalPStress(S);
            result.col(j * u.cols() + k) = res.first;
        }
    }

    return result;

}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->_computePoints(patch,u);

    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,2,2> E;
    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
            m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k)); // on point i, on height z(0,j)

            // E.setZero();
            E.block(0,0,2,2) = _E(0,out);
            result(0, j * u.cols() + k) = _Sij(0, 0, E);
            result(1, j * u.cols() + k) = _Sij(1, 1, E);
            result(2, j * u.cols() + k) = _Sij(0, 1, E);
        }
    }

    return result;

}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_detF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result(1, u.cols() * z.rows());
    result.setZero();
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k)); // on point i, on height z(0,j)
            result(0,j * u.cols() + k) = math::sqrt(m_data.mine().m_J0_sq*1.0);
        }
    }
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result = eval3D_stress(patch,u,z,out);
    index_t colIdx;
    T detF;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j * u.cols() + k;
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k)); // on point i, on height z(0,j)
            detF = math::sqrt(m_data.mine().m_J0_sq*1.0);
            result.col(colIdx) /= detF;
        }
    }
    return result;
}

/*
    Available class members:
        - m_data.mine().m_parvals
        - m_metric
        - m_metric_def
*/
template <short_t dim, class T>
T gsMaterialMatrixLinear<dim,T>::_Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ENSURE( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);


    // --------------------------
    // Saint Venant Kirchhoff
    // --------------------------
    T lambda, mu, Cconstant;

    mu = m_data.mine().m_parvals.at(0) / (2.*(1. + m_data.mine().m_parvals.at(1)));
    GISMO_ENSURE((1.-2.*m_data.mine().m_parvals.at(1)) != 0, "Division by zero in construction of SvK material parameters! (1.-2.*nu) = "<<(1.-2.*m_data.mine().m_parvals.at(1))<<"; m_data.mine().m_parvals.at(1) = "<<m_data.mine().m_parvals.at(1));
    lambda = m_data.mine().m_parvals.at(0) * m_data.mine().m_parvals.at(1) / ( (1. + m_data.mine().m_parvals.at(1))*(1.-2.*m_data.mine().m_parvals.at(1))) ;
    Cconstant = 2*lambda*mu/(lambda+2*mu);

    return Cconstant*m_data.mine().m_Acon_ori(i,j)*m_data.mine().m_Acon_ori(k,l) + mu*(m_data.mine().m_Acon_ori(i,k)*m_data.mine().m_Acon_ori(j,l) + m_data.mine().m_Acon_ori(i,l)*m_data.mine().m_Acon_ori(j,k));
}

template <short_t dim, class T>
T gsMaterialMatrixLinear<dim,T>::_Sij(const index_t i, const index_t j, const gsMatrix<T> & strain) const
{
    T result =  _Cijkl(i,j,0,0) * strain(0,0) + _Cijkl(i,j,0,1) * strain(0,1)
                + _Cijkl(i,j,1,0) * strain(1,0) + _Cijkl(i,j,1,1) * strain(1,1);
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::S(const gsMatrix<T> & strain) const
{
    gsWarn<<"This is dangerous, since it does not compute material parameters in advance!\n";
    GISMO_ASSERT(strain.rows()==strain.cols() && strain.rows()==2, "Strain tensor must be 2x2!");
    gsMatrix<T> result(2,2);
    for (index_t i = 0; i!=2; i++)
        for (index_t j = 0; j!=2; j++)
            result(i,j) = _Cijkl(i,j,0,0) * strain(0,0) + _Cijkl(i,j,0,1) * strain(0,1)
                        + _Cijkl(i,j,1,0) * strain(1,0) + _Cijkl(i,j,1,1) * strain(1,1);

    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::C(const gsMatrix<T> &) const
{
    gsMatrix<T> result(2,2);
    result(0,0)                 = _Cijkl(0,0,0,0); // C1111
    result(1,1)                 = _Cijkl(1,1,1,1); // C2222
    result(2,2)                 = _Cijkl(0,1,0,1); // C1212
    result(1,0) = result(0,1)   = _Cijkl(0,0,1,1); // C1122
    result(2,0) = result(0,2)   = _Cijkl(0,0,0,1); // C1112
    result(2,1) = result(1,2)   = _Cijkl(1,1,0,1); // C2212
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::_E(const T z, enum MaterialOutput out) const
{
    gsMatrix<T> strain;
    if      (   out == MaterialOutput::VectorN          ||
                out == MaterialOutput::CauchyVectorN    ||
                out == MaterialOutput::StrainN          ||
                out == MaterialOutput::StressN          ||
                out == MaterialOutput::CauchyStressN    ||
                out == MaterialOutput::PStrainN         ||
                out == MaterialOutput::PStressN         ||
                out == MaterialOutput::PCauchyStressN   ||
                out == MaterialOutput::TensionField     ||
                out == MaterialOutput::StrainN              ) // To be used with multiplyZ_into
        strain = 0.5*(m_data.mine().m_Acov_def - m_data.mine().m_Acov_ori);
    else if (   out == MaterialOutput::VectorM          ||
                out == MaterialOutput::CauchyVectorM    ||
                out == MaterialOutput::StressM          ||
                out == MaterialOutput::CauchyStressM    ||
                out == MaterialOutput::StrainM          ||
                out == MaterialOutput::PStrainM         ||
                out == MaterialOutput::PStressM         ||
                out == MaterialOutput::PCauchyStressM   ||
                out == MaterialOutput::StrainM             )
        strain = (m_data.mine().m_Bcov_ori - m_data.mine().m_Bcov_def);
    else if (   out == MaterialOutput::Generic          ||
                out == MaterialOutput::Strain           ||
                out == MaterialOutput::Stress           ||
                out == MaterialOutput::PStress             ) // To be used with multiplyLinZ_into or integrateZ_into
        strain = 0.5*(m_data.mine().m_Gcov_def.block(0,0,2,2) - m_data.mine().m_Gcov_ori.block(0,0,2,2));
    else
        GISMO_ERROR("Output type MaterialOutput::" + std::to_string((short_t)(out)) + " not understood. See gsMaterialMatrixUtils.h");

    return strain;
}

// template <short_t dim, class T>
// gsMatrix<T> gsMaterialMatrixLinear<dim,T>::_E(const gsMatrix<T> & C, const T z, enum MaterialOutput out) const
// {
//     gsMatrix<T> strain;
//     if      (out == MaterialOutput::VectorN || out == MaterialOutput::PStrainN || out == MaterialOutput::PStressN || out == MaterialOutput::TensionField || out == MaterialOutput::StrainN) // To be used with multiplyZ_into
//         strain = 0.5*(C - m_data.mine().m_Acov_ori);
//     else if (out == MaterialOutput::VectorM || out == MaterialOutput::PStrainM || out == MaterialOutput::PStressM || out == MaterialOutput::TensionField || out == MaterialOutput::StrainM) // To be used with multiplyZ_into
//         strain = (C - m_data.mine().m_Bcov_def);
//     else if (out == MaterialOutput::Generic || out == MaterialOutput::Strain) // To be used with multiplyLinZ_into or integrateZ_into
//         strain = 0.5*(C - m_data.mine().m_Gcov_ori);
//     else
//         GISMO_ERROR("Output type is not VectorN, PStrainN, PStressN, VectorM, PStrainM. PStressM, TensionField or Generic!");

//     return strain;
// }

//--------------------------------------------------------------------------------------------------------------------------------------

namespace internal
{

/// @brief get a Linear Material Matrix from XML data
///
/// \ingroup KLShell
template<short_t d, class T>
class gsXml< gsMaterialMatrixLinear<d,T> >
{
private:
    gsXml() { }
    typedef gsMaterialMatrixLinear<d,T> Object;

public:
    GSXML_COMMON_FUNCTIONS(gsMaterialMatrixLinear<TMPLA2(d,T)>);
    static std::string tag ()  { return "MaterialMatrix"; }
    static std::string type () { return "Linear" +  to_string(d); }

    GSXML_GET_POINTER(Object);

    // static Object * get(gsXmlNode * node)
    // {
    //     Object result;
    //     get_into(node, result);
    //     return result.clone().release();
    // }

    static void get_into(gsXmlNode * node,Object & obj)
    {
        obj = getMaterialMatrixFromXml< Object >( node );
    }

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data)
    {
        return putMaterialMatrixToXml< Object >( obj,data );
        // GISMO_NO_IMPLEMENTATION;
        // return putGeometryToXml(obj,data);
    }
};

}// namespace internal

} // end namespace
