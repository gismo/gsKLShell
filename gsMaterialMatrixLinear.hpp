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

using namespace gismo;

template <short_t d, typename T>
class Cfun : public gsFunction<T>
{
public:
    Cfun(const gsMaterialMatrixBaseDim<d,T> * materialMat, index_t patch, const gsVector<T> & u, const T & z)
    :
    m_materialMat(materialMat),
    m_patchID(patch),
    m_points(u),
    m_z(z)
    {

    }

    Cfun(const gsMaterialMatrixBaseDim<d,T> * materialMat, index_t patch, const gsVector<T> & u)
    :
    m_materialMat(materialMat),
    m_patchID(patch),
    m_points(u)
    {
        m_z.resize(1,1);
        m_z.setZero();
    }

    void eval_into(const gsMatrix<T>& C, gsMatrix<T>& result) const
    {
        result.resize(targetDim(),C.cols());
        for (index_t k=0; k!=C.cols(); k++)
            result.col(k) = m_materialMat->eval3D_matrix_C(C.col(k),m_patchID,m_points,m_z,MaterialOutput::MatrixA);
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 9;
    }

private:
    const gsMaterialMatrixBaseDim<d,T> * m_materialMat;
    mutable index_t m_patchID;
    mutable gsVector<T> m_points;
    mutable T m_z;
};


template <short_t d, typename T>
class Sfun : public gsFunction<T>
{
public:
    Sfun(const gsMaterialMatrixLinear<d,T> * materialMat, index_t patch, const gsVector<T> & u, const T & z)
    :
    m_materialMat(materialMat),
    m_patchID(patch),
    m_points(u),
    m_z(z)
    {

    }

    Sfun(const gsMaterialMatrixLinear<d,T> * materialMat, index_t patch, const gsVector<T> & u)
    :
    m_materialMat(materialMat),
    m_patchID(patch),
    m_points(u),
    m_z(0)
    {
    }

    void eval_into(const gsMatrix<T>& C, gsMatrix<T>& result) const
    {
        result.resize(targetDim(),C.cols());
        for (index_t k=0; k!=C.cols(); k++)
            result.col(k) = m_materialMat->eval3D_vector_C(C.col(k),m_patchID,m_points,m_z,MaterialOutput::VectorN);
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 3;
    }

private:
    const gsMaterialMatrixLinear<d,T> * m_materialMat;
    mutable index_t m_patchID;
    mutable gsVector<T> m_points;
    mutable T m_z;
};


namespace gismo
{

template <short_t dim, class T>
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

template <short_t dim, class T>
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

template <short_t dim, class T>
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear( const gsMaterialMatrixLinear<dim,T> & other)
{
    *this =(other);
}

template <short_t dim, class T>
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

template <short_t dim, class T>
void gsMaterialMatrixLinear<dim,T>::_defaultOptions()
{
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T>
void gsMaterialMatrixLinear<dim,T>::_initialize()
{
    // Set default options
    this->_defaultOptions();
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->_computePoints(patch,u,false);

    gsMatrix<T> result(9, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
            m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
                // _getMetric(k, z(j, k), false); // on point i, on height z(0,j)
                this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k), false); // on point i, on height z(0,j)

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

    this->_computePoints(patch,u,false);

    gsMatrix<T> result(9, 1);
    result.setZero();

    // Evaluate material properties on the quadrature point
    for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
        m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,0);

    // Get metric
    this->_getMetricUndeformed(0, z * m_data.mine().m_Tmat(0, 0), false); // on point i, on height z(0,j)
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
    // return gsMatrix<T>::Zero(27,u.cols()*z.rows());

    this->_computePoints(patch,u,false);

    gsMatrix<T> result(27,u.cols()*z.rows());
    result.setZero();

    gsMatrix<T> Cvoight(3,1), Ctmp;
    for (index_t k=0; k!=u.cols(); k++)
    {

        // // Evaluate material properties on the quadrature point
        // for (index_t v=0; v!=m_parmat.rows(); v++)
        //     m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k), false); // on point i, on height z(0,j)

            Cfun<dim,T> Cfunc(this,patch,u.col(k),z(j,k));
            gsMatrix<T> tmpresult;

            Cvoight(0,0) = m_data.mine().m_Gcov_def(0,0);
            Cvoight(1,0) = m_data.mine().m_Gcov_def(1,1);
            Cvoight(2,0) = m_data.mine().m_Gcov_def(0,1);


            Cfunc.deriv_into(Cvoight,tmpresult);
            result.col(j*u.cols()+k) = tmpresult;
        }
    }

    return result;
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

//     this->_computePoints(patch,u,false);

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
//             this->_getMetric(k, z(j, k) * m_Tmat(0, k), false); // on point i, on height z(0,j)

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

    this->_computePoints(patch,u,false);

    gsMatrix<T> result(3, 1);
    result.setZero();

    gsMatrix<T> E;
    // Evaluate material properties on the quadrature point
    for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
        m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,0);

    // Get metric
    this->_getMetricUndeformed(0, z * m_data.mine().m_Tmat(0, 0), false); // on point i, on height z(0,j)
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
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    gsMatrix<T> Emat = eval3D_strain(patch,u,z,out);

    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> E;
    index_t colIdx;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j * u.cols() + k;
            E.setZero();
            E(0,0) = Emat(0,colIdx);
            E(1,1) = Emat(1,colIdx);
            E(0,1) = E(1,0) = 0.5*Emat(2,colIdx);
            E(2,2) = 0;
            res = this->_evalPStrain(E);
            result.col(j * u.cols() + k) = res.first;
        }
    }

    return result;

}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    this->_computePoints(patch,u,true);

    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,2,2> E;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // this->_getMetric(k, z(j, k), true); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k), true); // on point i, on height z(0,j)

            // E.setZero();
            E.block(0,0,2,2) = _E(0,out);
            result(0,j*u.cols() + k) = E(0,0);
            result(1,j*u.cols() + k) = E(1,1);
            result(2,j*u.cols() + k) = E(0,1) + E(1,0);
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

    this->_computePoints(patch,u,true);

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
            // this->_getMetric(k, z(j, k), true); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k), true); // on point i, on height z(0,j)

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
    this->_computePoints(patch,u,true);
    gsMatrix<T> result(1, u.cols() * z.rows());
    result.setZero();
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k), true); // on point i, on height z(0,j)
            result(0,j * u.cols() + k) = math::sqrt(m_data.mine().m_J0_sq*1.0);
        }
    }
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    this->_computePoints(patch,u,true);
    gsMatrix<T> Smat = eval3D_stress(patch,u,z,out);

    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();
    index_t colIdx;
    T detF;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j * u.cols() + k;
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k), true); // on point i, on height z(0,j)
            detF = math::sqrt(m_data.mine().m_J0_sq*1.0);
            Smat.col(colIdx) /= detF;
        }
    }
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    this->_computePoints(patch,u,true);

    gsMatrix<T> result(1, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> S;
    gsMatrix<T,3,3> E;
    gsVector<T> Sp, Ep;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_data.mine().m_parmat.rows(); v++)
            m_data.mine().m_parvals.at(v) = m_data.mine().m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // this->_getMetric(k, z(j, k), true); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_data.mine().m_Tmat(0, k), true); // on point i, on height z(0,j)

            E.setZero();
            E.block(0,0,2,2) = _E(0,out);
            E(2, 2) = 0;
            res = this->_evalPStrain(E);
            Ep = res.first;

            S.setZero();
            S(0, 0) = _Sij(0, 0, E.block(0,0,2,2));
            S(0, 1) = _Sij(0, 1, E.block(0,0,2,2));
            S(1, 0) = _Sij(1, 0, E.block(0,0,2,2));
            S(1, 1) = _Sij(1, 1, E.block(0,0,2,2));
            S(2, 2) = 0;
            res = this->_evalPStress(S);
            Sp = res.first;

            // See Nakashino 2020
            // Smin = Sp[0], Smax = Sp[1], S33 = Sp[2]
            // Emin = Ep[0], Emax = Ep[1], E33 = Ep[2]
            if (Sp[0] >= 0) // taut
                result.col(j * u.cols() + k) << 1;
            else if (Ep[1] < 0) // slack
                result.col(j * u.cols() + k) << -1;
            else // wrinkled
                result.col(j * u.cols() + k) << 0;
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
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::S(const gsMatrix<T> & C, const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    this->_computePoints(patch,u,false);

    gsMatrix<T> tmpresult;
    Sfun<dim,T> Sfunc(this,patch,u.col(0),z(0,0));
    Sfunc.eval_into(C,tmpresult);
    gsDebugVar(tmpresult.array());

    gsMatrix<T> result;
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::C(const gsMatrix<T> &  C, const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    GISMO_ASSERT(C.cols()==1,"Expect C in voight notation, C = [C11, C22, C12]");
    gsWarn<<"This function needs to compute the material parameters! Therefore needs compute points\n";

    this->_computePoints(patch,u,false);

    gsMatrix<T> tmpresult;
    Sfun<dim,T> Sfunc(this,patch,u.col(0),z(0,0));
    Sfunc.deriv_into(C,tmpresult);
    // Because the derivative of S is taken w.r.t. the vector C = C11, C22, C12, since E = 0.5 ( C - G), and since voight notation is used, we need to multiply the last column of the derivative with 0.5 to revover the correct matrix. This is because the last column is dS_ij / dE_12 = 2C_1112; 2C_2212; 2C_1212.
    gsAsMatrix<T> Cmat = tmpresult.reshapeCol(0,3,3);
    Cmat.col(2) *= 0.5;
    gsDebugVar(tmpresult.array()*2);

    Cfun<dim,T> Cfunc(this,patch,u.col(0),z(0,0));
    Cfunc.eval_into(C,tmpresult);
    gsDebugVar(tmpresult.array());

    Cfunc.deriv_into(C,tmpresult);
    gsDebugVar(tmpresult.array());

    gsMatrix<T> result(3,3);
    result(0,0)                 = _Cijkl(0,0,0,0); // C1111
    result(1,1)                 = _Cijkl(1,1,1,1); // C2222
    result(2,2)                 = _Cijkl(0,1,0,1); // C1212
    result(1,0) = result(0,1)   = _Cijkl(0,0,1,1); // C1122
    result(2,0) = result(0,2)   = _Cijkl(0,0,0,1); // C1112
    result(2,1) = result(1,2)   = _Cijkl(1,1,0,1); // C2212
    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::dC(const gsMatrix<T> & C, const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    GISMO_ASSERT(C.cols()==1,"Expect C in voight notation, C = [C11, C22, C12]");

    // Put C in matrix notation
    gsMatrix<T> Cmat(3,3);
    Cmat.setZero();
    Cmat(0,0) = C(0,0);
    Cmat(1,1) = C(1,0);
    Cmat(0,1) = Cmat(1,0) = C(2,0);
    Cmat(2,2) = 1;

    // gsMatrix<T> tmpresult;
    // Cfun<dim,T> fun(this,patch,u,z);
    // fun.eval_into(Cmat,tmpresult);


    // Cfun<dim,T> Cfunc(this,patch,u,z);
    // Cfunc.eval_into(Cmat,tmpresult);
    // gsDebugVar(tmpresult);

    gsMatrix<T> result;
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
    else if (out == MaterialOutput::Generic || out == MaterialOutput::Strain || out == MaterialOutput::Stress) // To be used with multiplyLinZ_into or integrateZ_into
        strain = 0.5*(m_data.mine().m_Acov_def - m_data.mine().m_Acov_ori) + z*(m_data.mine().m_Bcov_ori - m_data.mine().m_Bcov_def);
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

} // end namespace
