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
                                        const gsFunction<T> & thickness,
                                        const std::vector<gsFunction<T>*> &pars
                                        )
                                        :
                                        Base(mp),
                                        m_thickness(&thickness),
                                        m_pars(pars)
{
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
                                    Base(mp),
                                    m_thickness(&thickness),
                                    m_pars(pars),
                                    m_density(&density)
{
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

    // Initialize some parameters
    m_numPars = m_pars.size();
}


template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::_computePoints(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMatrix<T> tmp;

    this->_computeMetricUndeformed(patch,u,basis);

    if (Base::m_defpatches->nPieces()!=0)
        this->_computeMetricDeformed(patch,u,basis);

    m_thickness->eval_into(m_map.mine().values[0], m_Tmat);

    m_parmat.resize(m_numPars,m_map.mine().values[0].cols());
    m_parmat.setZero();

    for (size_t v=0; v!=m_pars.size(); v++)
    {
        m_pars[v]->eval_into(m_map.mine().values[0], tmp);
        m_parmat.row(v) = tmp;
    }

    m_parvals.resize(m_numPars);
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.mine().flags = NEED_VALUE;
    m_map.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map);

    result.resize(1, u.cols());
    m_thickness->eval_into(m_map.mine().values[0], m_Tmat);
    m_density->eval_into(m_map.mine().values[0], m_rhomat);
    for (index_t i = 0; i != u.cols(); ++i) // points
    {
        result(0,i) = m_Tmat(0,i)*m_rhomat(0,i);
    }

}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map);

    this->_computePoints(patch,u,true);

    result.resize(3, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0,true); // on point i, with height 0.0

        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_J0_sq;

        res = this->_evalStretch(C);
        result.col(i) = res.first;
    }
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map);

    this->_computePoints(patch,u,true);

    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0,true); // on point i, with height 0.0

        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_J0_sq;

        res = this->_evalStretch(C);
        result.col(i) = res.second.reshape(9,1);
    }
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map.mine().flags = NEED_VALUE;
    m_map.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map);
    m_thickness->eval_into(m_map.mine().values[0], result);
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) covariant basis
template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::covtransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.resize(9, u.cols());
    gsMatrix<T> tmp, conbasis,sbasis;
    this->stretchDir_into(patch,u,tmp);
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0,true); // on point i, with height 0.0
        sbasis = tmp.reshapeCol(i,3,3);
        conbasis = m_gcov_ori;
        result.col(i) = this->_transformation(conbasis,sbasis).reshape(9,1);
    }
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::contransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.resize(9, u.cols());
    gsMatrix<T> tmp, conbasis,sbasis;
    this->stretchDir_into(patch,u,tmp);
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0,true); // on point i, with height 0.0
        sbasis = tmp.reshapeCol(i,3,3);
        conbasis = m_gcon_ori;
        result.col(i) = this->_transformation(conbasis,sbasis).reshape(9,1);
    }
}

template <short_t dim, class T >
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
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
                // this->_getMetric(k, z(j, k), false); // on point i, on height z(0,j)
                this->_getMetric(k, z(j, k) * m_Tmat(0, k), false); // on point i, on height z(0,j)

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

    this->_computePoints(patch,u,false);
    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // this->_getMetric(k, z(j, k), false); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_Tmat(0, k), false); // on point i, on height z(0,j)

            result(0, j * u.cols() + k) = _Sij(0, 0, z(j, k), out);
            result(1, j * u.cols() + k) = _Sij(1, 1, z(j, k), out);
            result(2, j * u.cols() + k) = _Sij(0, 1, z(j, k), out);
        }
    }

    return result;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{

    this->_computePoints(patch,u,true);
    gsMatrix<T> result(3, u.cols() * z.rows());
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
            // this->_getMetric(k, z(j, k), true); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_Tmat(0, k), true); // on point i, on height z(0,j)

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

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{

    this->_computePoints(patch,u,true);
    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> E;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // this->_getMetric(k, z(j, k), true); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_Tmat(0, k), true); // on point i, on height z(0,j)

            // E.setZero();
            // E.block(0,0,2,2) = _E(0,out);
            // E(2, 2) = 0;
            // res = _evalPStrain(E);
            // result.col(j * u.cols() + k) = res.first;
            //
            //
            gsMatrix<> E = _E(0,out);
            gsDebugVar("Temporary");
            result(0,j*u.cols() + k) = E(0,0);
            result(1,j*u.cols() + k) = E(1,1);
            result(2,j*u.cols() + k) = E(0,1) + E(1,0);
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
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // this->_getMetric(k, z(j, k), true); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_Tmat(0, k), true); // on point i, on height z(0,j)

            S.setZero();
            S(0, 0) = _Sij(0, 0, 0, out);
            S(0, 1) = _Sij(0, 1, 0, out);
            S(1, 0) = _Sij(1, 0, 0, out);
            S(1, 1) = _Sij(1, 1, 0, out);
            S(2, 2) = 0;
            res = _evalPStress(S);
            Sp = res.first;

            E.setZero();
            E.block(0,0,2,2) = _E(0,out);
            E(2, 2) = 0;
            res = _evalPStrain(E);
            Ep = res.first;

            // See Nakashino 2020
            // Smin = Sp[0], Smax = Sp[1], S33 = Sp[2]
            // Emin = Ep[0], Emax = Ep[1], E33 = Ep[2]
            if (Sp[0] > 0) // taut
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
    gsMatrix<T> strain = _E(z,out);
    T result =  _Cijkl(i,j,0,0) * strain(0,0) + _Cijkl(i,j,0,1) * strain(0,1)
                + _Cijkl(i,j,1,0) * strain(1,0) + _Cijkl(i,j,1,1) * strain(1,1);

    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::_E(const T z, enum MaterialOutput out) const
{
    gsMatrix<T> strain;
    if      (out == MaterialOutput::VectorN || out == MaterialOutput::PStrainN || out == MaterialOutput::PStressN || out == MaterialOutput::TensionField) // To be used with multiplyZ_into
        strain = 0.5*(m_Acov_def - m_Acov_ori);
    else if (out == MaterialOutput::VectorM || out == MaterialOutput::PStrainM || out == MaterialOutput::PStressM || out == MaterialOutput::TensionField) // To be used with multiplyZ_into
        strain = (m_Bcov_ori - m_Bcov_def);
    else if (out == MaterialOutput::Generic) // To be used with multiplyLinZ_into or integrateZ_into
        strain = 0.5*(m_Acov_def - m_Acov_ori) + z*(m_Bcov_ori - m_Bcov_def);
    else
        GISMO_ERROR("Output type is not VectorN, PStrainN, PStressN, VectorM, PStrainM. PStressM, TensionField or Generic!");

    return strain;
}

//--------------------------------------------------------------------------------------------------------------------------------------


template <short_t dim, class T >
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixLinear<dim,T>::_evalPStress(const gsMatrix<T> & S) const
{
    gsVector<T> pstresses;
    gsMatrix<T> pstressvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    pstresses.resize(3,1);    pstresses.setZero();
    pstressvec.resize(3,3);   pstressvec.setZero();

    Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

    gsMatrix<T> B(3,3);
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += S(k,l) * m_gcov_ori.col(k) * m_gcov_ori.col(l).transpose();

    eigSolver.compute(B);

    index_t zeroIdx = -1;
    real_t tol = 1e-14;
    real_t max = eigSolver.eigenvalues().array().abs().maxCoeff();
    max = (max==0) ? 1 : max;
    for (index_t k=0; k!=3; k++)
        zeroIdx = std::abs(eigSolver.eigenvalues()[k] ) / max < tol ? k : zeroIdx;

    GISMO_ASSERT(zeroIdx!=-1,"No zero found?");

    index_t count = 0;
    pstressvec.col(2) = m_gcon_ori.col(2);
    pstresses(2,0) = S(2,2);

    for (index_t k=0; k!=3; k++)
    {
        if (k==zeroIdx) continue;
        pstressvec.col(count) = eigSolver.eigenvectors().col(k);
        pstresses(count,0) = eigSolver.eigenvalues()(k,0);
        count++;
    }

    result.first = pstresses;
    result.second = pstressvec;

    return result;
}

template <short_t dim, class T >
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixLinear<dim,T>::_evalPStrain(const gsMatrix<T> & S) const
{
    gsVector<T> pstrains;
    gsMatrix<T> pstrainvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    pstrains.resize(3,1);    pstrains.setZero();
    pstrainvec.resize(3,3);   pstrainvec.setZero();

    Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

    gsMatrix<T> B(3,3);
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += S(k,l) * m_gcon_ori.col(k) * m_gcon_ori.col(l).transpose();

    eigSolver.compute(B);

    index_t zeroIdx = -1;
    real_t tol = 1e-14;
    real_t max = eigSolver.eigenvalues().array().abs().maxCoeff();
    max = (max==0) ? 1 : max;
    for (index_t k=0; k!=3; k++)
        zeroIdx = std::abs(eigSolver.eigenvalues()[k] ) / max < 1e-14 ? k : zeroIdx;

    GISMO_ASSERT(zeroIdx!=-1,"No zero found?");

    index_t count = 0;
    pstrainvec.col(2) = m_gcon_ori.col(2);
    pstrains(2,0) = S(2,2);

    for (index_t k=0; k!=3; k++)
    {
        if (k==zeroIdx) continue;
        pstrainvec.col(count) = eigSolver.eigenvectors().col(k);
        pstrains(count,0) = eigSolver.eigenvalues()(k,0);
        count++;
    }

    result.first = pstrains;
    result.second = pstrainvec;

    return result;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::_computePStress(const gsMatrix<T> & S) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalPStress(S);
    m_pstress = result.first;
    m_pstressvec = result.second;
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::_computePStrain(const gsMatrix<T> & E) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalPStrain(E);
    m_pstrain = result.first;
    m_pstrainvec = result.second;
}

} // end namespace
