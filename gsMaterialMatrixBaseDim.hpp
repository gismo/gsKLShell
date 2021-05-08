/** @file gsMaterialMatrix.hpp

    @brief Provides hyperelastic material matrices

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
    - Make beta (compressible materials) and material parameters universal for all integration points over the thickness. So get them out of the _dPsi functions etc and move them into the integration loops as global variables.

*/



#pragma once

#include <gsKLShell/gsMaterialMatrixBaseDim.h>

namespace gismo
{

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    m_map_def.flags = m_map.flags;
    m_map_def.points = u;
    static_cast<const gsFunction<T>&>(Base::m_defpatches->piece(patch)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    _computeMetricDeformed_impl<dim>(patch,u,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.points.cols());    m_Acon_def_mat.setZero();
    m_Bcov_def_mat.resize(4,m_map_def.points.cols());    m_Bcov_def_mat.setZero();

    m_acov_def_mat.resize(2*3,m_map_def.points.cols());    m_acov_def_mat.setZero();
    if (basis)
    {
        m_acon_def_mat.resize(2*3,m_map_def.points.cols());    m_acon_def_mat.setZero();
        m_ncov_def_mat.resize(2*3,m_map_def.points.cols());    m_ncov_def_mat.setZero();
    }

    for (index_t k=0; k!= m_map_def.points.cols(); k++)
    {
        m_acov_def_mat.reshapeCol(k,3,2)   = m_map_def.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,3,2);
        acov = m_map_def.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map_def.deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map_def.normal(k).normalized();

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_Bcov_def_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_def_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
            gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,3,2);
            for (index_t i=0; i < 2; i++)
                acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

            // Mixed tensor
            for (index_t i=0; i < 2; i++)
                for (index_t j=0; j < 2; j++)
                    mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

            gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_def_mat.reshapeCol(k,3,2);
            for (index_t i=0; i < 2; i++)
                ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);
        }
    }
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.points.cols());    m_Acon_def_mat.setZero();

    m_acov_def_mat.resize(2*2,m_map_def.points.cols());    m_acov_def_mat.setZero();
    if (basis)
    {
        m_acon_def_mat.resize(2*2,m_map_def.points.cols());    m_acon_def_mat.setZero();
    }

    for (index_t k=0; k!= m_map_def.points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,2,2);
        acov = m_map_def.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
            gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,2,2);
            for (index_t i=0; i < 2; i++)
                acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    m_map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    _computeMetricUndeformed_impl<dim>(patch,u,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed_impl(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map.points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map.points.cols());    m_Acon_ori_mat.setZero();
    m_Bcov_ori_mat.resize(4,m_map.points.cols());    m_Bcov_ori_mat.setZero();

    m_acov_ori_mat.resize(2*3,m_map.points.cols());    m_acov_ori_mat.setZero();
    if (basis)
    {
        m_acon_ori_mat.resize(2*3,m_map.points.cols());    m_acon_ori_mat.setZero();
        m_ncov_ori_mat.resize(2*3,m_map.points.cols());    m_ncov_ori_mat.setZero();
    }

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        m_acov_ori_mat.reshapeCol(k,3,2)   = m_map.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,3,2);
        acov = m_map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map.deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map.normal(k).normalized();

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_Bcov_ori_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_ori_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
           gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_ori_mat.reshapeCol(k,3,2);
           for (index_t i=0; i < 2; i++)
               acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

           // Mixed tensor
           for (index_t i=0; i < 2; i++)
               for (index_t j=0; j < 2; j++)
                   mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

           gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_ori_mat.reshapeCol(k,3,2);
           for (index_t i=0; i < 2; i++)
               ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);
       }
    }
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed_impl(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map.points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map.points.cols());    m_Acon_ori_mat.setZero();

    m_acov_ori_mat.resize(2*2,m_map.points.cols());    m_acov_ori_mat.setZero();
    if (basis)
    {
        m_acon_ori_mat.resize(2*2,m_map.points.cols());    m_acon_ori_mat.setZero();
    }

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,2,2);
        acov = m_map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
            gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_ori_mat.reshapeCol(k,2,2);
            for (index_t i=0; i < 2; i++)
                acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetric(index_t k, T z, bool basis) const
{
    this->_getMetricDeformed(k,z,basis);
    this->_getMetricUndeformed(k,z,basis);

    T ratio = m_Gcov_def.determinant() / m_Gcov_ori.determinant();
    GISMO_ENSURE(ratio > 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_Gcov_ori.determinant());
    m_J0_sq = ratio;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed(index_t k, T z, bool basis) const
{
    _getMetricDeformed_impl<dim>(k,z,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed_impl(index_t k, T z, bool basis) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    if (basis)
    {
        GISMO_ENSURE(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_acon_def_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_ncov_def_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    m_Acon_def = m_Acon_def_mat.reshapeCol(k,2,2);
    m_Bcov_def = m_Bcov_def_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def - 2.0 * z * m_Bcov_def + z*z * m_ncov_def.transpose()*m_ncov_def;
    m_Gcov_def(2,2) = 1.0;
    m_Gcon_def = m_Gcov_def.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_acov_def = m_acov_def_mat.reshapeCol(k,3,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,3,2);
    m_ncov_def = m_ncov_def_mat.reshapeCol(k,3,2);
    // g
    gsMatrix<T,3,1> normal = m_map_def.normal(k).normalized();
    m_gcov_def.leftCols(2) = m_acov_def + z * m_ncov_def;
    m_gcov_def.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_def.col(c) = m_Gcon_def(c,0) * m_gcov_def.col(0)
                            + m_Gcon_def(c,1) * m_gcov_def.col(1)
                            + m_Gcon_def(c,2) * m_gcov_def.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_def_L = m_Gcov_def;
    // m_Gcov_def.block(0,0,2,2) -= z*z * m_ncov_def.transpose()*m_ncov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed_impl(index_t k, T z, bool basis) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");

    if (basis)
    {
        GISMO_ENSURE(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    m_Acon_def = m_Acon_def_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def;
    m_Gcov_def(2,2) = 1.0;
    m_Gcon_def = m_Gcov_def.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_acov_def = m_acov_def_mat.reshapeCol(k,2,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,2,2);
    // g
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_gcov_def.setZero();
    m_gcov_def.block(0,0,2,2) = m_acov_def;
    m_gcov_def.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_def.col(c) = m_Gcon_def(c,0) * m_gcov_def.col(0)
                            + m_Gcon_def(c,1) * m_gcov_def.col(1)
                            + m_Gcon_def(c,2) * m_gcov_def.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_def_L = m_Gcov_def;
    // m_Gcov_def.block(0,0,2,2) -= z*z * m_ncov_def.transpose()*m_ncov_def;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed(index_t k, T z, bool basis) const
{
    _getMetricUndeformed_impl<dim>(k,z,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed_impl(index_t k, T z, bool basis) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");

    if (basis)
    {
        GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    m_Acon_ori = m_Acon_ori_mat.reshapeCol(k,2,2);
    m_Bcov_ori = m_Bcov_ori_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori - 2.0 * z * m_Bcov_ori + z*z * m_ncov_ori.transpose()*m_ncov_ori;
    m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori = m_Gcov_ori.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_acov_ori = m_acov_ori_mat.reshapeCol(k,3,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,3,2);
    m_ncov_ori = m_ncov_ori_mat.reshapeCol(k,3,2);
    // g
    gsMatrix<T,3,1> normal = m_map.normal(k).normalized();
    m_gcov_ori.block(0,0,3,2) = m_acov_ori + z * m_ncov_ori;
    m_gcov_ori.col(2) = normal;
    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_ori.col(c) = m_Gcon_ori(c,0) * m_gcov_ori.col(0)
                            + m_Gcon_ori(c,1) * m_gcov_ori.col(1)
                            + m_Gcon_ori(c,2) * m_gcov_ori.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_ori_L = m_Gcov_ori;
    // m_Gcov_ori.block(0,0,2,2) -= z*z * m_ncov_ori.transpose()*m_ncov_ori;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed_impl(index_t k, T z, bool basis) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");

    if (basis)
    {
        GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    m_Acon_ori = m_Acon_ori_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori;
    m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori = m_Gcov_ori.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_acov_ori = m_acov_ori_mat.reshapeCol(k,2,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,2,2);
    // g
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_gcov_ori.setZero();
    m_gcov_ori.block(0,0,2,2) = m_acov_ori;
    m_gcov_ori.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_ori.col(c) = m_Gcon_ori(c,0) * m_gcov_ori.col(0)
                            + m_Gcon_ori(c,1) * m_gcov_ori.col(1)
                            + m_Gcon_ori(c,2) * m_gcov_ori.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_ori_L = m_Gcov_ori;
    // m_Gcov_ori.block(0,0,2,2) -= z*z * m_ncov_ori.transpose()*m_ncov_ori;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixBaseDim<dim,T>::_evalStretch(const gsMatrix<T> & C) const
{
    gsVector<T> stretches;
    gsMatrix<T> stretchvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    stretches.resize(3,1);    stretches.setZero();
    stretchvec.resize(3,3);   stretchvec.setZero();

    Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;


    GISMO_ENSURE(m_gcon_ori.cols()!=0,"Is the basis initialized?");

    gsMatrix<T,3,3> B;
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += C(k,l) * m_gcon_ori.col(k) * m_gcon_ori.col(l).transpose();

    eigSolver.compute(B);

    stretchvec.leftCols(2) = eigSolver.eigenvectors().rightCols(2);
    stretchvec.col(2) = m_gcon_ori.col(2); // replace with: stretchvec.col(0).template head<3>().cross(stretchvec.col(1).template head<3>())
    stretches.block(0,0,2,1) = eigSolver.eigenvalues().block(1,0,2,1); // the eigenvalues are a 3x1 matrix, so we need to use matrix block-operations

    // m_stretches.at(2) = 1/m_J0_sq;
    stretches.at(2) = C(2,2);

    for (index_t k=0; k!=3; k++)
        stretches.at(k) = math::sqrt(stretches.at(k));

    result.first = stretches;
    result.second = stretchvec;

    // // DEBUGGING ONLY!
    // gsMatrix<T> ones(3,1);
    // ones.setOnes();
    // gsDebugVar(m_stretchvec);
    // gsDebugVar(result.first);

    return result;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeStretch(const gsMatrix<T> & C) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalStretch(C);
    m_stretches = result.first;
    m_stretchvec = result.second;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixBaseDim<dim,T>::_transformation(const gsMatrix<T> & basis1, const gsMatrix<T> & basis2) const
{
    // Transformation of a quantity FROM basis2 TO basis1
    gsMatrix<T> Tmat(3,3);

    for (index_t i = 0; i!=2; i++)
        for (index_t j = 0; j!=2; j++)
            Tmat(i,j) = math::pow(basis2.col(i).dot(basis1.col(j)),2);

    Tmat(2,0)   = basis2.col(1).dot(basis1.col(0)) * basis2.col(0).dot(basis1.col(0));
    Tmat(2,1)   = basis2.col(0).dot(basis1.col(1)) * basis2.col(1).dot(basis1.col(1));

    Tmat(0,2)   = 2*basis2.col(0).dot(basis1.col(1)) * basis2.col(0).dot(basis1.col(0));
    Tmat(1,2)   = 2*basis2.col(1).dot(basis1.col(0)) * basis2.col(1).dot(basis1.col(1));

    Tmat(2,2)   = basis2.col(0).dot(basis1.col(1)) * basis2.col(1).dot(basis1.col(0))
                 +basis2.col(0).dot(basis1.col(0)) * basis2.col(1).dot(basis1.col(1));

    return Tmat;
}

} // end namespace
