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
#include <gsCore/gsGeometry.h>

namespace gismo
{

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    gsMapData<T> map;
    map.flags = NEED_VALUE;
    map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    result.resize(1, u.cols());
    m_thickness->piece(patch).eval_into(map.values[0], m_data.mine().m_Tmat);
    m_density->piece(patch).eval_into(map.values[0], m_data.mine().m_rhomat);
    for (index_t i = 0; i != u.cols(); ++i) // points
    {
        result(0,i) = m_data.mine().m_Tmat(0,i)*m_data.mine().m_rhomat(0,i);
    }
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    gsMapData<T> map;
    map.flags = NEED_VALUE;
    map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);
    m_thickness->piece(patch).eval_into(map.values[0], result);
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(3, u.cols());
    gsMatrix<T> tmp, C;
    std::pair<gsMatrix<T>,gsMatrix<T>> res;
    deformation_into(patch,u,tmp);
    for (index_t i=0; i!= u.cols(); i++)
    {
        _getMetric(i,0.0,true); // on point i, with height 0.0
        C = tmp.reshapeCol(i,3,3);
        res = this->_evalStretch(C);
        result.col(i) = res.second.reshape(3,1);
    }
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(9, u.cols());
    gsMatrix<T> tmp, C;
    std::pair<gsMatrix<T>,gsMatrix<T>> res;
    deformation_into(patch,u,tmp);
    for (index_t i=0; i!= u.cols(); i++)
    {
        _getMetric(i,0.0,true); // on point i, with height 0.0
        C = tmp.reshapeCol(i,3,3);
        res = this->_evalStretch(C);
        result.col(i) = res.second.reshape(9,1);
    }
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.resize(9, u.cols());
    gsMatrix<T> tmp, conbasis,sbasis;
    stretchDir_into(patch,u,tmp);
    for (index_t i=0; i!= u.cols(); i++)
    {
        _getMetric(i,0.0,true); // on point i, with height 0.0
        sbasis = tmp.reshapeCol(i,3,3);
        conbasis = m_data.mine().m_gcon_ori;
        result.col(i) = _transformation(conbasis,sbasis).reshape(9,1);
    }
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::parameters_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    gsMapData<T> map;
    map.flags = NEED_VALUE;
    map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    gsMatrix<T> tmp;
    result.resize(m_pars.size(),map.values[0].cols());
    result.setZero();
    for (size_t v=0; v!=m_pars.size(); v++)
    {
        m_pars[v]->piece(patch).eval_into(map.values[0], tmp);
        result.row(v) = tmp;
    }
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::deformation_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    _computePoints(patch,u,true);

    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        _getMetric(i,0.0,true); // on point i, with height 0.0

        gsAsMatrix<T> C = result.reshapeCol(i,3,3);
        C.setZero();
        C.block(0,0,2,2) = m_data.mine().m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_data.mine().m_J0_sq;
    }
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::_computePoints(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    this->_computeMetricUndeformed(patch,u,basis);
    if (Base::m_defpatches->nPieces()!=0)
        this->_computeMetricDeformed(patch,u,basis);

    gsMapData<T> map;
    map.flags = NEED_VALUE;
    map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    m_thickness->piece(patch).eval_into(map.values[0], m_data.mine().m_Tmat);

    m_data.mine().m_parmat.resize(m_pars.size(),map.values[0].cols());
    m_data.mine().m_parmat.setZero();

    gsMatrix<T> tmp;
    for (size_t v=0; v!=m_pars.size(); v++)
    {
        m_pars[v]->piece(patch).eval_into(map.values[0], tmp);
        m_data.mine().m_parmat.row(v) = tmp;
    }

    m_data.mine().m_parvals.resize(m_pars.size());
}

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    _computeMetricDeformed_impl<dim>(patch,u,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_defpatches->piece(patch)   ).computeMap(map);

    gsMatrix<T> deriv2;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_data.mine().m_normal_def_mat = map.normals;
    m_data.mine().m_normal_def_mat.colwise().normalized();

    m_data.mine().m_Acov_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_def_mat.setZero();
    m_data.mine().m_Acon_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_def_mat.setZero();
    m_data.mine().m_Bcov_def_mat.resize(4,map.points.cols());    m_data.mine().m_Bcov_def_mat.setZero();

    m_data.mine().m_acov_def_mat.resize(2*3,map.points.cols());    m_data.mine().m_acov_def_mat.setZero();
    if (basis)
    {
        m_data.mine().m_acon_def_mat.resize(2*3,map.points.cols());    m_data.mine().m_acon_def_mat.setZero();
        m_data.mine().m_ncov_def_mat.resize(2*3,map.points.cols());    m_data.mine().m_ncov_def_mat.setZero();
    }

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        m_data.mine().m_acov_def_mat.reshapeCol(k,3,2)   = map.jacobian(k);
        gsMatrix<T> acov = m_data.mine().m_acov_def_mat.reshapeCol(k,3,2);

        tmp = acov.transpose() * acov;
        m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = map.deriv2(k);
        deriv2.resize(3,3);
        gsVector<T> normal = m_data.mine().m_normal_def_mat.col(k);

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_data.mine().m_Bcov_def_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_data.mine().m_Bcov_def_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
            gsAsMatrix<T,Dynamic,Dynamic> acon = m_data.mine().m_acon_def_mat.reshapeCol(k,3,2);
            for (index_t i=0; i < 2; i++)
                acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

            // Mixed tensor
            for (index_t i=0; i < 2; i++)
                for (index_t j=0; j < 2; j++)
                    mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

            gsAsMatrix<T,Dynamic,Dynamic> ncov = m_data.mine().m_ncov_def_mat.reshapeCol(k,3,2);
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
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_defpatches->piece(patch)   ).computeMap(map);

    gsMatrix<T,2,2> tmp;

    m_data.mine().m_Acov_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_def_mat.setZero();
    m_data.mine().m_Acon_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_def_mat.setZero();

    m_data.mine().m_acov_def_mat.resize(2*2,map.points.cols());    m_data.mine().m_acov_def_mat.setZero();
    if (basis)
    {
        m_data.mine().m_acon_def_mat.resize(2*2,map.points.cols());    m_data.mine().m_acon_def_mat.setZero();
    }

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_data.mine().m_acov_def_mat.reshapeCol(k,2,2);
        acov = map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
            gsAsMatrix<T,Dynamic,Dynamic> acon = m_data.mine().m_acon_def_mat.reshapeCol(k,2,2);
            for (index_t i=0; i < 2; i++)
                acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    _computeMetricUndeformed_impl<dim>(patch,u,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed_impl(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    gsMatrix<T> deriv2;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_data.mine().m_normal_ori_mat = map.normals;
    m_data.mine().m_normal_ori_mat.colwise().normalized();

    m_data.mine().m_Acov_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_ori_mat.setZero();
    m_data.mine().m_Acon_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_ori_mat.setZero();
    m_data.mine().m_Bcov_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Bcov_ori_mat.setZero();

    m_data.mine().m_acov_ori_mat.resize(2*3,map.points.cols());    m_data.mine().m_acov_ori_mat.setZero();
    if (basis)
    {
        m_data.mine().m_acon_ori_mat.resize(2*3,map.points.cols());    m_data.mine().m_acon_ori_mat.setZero();
        m_data.mine().m_ncov_ori_mat.resize(2*3,map.points.cols());    m_data.mine().m_ncov_ori_mat.setZero();
    }

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        m_data.mine().m_acov_ori_mat.reshapeCol(k,3,2)   = map.jacobian(k);
        gsMatrix<T> acov = m_data.mine().m_acov_ori_mat.reshapeCol(k,3,2);

        tmp = acov.transpose() * acov;
        m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = map.deriv2(k);
        deriv2.resize(3,3);
        gsVector<T> normal = m_data.mine().m_normal_ori_mat.col(k);

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_data.mine().m_Bcov_ori_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_data.mine().m_Bcov_ori_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
           gsAsMatrix<T,Dynamic,Dynamic> acon = m_data.mine().m_acon_ori_mat.reshapeCol(k,3,2);
           for (index_t i=0; i < 2; i++)
               acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

           // Mixed tensor
           for (index_t i=0; i < 2; i++)
               for (index_t j=0; j < 2; j++)
                   mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

           gsAsMatrix<T,Dynamic,Dynamic> ncov = m_data.mine().m_ncov_ori_mat.reshapeCol(k,3,2);
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
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    gsMatrix<T,2,2> tmp;

    m_data.mine().m_Acov_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_ori_mat.setZero();
    m_data.mine().m_Acon_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_ori_mat.setZero();

    m_data.mine().m_acov_ori_mat.resize(2*2,map.points.cols());    m_data.mine().m_acov_ori_mat.setZero();
    if (basis)
    {
        m_data.mine().m_acon_ori_mat.resize(2*2,map.points.cols());    m_data.mine().m_acon_ori_mat.setZero();
    }

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_data.mine().m_acov_ori_mat.reshapeCol(k,2,2);
        acov = map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();
        
        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2);

        // Construct basis
        if (basis)
        {
            gsAsMatrix<T,Dynamic,Dynamic> acon = m_data.mine().m_acon_ori_mat.reshapeCol(k,2,2);
            for (index_t i=0; i < 2; i++)
                acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetric(const index_t k, const T z, const bool basis) const
{
    this->_getMetricDeformed(k,z,basis);
    this->_getMetricUndeformed(k,z,basis);

    T ratio = m_data.mine().m_Gcov_def.determinant() / m_data.mine().m_Gcov_ori.determinant();
    GISMO_ENSURE(ratio >= 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_data.mine().m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_data.mine().m_Gcov_ori.determinant());
    m_data.mine().m_J0_sq = ratio;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed(const index_t k, const T z, const bool basis) const
{
    _getMetricDeformed_impl<dim>(k,z,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed_impl(const index_t k, const T z, const bool basis) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    if (basis)
    {
        GISMO_ENSURE(m_data.mine().m_acov_def_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_acon_def_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_ncov_def_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_normal_def_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_data.mine().m_Acov_def = m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2);
    m_data.mine().m_Acon_def = m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2);
    m_data.mine().m_Bcov_def = m_data.mine().m_Bcov_def_mat.reshapeCol(k,2,2);
    if (basis)
        m_data.mine().m_ncov_def = m_data.mine().m_ncov_def_mat.reshapeCol(k,3,2);

    // Compute full metric
    gsMatrix<T,3,3> Gcov_def;
    Gcov_def.setZero();
    Gcov_def.block(0,0,2,2)= m_data.mine().m_Acov_def - 2.0 * z * m_data.mine().m_Bcov_def + z*z * m_data.mine().m_ncov_def.transpose()*m_data.mine().m_ncov_def;
    Gcov_def(2,2) = 1.0;
    
    m_data.mine().m_Gcov_def = Gcov_def;
    m_data.mine().m_Gcon_def = Gcov_def.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_data.mine().m_acov_def = m_data.mine().m_acov_def_mat.reshapeCol(k,3,2);
    m_data.mine().m_acon_def = m_data.mine().m_acon_def_mat.reshapeCol(k,3,2);
    // g
    m_data.mine().m_gcov_def.leftCols(2) = m_data.mine().m_acov_def + z * m_data.mine().m_ncov_def;
    m_data.mine().m_gcov_ori.col(2) = m_data.mine().m_normal_def_mat.reshapeCol(k,3,1);

    for (index_t c = 0; c!=3; c++)
    {
        m_data.mine().m_gcon_def.col(c) = m_data.mine().m_Gcon_def(c,0) * m_data.mine().m_gcov_def.col(0)
                            + m_data.mine().m_Gcon_def(c,1) * m_data.mine().m_gcov_def.col(1)
                            + m_data.mine().m_Gcon_def(c,2) * m_data.mine().m_gcov_def.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_data.mine().m_Gcov_def_L = m_data.mine().m_Gcov_def;
    // m_data.mine().m_Gcov_def.block(0,0,2,2) -= z*z * m_data.mine().m_ncov_def.transpose()*m_data.mine().m_ncov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed_impl(const index_t k, const T z, const bool basis) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_def.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_def_mat.cols()!=0,"Is the metric initialized?");

    if (basis)
    {
        GISMO_ENSURE(m_data.mine().m_acov_def_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_data.mine().m_Acov_def = m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2);
    m_data.mine().m_Acon_def = m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_data.mine().m_Gcov_def.setZero();
    m_data.mine().m_Gcov_def.block(0,0,2,2)= m_data.mine().m_Acov_def;
    m_data.mine().m_Gcov_def(2,2) = 1.0;
    m_data.mine().m_Gcon_def = m_data.mine().m_Gcov_def.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_data.mine().m_acov_def = m_data.mine().m_acov_def_mat.reshapeCol(k,2,2);
    m_data.mine().m_acon_def = m_data.mine().m_acon_def_mat.reshapeCol(k,2,2);
    // g
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_data.mine().m_gcov_def.setZero();
    m_data.mine().m_gcov_def.block(0,0,2,2) = m_data.mine().m_acov_def;
    m_data.mine().m_gcov_def.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_data.mine().m_gcon_def.col(c) = m_data.mine().m_Gcon_def(c,0) * m_data.mine().m_gcov_def.col(0)
                            + m_data.mine().m_Gcon_def(c,1) * m_data.mine().m_gcov_def.col(1)
                            + m_data.mine().m_Gcon_def(c,2) * m_data.mine().m_gcov_def.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_data.mine().m_Gcov_def_L = m_data.mine().m_Gcov_def;
    // m_data.mine().m_Gcov_def.block(0,0,2,2) -= z*z * m_data.mine().m_ncov_def.transpose()*m_data.mine().m_ncov_def;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed(const index_t k, const T z, const bool basis) const
{
    _getMetricUndeformed_impl<dim>(k,z,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed_impl(const index_t k, const T z, const bool basis) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");

    if (basis)
    {
        GISMO_ENSURE(m_data.mine().m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_normal_ori_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_data.mine().m_Acov_ori = m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2);
    m_data.mine().m_Acon_ori = m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2);
    m_data.mine().m_Bcov_ori = m_data.mine().m_Bcov_ori_mat.reshapeCol(k,2,2);
    if (basis)
        m_data.mine().m_ncov_ori = m_data.mine().m_ncov_ori_mat.reshapeCol(k,3,2);

    // Compute full metric
    m_data.mine().m_Gcov_ori.setZero();
    m_data.mine().m_Gcov_ori.block(0,0,2,2)= m_data.mine().m_Acov_ori - 2.0 * z * m_data.mine().m_Bcov_ori + z*z * m_data.mine().m_ncov_ori.transpose()*m_data.mine().m_ncov_ori;
    m_data.mine().m_Gcov_ori(2,2) = 1.0;
    m_data.mine().m_Gcon_ori = m_data.mine().m_Gcov_ori.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_data.mine().m_acov_ori = m_data.mine().m_acov_ori_mat.reshapeCol(k,3,2);
    m_data.mine().m_acon_ori = m_data.mine().m_acon_ori_mat.reshapeCol(k,3,2);
    // g
    m_data.mine().m_gcov_ori.block(0,0,3,2) = m_data.mine().m_acov_ori + z * m_data.mine().m_ncov_ori;
    m_data.mine().m_gcov_ori.col(2) = m_data.mine().m_normal_ori_mat.reshapeCol(k,3,1);
    for (index_t c = 0; c!=3; c++)
    {
        m_data.mine().m_gcon_ori.col(c) = m_data.mine().m_Gcon_ori(c,0) * m_data.mine().m_gcov_ori.col(0)
                            + m_data.mine().m_Gcon_ori(c,1) * m_data.mine().m_gcov_ori.col(1)
                            + m_data.mine().m_Gcon_ori(c,2) * m_data.mine().m_gcov_ori.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_data.mine().m_Gcov_ori_L = m_data.mine().m_Gcov_ori;
    // m_data.mine().m_Gcov_ori.block(0,0,2,2) -= z*z * m_data.mine().m_ncov_ori.transpose()*m_data.mine().m_ncov_ori;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed_impl(const index_t k, const T z, const bool basis) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");

    if (basis)
    {
        GISMO_ENSURE(m_data.mine().m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
        GISMO_ENSURE(m_data.mine().m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    }

    // metrics
    m_data.mine().m_Acov_ori = m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2);
    m_data.mine().m_Acon_ori = m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_data.mine().m_Gcov_ori.setZero();
    m_data.mine().m_Gcov_ori.block(0,0,2,2)= m_data.mine().m_Acov_ori;
    m_data.mine().m_Gcov_ori(2,2) = 1.0;
    m_data.mine().m_Gcon_ori = m_data.mine().m_Gcov_ori.inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_data.mine().m_acov_ori = m_data.mine().m_acov_ori_mat.reshapeCol(k,2,2);
    m_data.mine().m_acon_ori = m_data.mine().m_acon_ori_mat.reshapeCol(k,2,2);
    // g
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_data.mine().m_gcov_ori.setZero();
    m_data.mine().m_gcov_ori.block(0,0,2,2) = m_data.mine().m_acov_ori;
    m_data.mine().m_gcov_ori.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_data.mine().m_gcon_ori.col(c) = m_data.mine().m_Gcon_ori(c,0) * m_data.mine().m_gcov_ori.col(0)
                            + m_data.mine().m_Gcon_ori(c,1) * m_data.mine().m_gcov_ori.col(1)
                            + m_data.mine().m_Gcon_ori(c,2) * m_data.mine().m_gcov_ori.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_data.mine().m_Gcov_ori_L = m_data.mine().m_Gcov_ori;
    // m_data.mine().m_Gcov_ori.block(0,0,2,2) -= z*z * m_data.mine().m_ncov_ori.transpose()*m_data.mine().m_ncov_ori;
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


    GISMO_ENSURE(m_data.mine().m_gcon_ori.cols()!=0,"Is the basis initialized?");

    gsMatrix<T,3,3> B;
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += C(k,l) * m_data.mine().m_gcon_ori.col(k) * m_data.mine().m_gcon_ori.col(l).transpose();

    eigSolver.compute(B);

    stretchvec.leftCols(2) = eigSolver.eigenvectors().rightCols(2);
    stretchvec.col(2) = m_data.mine().m_gcon_ori.col(2); // replace with: stretchvec.col(0).template head<3>().cross(stretchvec.col(1).template head<3>())
    stretches.block(0,0,2,1) = eigSolver.eigenvalues().block(1,0,2,1); // the eigenvalues are a 3x1 matrix, so we need to use matrix block-operations

    // m_data.mine().m_stretches.at(2) = 1/m_data.mine().m_J0_sq;
    stretches.at(2) = C(2,2);

    for (index_t k=0; k!=3; k++)
        stretches.at(k) = math::sqrt(stretches.at(k));

    result.first = stretches;
    result.second = stretchvec;

    // // DEBUGGING ONLY!
    // gsMatrix<T> ones(3,1);
    // ones.setOnes();
    // gsDebugVar(m_data.mine().m_stretchvec);
    // gsDebugVar(result.first);

    return result;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeStretch(const gsMatrix<T> & C) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalStretch(C);
    m_data.mine().m_stretches = result.first;
    m_data.mine().m_stretchvec = result.second;
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
