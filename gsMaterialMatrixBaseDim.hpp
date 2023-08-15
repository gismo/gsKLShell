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
        _getMetric(i,0.0); // on point i, with height 0.0
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
        _getMetric(i,0.0); // on point i, with height 0.0
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
        _getMetric(i,0.0); // on point i, with height 0.0
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
    _computePoints(patch,u);

    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        _getMetric(i,0.0); // on point i, with height 0.0

        gsAsMatrix<T> C = result.reshapeCol(i,3,3);
        C.setZero();
        C.block(0,0,2,2) = m_data.mine().m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_data.mine().m_J0_sq;
    }
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::_computePoints(const index_t patch, const gsMatrix<T> & u) const
{
    gsMatrix<T> tmp;

    this->_computeMetricUndeformed(patch,u);

    if (Base::m_defpatches->nPieces()!=0)
        this->_computeMetricDeformed(patch,u);

    gsMapData<T> map;
    map.flags = NEED_VALUE;
    map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    m_thickness->piece(patch).eval_into(map.values[0], m_data.mine().m_Tmat);

    m_data.mine().m_parmat.resize(m_pars.size(),map.values[0].cols());
    m_data.mine().m_parmat.setZero();

    for (size_t v=0; v!=m_pars.size(); v++)
    {
        m_pars[v]->piece(patch).eval_into(map.values[0], tmp);
        m_data.mine().m_parmat.row(v) = tmp;
    }

    m_data.mine().m_parvals.resize(m_pars.size());
}

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed(const index_t patch, const gsMatrix<T> & u) const
{
    _computeMetricDeformed_impl<dim>(patch,u);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, const gsMatrix<T> & u) const
{
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_defpatches->piece(patch)   ).computeMap(map);

    gsMatrix<T> deriv2(3,3), mixedB(2,2), acov(3,2), acon(3,2), ncov(3,2), Acov(2,2), Acon(2,2), Bcov(2,2);
    gsMatrix<T> normals;
    gsVector<T> normal;

    normals = map.normals;
    normals.colwise().normalize();
    m_data.mine().m_normal_def_mat = normals;

    m_data.mine().m_Acov_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_def_mat.setZero();
    m_data.mine().m_Acon_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_def_mat.setZero();
    m_data.mine().m_Bcov_def_mat.resize(4,map.points.cols());    m_data.mine().m_Bcov_def_mat.setZero();

    m_data.mine().m_acov_def_mat.resize(2*3,map.points.cols());    m_data.mine().m_acov_def_mat.setZero();
    m_data.mine().m_acon_def_mat.resize(2*3,map.points.cols());    m_data.mine().m_acon_def_mat.setZero();
    m_data.mine().m_ncov_def_mat.resize(2*3,map.points.cols());    m_data.mine().m_ncov_def_mat.setZero();

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2 = map.deriv2(k).reshaped(3,3);
        normal = normals.col(k);

        Bcov(0,0) = deriv2.row(0).dot(normal);
        Bcov(1,1) = deriv2.row(1).dot(normal);
        Bcov(0,1) = Bcov(1,0) = deriv2.row(2).dot(normal);

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = Acon(i,0)*Bcov(0,j) + Acon(i,1)*Bcov(1,j);

        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);

        // Assign members
        m_data.mine().m_acov_def_mat.reshapeCol(k,3,2) = acov;
        m_data.mine().m_acon_def_mat.reshapeCol(k,3,2) = acon;
        m_data.mine().m_ncov_def_mat.reshapeCol(k,3,2) = ncov;
        m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2) = Acov;
        m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2) = Acon;
        m_data.mine().m_Bcov_def_mat.reshapeCol(k,2,2) = Bcov;
    }
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, const gsMatrix<T> & u) const
{
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_defpatches->piece(patch)   ).computeMap(map);

    gsMatrix<T> tmp;
    gsMatrix<T> acov(2,2), acon(2,2), Acov(2,2), Acon(2,2);

    m_data.mine().m_Acov_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_def_mat.setZero();
    m_data.mine().m_Acon_def_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_def_mat.setZero();

    m_data.mine().m_acov_def_mat.resize(2*2,map.points.cols());    m_data.mine().m_acov_def_mat.setZero();
    m_data.mine().m_acon_def_mat.resize(2*2,map.points.cols());    m_data.mine().m_acon_def_mat.setZero();

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Assign members
        m_data.mine().m_acov_def_mat.reshapeCol(k,2,2) = acov;
        m_data.mine().m_acon_def_mat.reshapeCol(k,2,2) = acon;
        m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2) = Acov;
        m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2) = Acon;
    }
}


//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u) const
{
    _computeMetricUndeformed_impl<dim>(patch,u);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed_impl(const index_t patch, const gsMatrix<T> & u) const
{
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    gsMatrix<T> deriv2(3,3), mixedB(2,2), acov(3,2), acon(3,2), ncov(3,2), Acov(2,2), Acon(2,2), Bcov(2,2);
    gsMatrix<T> normals;
    gsVector<T> normal;

    normals = map.normals;
    normals.colwise().normalize();
    m_data.mine().m_normal_ori_mat = normals;

    m_data.mine().m_Acov_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_ori_mat.setZero();
    m_data.mine().m_Acon_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_ori_mat.setZero();
    m_data.mine().m_Bcov_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Bcov_ori_mat.setZero();

    m_data.mine().m_acov_ori_mat.resize(2*3,map.points.cols());    m_data.mine().m_acov_ori_mat.setZero();
    m_data.mine().m_acon_ori_mat.resize(2*3,map.points.cols());    m_data.mine().m_acon_ori_mat.setZero();
    m_data.mine().m_ncov_ori_mat.resize(2*3,map.points.cols());    m_data.mine().m_ncov_ori_mat.setZero();

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2 = map.deriv2(k).reshaped(3,3);
        normal = normals.col(k);

        Bcov(0,0) = deriv2.row(0).dot(normal);
        Bcov(1,1) = deriv2.row(1).dot(normal);
        Bcov(0,1) = Bcov(1,0) = deriv2.row(2).dot(normal);

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = Acon(i,0)*Bcov(0,j) + Acon(i,1)*Bcov(1,j);

        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);

        // Assign members
        m_data.mine().m_acov_ori_mat.reshapeCol(k,3,2) = acov;
        m_data.mine().m_acon_ori_mat.reshapeCol(k,3,2) = acon;
        m_data.mine().m_ncov_ori_mat.reshapeCol(k,3,2) = ncov;
        m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2) = Acov;
        m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2) = Acon;
        m_data.mine().m_Bcov_ori_mat.reshapeCol(k,2,2) = Bcov;
    }
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed_impl(const index_t patch, const gsMatrix<T> & u) const
{
    gsMapData<T> map;
    map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    map.points = u;
    dynamic_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map);

    gsMatrix<T> tmp;
    gsMatrix<T> acov(2,2), acon(2,2), Acov(2,2), Acon(2,2);

    m_data.mine().m_Acov_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acov_ori_mat.setZero();
    m_data.mine().m_Acon_ori_mat.resize(4,map.points.cols());    m_data.mine().m_Acon_ori_mat.setZero();

    m_data.mine().m_acov_ori_mat.resize(2*2,map.points.cols());    m_data.mine().m_acov_ori_mat.setZero();
    m_data.mine().m_acon_ori_mat.resize(2*2,map.points.cols());    m_data.mine().m_acon_ori_mat.setZero();

    for (index_t k=0; k!= map.points.cols(); k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Assign members
        m_data.mine().m_acov_ori_mat.reshapeCol(k,2,2) = acov;
        m_data.mine().m_acon_ori_mat.reshapeCol(k,2,2) = acon;
        m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2) = Acov;
        m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2) = Acon;
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetric(const index_t k, const T z) const
{
    this->_getMetricDeformed(k,z);
    this->_getMetricUndeformed(k,z);

    T ratio;
    T det_ori = m_data.mine().m_Gcov_ori.determinant();
    T det_def = m_data.mine().m_Gcov_def.determinant();

    if (det_ori==0 && det_def==0 || math::isnan(det_ori) || math::isnan(det_def))
    {
        gsWarn<<"Jacobian determinant is undefined: J^2 = det(Gcov_def) / det(Gcov_ori) = "<<det_def<<"/"<<det_ori<<"! J^2 is set to 1";
        ratio = 1;
    }
    else
        ratio = det_def / det_ori;

    GISMO_ENSURE(ratio >= 0, "Jacobian determinant is negative! det(Gcov_def) = "<<det_def<<"; det(Gcov_ori) = "<<det_ori);
    m_data.mine().m_J0_sq = ratio;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed(const index_t k, const T z) const
{
    _getMetricDeformed_impl<dim>(k,z);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed_impl(const index_t k, const T z) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_ncov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_normal_def_mat.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> Acov_def, Acon_def, Bcov_def, ncov_def, Gcov_def(3,3), Gcon_def(3,3),
                acov_def, acon_def, normal(3,1), gcov_def(3,3), gcon_def(3,3);
    // Get metric information
    Acov_def = m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2);
    Acon_def = m_data.mine().m_Acon_def_mat.reshapeCol(k,2,2);
    Bcov_def = m_data.mine().m_Bcov_def_mat.reshapeCol(k,2,2);
    ncov_def = m_data.mine().m_ncov_def_mat.reshapeCol(k,3,2);

    // Compute full metric
    Gcov_def.setZero();
    Gcov_def.block(0,0,2,2)= Acov_def - 2.0 * z * Bcov_def + z*z * ncov_def.transpose()*ncov_def;
    Gcov_def(2,2) = 1.0;
    Gcon_def = Gcov_def.inverse();

    // Assign members
    m_data.mine().m_Acov_def = Acov_def;
    m_data.mine().m_Acon_def = Acon_def;
    m_data.mine().m_Bcov_def = Bcov_def;
    m_data.mine().m_ncov_def = ncov_def;
    m_data.mine().m_Gcov_def = Gcov_def;
    m_data.mine().m_Gcon_def = Gcon_def;

    // Get basis vectors
    acov_def = m_data.mine().m_acov_def_mat.reshapeCol(k,3,2);
    acon_def = m_data.mine().m_acon_def_mat.reshapeCol(k,3,2);
    normal   = m_data.mine().m_normal_def_mat.reshapeCol(k,3,1);

    // Compute g_cov
    gcov_def.setZero();
    gcov_def.leftCols(2) = acov_def + z * ncov_def;
    gcov_def.col(2) = normal;

    // Compute g_con
    for (index_t c = 0; c!=3; c++)
        gcon_def.col(c) = Gcon_def(c,0) * gcov_def.col(0) + Gcon_def(c,1) * gcov_def.col(1) + Gcon_def(c,2) * gcov_def.col(2);

    // Assign members
    m_data.mine().m_acov_def = acov_def;
    m_data.mine().m_acon_def = acov_def;
    m_data.mine().m_gcov_def = gcov_def;
    m_data.mine().m_gcon_def = gcon_def;

    // // create a Linearised covariant tensor for SvK models
    // m_data.mine().m_Gcov_def_L = m_data.mine().m_Gcov_def;
    // m_data.mine().m_Gcov_def.block(0,0,2,2) -= z*z * m_data.mine().m_ncov_def.transpose()*m_data.mine().m_ncov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed_impl(const index_t k, const T z) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_acon_def_mat.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> Acov_def, Acon_def, Gcov_def(3,3), Gcon_def(3,3),
                acov_def, acon_def, normal(3,1), gcov_def(3,3), gcon_def(3,3);

    // Get metric information
    Acov_def = m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2);
    Acon_def = m_data.mine().m_Acov_def_mat.reshapeCol(k,2,2);
    
    // Compute full metric
    Gcov_def.setZero();
    Gcov_def.block(0,0,2,2)= Acov_def;;
    Gcov_def(2,2) = 1.0;
    Gcon_def = Gcov_def.inverse();

    // Assign members
    m_data.mine().m_Acov_def = Acov_def;
    m_data.mine().m_Acon_def = Acon_def;
    m_data.mine().m_Gcov_def = Gcov_def;
    m_data.mine().m_Gcon_def = Gcon_def;

    // Get basis vectors
    acov_def = m_data.mine().m_acov_def_mat.reshapeCol(k,2,2);
    acon_def = m_data.mine().m_acon_def_mat.reshapeCol(k,2,2);
    normal << 0,0,1;

    // Compute g_cov
    gcov_def.setZero();
    gcov_def.block(0,0,2,2) = acov_def;
    gcov_def.col(2) = normal;

    // Compute g_con
    for (index_t c = 0; c!=3; c++)
        gcon_def.col(c) = Gcon_def(c,0) * gcov_def.col(0) + Gcon_def(c,1) * gcov_def.col(1) + Gcon_def(c,2) * gcov_def.col(2);

    // Assign members
    m_data.mine().m_acov_def = acov_def;
    m_data.mine().m_acon_def = acon_def;
    m_data.mine().m_gcov_def = gcov_def;
    m_data.mine().m_gcon_def = gcon_def;

    // // create a Linearised covariant tensor for SvK models
    // m_data.mine().m_Gcov_def_L = m_data.mine().m_Gcov_def;
    // m_data.mine().m_Gcov_def.block(0,0,2,2) -= z*z * m_data.mine().m_ncov_def.transpose()*m_data.mine().m_ncov_def;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed(const index_t k, const T z) const
{
    _getMetricUndeformed_impl<dim>(k,z);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed_impl(const index_t k, const T z) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_normal_ori_mat.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> Acov_ori, Acon_ori, Bcov_ori, ncov_ori, Gcov_ori(3,3), Gcon_ori(3,3),
                acov_ori, acon_ori, normal(3,1), gcov_ori(3,3), gcon_ori(3,3);

    // Get metric information
    Acov_ori = m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2);
    Acon_ori = m_data.mine().m_Acon_ori_mat.reshapeCol(k,2,2);
    Bcov_ori = m_data.mine().m_Bcov_ori_mat.reshapeCol(k,2,2);
    ncov_ori = m_data.mine().m_ncov_ori_mat.reshapeCol(k,3,2);

    // Compute full metric
    Gcov_ori.setZero();
    Gcov_ori.block(0,0,2,2)= Acov_ori - 2.0 * z * Bcov_ori + z*z * ncov_ori.transpose()*ncov_ori;
    Gcov_ori(2,2) = 1.0;
    Gcon_ori = Gcov_ori.inverse();

    // Assign members
    m_data.mine().m_Acov_ori = Acov_ori;
    m_data.mine().m_Acon_ori = Acon_ori;
    m_data.mine().m_Bcov_ori = Bcov_ori;
    m_data.mine().m_ncov_ori = ncov_ori;
    m_data.mine().m_Gcov_ori = Gcov_ori;
    m_data.mine().m_Gcon_ori = Gcon_ori;

    // Get basis vectors
    acov_ori = m_data.mine().m_acov_ori_mat.reshapeCol(k,3,2);
    acon_ori = m_data.mine().m_acon_ori_mat.reshapeCol(k,3,2);
    normal   = m_data.mine().m_normal_ori_mat.reshapeCol(k,3,1);

    // Compute g_cov
    gcov_ori.setZero();
    gcov_ori.leftCols(2) = acov_ori + z * ncov_ori;
    gcov_ori.col(2) = normal;

    // Compute g_con
    for (index_t c = 0; c!=3; c++)
        gcon_ori.col(c) = Gcon_ori(c,0) * gcov_ori.col(0) + Gcon_ori(c,1) * gcov_ori.col(1) + Gcon_ori(c,2) * gcov_ori.col(2);

    // Assign members
    m_data.mine().m_acov_ori = acov_ori;
    m_data.mine().m_acon_ori = acov_ori;
    m_data.mine().m_gcov_ori = gcov_ori;
    m_data.mine().m_gcon_ori = gcon_ori;

    // // create a Linearised covariant tensor for SvK models
    // m_data.mine().m_Gcov_ori_L = m_data.mine().m_Gcov_ori;
    // m_data.mine().m_Gcov_ori.block(0,0,2,2) -= z*z * m_data.mine().m_ncov_ori.transpose()*m_data.mine().m_ncov_ori;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixBaseDim<dim,T>::_getMetricUndeformed_impl(const index_t k, const T z) const
{
    GISMO_ENSURE(m_data.mine().m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_data.mine().m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_data.mine().m_acon_ori_mat.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> Acov_ori, Acon_ori, Gcov_ori(3,3), Gcon_ori(3,3),
                acov_ori, acon_ori, normal(3,1), gcov_ori(3,3), gcon_ori(3,3);

    // Get metric information
    Acov_ori = m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2);
    Acon_ori = m_data.mine().m_Acov_ori_mat.reshapeCol(k,2,2);
    
    // Compute full metric
    Gcov_ori.setZero();
    Gcov_ori.block(0,0,2,2)= Acov_ori;;
    Gcov_ori(2,2) = 1.0;
    Gcon_ori = Gcov_ori.inverse();

    // Assign members
    m_data.mine().m_Acov_ori = Acov_ori;
    m_data.mine().m_Acon_ori = Acon_ori;
    m_data.mine().m_Gcov_ori = Gcov_ori;
    m_data.mine().m_Gcon_ori = Gcon_ori;

    // Get basis vectors
    acov_ori = m_data.mine().m_acov_ori_mat.reshapeCol(k,2,2);
    acon_ori = m_data.mine().m_acon_ori_mat.reshapeCol(k,2,2);
    normal << 0,0,1;

    // Compute g_cov
    gcov_ori.setZero();
    gcov_ori.block(0,0,2,2) = acov_ori;
    gcov_ori.col(2) = normal;

    // Compute g_con
    for (index_t c = 0; c!=3; c++)
        gcon_ori.col(c) = Gcon_ori(c,0) * gcov_ori.col(0) + Gcon_ori(c,1) * gcov_ori.col(1) + Gcon_ori(c,2) * gcov_ori.col(2);

    // Assign members
    m_data.mine().m_acov_ori = acov_ori;
    m_data.mine().m_acon_ori = acon_ori;
    m_data.mine().m_gcov_ori = gcov_ori;
    m_data.mine().m_gcon_ori = gcon_ori;

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

    gsEigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;


    GISMO_ENSURE(m_data.mine().m_gcon_ori.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> B(3,3);
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
