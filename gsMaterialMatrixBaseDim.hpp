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

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map_ori.mine().flags = NEED_VALUE;
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);

    result.resize(1, u.cols());
    m_thickness->piece(patch).eval_into(m_map_ori.mine().values[0], m_Tmat);
    m_density->piece(patch).eval_into(m_map_ori.mine().values[0], m_rhomat);
    for (index_t i = 0; i != u.cols(); ++i) // points
    {
        result(0,i) = m_Tmat(0,i)*m_rhomat(0,i);
    }

}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map_ori.mine().flags = NEED_VALUE;
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);
    m_thickness->piece(patch).eval_into(m_map_ori.mine().values[0], result);
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
        res = this->_evalStretch(C,m_gcon_ori);
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
        res = this->_evalStretch(C,m_gcon_ori);
        result.col(i) = res.second.reshape(9,1);
    }
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) covariant basis
template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::spec2cov_transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.resize(9, u.cols());
    gsMatrix<T> tmp, covbasis,sbasis;
    this->stretchDir_into(patch,u,tmp);
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0,true); // on point i, with height 0.0
        sbasis = tmp.reshapeCol(i,3,3);
        covbasis = m_gcov_ori;
        result.col(i) = this->_transformation(covbasis,sbasis).reshape(9,1);
    }
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::spec2con_transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
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

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) covariant basis
template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::cov2cart_transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.resize(9, u.cols());
    gsMatrix<T> tmp, covbasis,cartbasis(3,3);
    cartbasis.setIdentity();
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0,true); // on point i, with height 0.0
        covbasis = m_gcov_ori;
        result.col(i) = this->_transformation(cartbasis,covbasis).reshape(9,1);
    }
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::con2cart_transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.resize(9, u.cols());
    gsMatrix<T> tmp, conbasis,cartbasis(3,3);
    cartbasis.setIdentity();
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0,true); // on point i, with height 0.0
        conbasis = m_gcon_ori;
        result.col(i) = this->_transformation(cartbasis,conbasis).reshape(9,1);
    }
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::parameters_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map_ori.mine().flags = NEED_VALUE;
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);

    gsMatrix<T> tmp;
    result.resize(m_pars.size(),m_map_ori.mine().values[0].cols());
    result.setZero();
    for (size_t v=0; v!=m_pars.size(); v++)
    {
        m_pars[v]->piece(patch).eval_into(m_map_ori.mine().values[0], tmp);
        result.row(v) = tmp;
    }
}

// Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::deformation_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);

    _computePoints(patch,u,true);

    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        _getMetric(i,0.0,true); // on point i, with height 0.0

        gsAsMatrix<T> C = result.reshapeCol(i,3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_J0_sq;
    }
}

template <short_t dim, class T >
void gsMaterialMatrixBaseDim<dim,T>::_computePoints(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    gsMatrix<T> tmp;

    this->_computeMetricUndeformed(patch,u,basis);
    if (Base::m_defpatches->nPieces()!=0)
        this->_computeMetricDeformed(patch,u,basis);

    m_thickness->piece(patch).eval_into(m_map_ori.mine().values[0], m_Tmat);

    m_parmat.resize(m_pars.size(),m_map_ori.mine().values[0].cols());
    m_parmat.setZero();

    for (size_t v=0; v!=m_pars.size(); v++)
    {
        m_pars[v]->piece(patch).eval_into(m_map_ori.mine().values[0], tmp);
        m_parmat.row(v) = tmp;
    }

    m_parvals.resize(m_pars.size());
}

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    m_map_def.mine().flags = m_map_ori.mine().flags;
    m_map_def.mine().points = u;
    static_cast<const gsFunction<T>&>(Base::m_defpatches->piece(patch)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    _computeMetricDeformed_impl<dim>(patch,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, bool basis) const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acon_def_mat.setZero();
    m_Bcov_def_mat.resize(4,m_map_def.mine().points.cols());    m_Bcov_def_mat.setZero();

    m_acov_def_mat.resize(2*3,m_map_def.mine().points.cols());    m_acov_def_mat.setZero();
    if (basis)
    {
        m_acon_def_mat.resize(2*3,m_map_def.mine().points.cols());    m_acon_def_mat.setZero();
        m_ncov_def_mat.resize(2*3,m_map_def.mine().points.cols());    m_ncov_def_mat.setZero();
    }

    for (index_t k=0; k!= m_map_def.mine().points.cols(); k++)
    {
        m_acov_def_mat.reshapeCol(k,3,2)   = m_map_def.mine().jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,3,2);
        acov = m_map_def.mine().jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map_def.mine().deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map_def.mine().normal(k).normalized();

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
gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, bool basis) const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acon_def_mat.setZero();

    m_acov_def_mat.resize(2*2,m_map_def.mine().points.cols());    m_acov_def_mat.setZero();
    if (basis)
    {
        m_acon_def_mat.resize(2*2,m_map_def.mine().points.cols());    m_acon_def_mat.setZero();
    }

    for (index_t k=0; k!= m_map_def.mine().points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,2,2);
        acov = m_map_def.mine().jacobian(k);

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


// template <short_t dim, class T>
// template <short_t _dim>
// typename std::enable_if<_dim==3 &&  _TFT, void>::type
// gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, bool basis) const
// {
//     GISMO_ASSERT(m_map_def.mine().points.cols() == m_thetas.size(),"The number of thetas is not equal to the number of evaluation points, points.size() != thetas.size() ("<<m_map_def.mine().points.cols()<<"!="<<m_thetas.size()<<")");
//     GISMO_ASSERT(m_map_def.mine().points.cols() == m_gammas.size(),"The number of gammas is not equal to the number of evaluation points, points.size() != gammas.size() ("<<m_map_def.mine().points.cols()<<"!="<<m_gammas.size()<<")");

//     gsMatrix<T> deriv2;
//     gsMatrix<T,3,1> normal;
//     gsMatrix<T,2,2> mixedB;
//     gsMatrix<T,2,2> tmp, n_mat;

//     m_Acov_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acov_def_mat.setZero();
//     m_Acon_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acon_def_mat.setZero();
//     m_Bcov_def_mat.resize(4,m_map_def.mine().points.cols());    m_Bcov_def_mat.setZero();

//     m_acov_def_mat.resize(2*3,m_map_def.mine().points.cols());    m_acov_def_mat.setZero();
//     if (basis)
//     {
//         m_acon_def_mat.resize(2*3,m_map_def.mine().points.cols());    m_acon_def_mat.setZero();
//         m_ncov_def_mat.resize(2*3,m_map_def.mine().points.cols());    m_ncov_def_mat.setZero();
//     }

//     for (index_t k=0; k!= m_map_def.mine().points.cols(); k++)
//     {
//         m_acov_def_mat.reshapeCol(k,3,2)   = m_map_def.mine().jacobian(k);
//         gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,3,2);
//         acov = m_map_def.mine().jacobian(k);

//         n_mat<< math::sin(m_thetas.at(k))*math::sin(m_thetas.at(k)), math::sin(m_thetas.at(k))*math::cos(m_thetas.at(k)),
//                 math::cos(m_thetas.at(k))*math::sin(m_thetas.at(k)), math::cos(m_thetas.at(k))*math::cos(m_thetas.at(k));

//         tmp = acov.transpose() * acov + 2*m_gammas.at(k) * n_mat; // this line is the big change!
//         m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
//         m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

//         gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

//         // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
//         deriv2    = m_map_def.mine().deriv2(k);
//         deriv2.resize(3,3);
//         normal    = m_map_def.mine().normal(k).normalized();

//         tmp(0,0) = deriv2.row(0).dot(normal);
//         tmp(1,1) = deriv2.row(1).dot(normal);
//         tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

//         m_Bcov_def_mat.reshapeCol(k,2,2) = tmp;
//         gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_def_mat.reshapeCol(k,2,2);

//         // Construct basis
//         if (basis)
//         {
//             gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,3,2);
//             for (index_t i=0; i < 2; i++)
//                 acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

//             // Mixed tensor
//             for (index_t i=0; i < 2; i++)
//                 for (index_t j=0; j < 2; j++)
//                     mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

//             gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_def_mat.reshapeCol(k,3,2);
//             for (index_t i=0; i < 2; i++)
//                 ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);
//         }
//     }
// }

// template <short_t dim, class T>
// template <short_t _dim>
// typename std::enable_if<_dim==2 && _TFT, void>::type
// gsMaterialMatrixBaseDim<dim,T>::_computeMetricDeformed_impl(const index_t patch, bool basis) const
// {
//     GISMO_ASSERT(m_map_def.mine().points.cols() == m_thetas.size(),"The number of thetas is not equal to the number of evaluation points, points.size() != thetas.size() ("<<m_map_def.mine().points.cols()<<"!="<<m_thetas.size()<<")");
//     GISMO_ASSERT(m_map_def.mine().points.cols() == m_gammas.size(),"The number of gammas is not equal to the number of evaluation points, points.size() != gammas.size() ("<<m_map_def.mine().points.cols()<<"!="<<m_gammas.size()<<")");

//     gsMatrix<T,2,2> tmp, n_mat;

//     m_Acov_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acov_def_mat.setZero();
//     m_Acon_def_mat.resize(4,m_map_def.mine().points.cols());    m_Acon_def_mat.setZero();

//     m_acov_def_mat.resize(2*2,m_map_def.mine().points.cols());    m_acov_def_mat.setZero();
//     if (basis)
//     {
//         m_acon_def_mat.resize(2*2,m_map_def.mine().points.cols());    m_acon_def_mat.setZero();
//     }

//     for (index_t k=0; k!= m_map_def.mine().points.cols(); k++)
//     {
//         gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,2,2);
//         acov = m_map_def.mine().jacobian(k);

//         n_mat<< math::cos(m_thetas.at(k))*math::cos(m_thetas.at(k)), math::sin(m_thetas.at(k))*math::cos(m_thetas.at(k)),
//                 math::cos(m_thetas.at(k))*math::sin(m_thetas.at(k)), math::sin(m_thetas.at(k))*math::sin(m_thetas.at(k));

//         tmp = acov.transpose() * acov + 2*m_gammas.at(k) * n_mat; // this line is the big change!
//         m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
//         m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

//         gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

//         // Construct basis
//         if (basis)
//         {
//             gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,2,2);
//             for (index_t i=0; i < 2; i++)
//                 acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
//         }
//     }
// }

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed(const index_t patch, const gsMatrix<T> & u, bool basis) const
{
    m_map_ori.mine().flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    _computeMetricUndeformed_impl<dim>(patch,basis);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed_impl(const index_t patch, bool basis) const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map_ori.mine().points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map_ori.mine().points.cols());    m_Acon_ori_mat.setZero();
    m_Bcov_ori_mat.resize(4,m_map_ori.mine().points.cols());    m_Bcov_ori_mat.setZero();

    m_acov_ori_mat.resize(2*3,m_map_ori.mine().points.cols());    m_acov_ori_mat.setZero();
    if (basis)
    {
        m_acon_ori_mat.resize(2*3,m_map_ori.mine().points.cols());    m_acon_ori_mat.setZero();
        m_ncov_ori_mat.resize(2*3,m_map_ori.mine().points.cols());    m_ncov_ori_mat.setZero();
    }

    for (index_t k=0; k!= m_map_ori.mine().points.cols(); k++)
    {
        m_acov_ori_mat.reshapeCol(k,3,2)   = m_map_ori.mine().jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,3,2);
        acov = m_map_ori.mine().jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map_ori.mine().deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map_ori.mine().normal(k).normalized();

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
gsMaterialMatrixBaseDim<dim,T>::_computeMetricUndeformed_impl(const index_t patch, const bool basis) const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map_ori.mine().points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map_ori.mine().points.cols());    m_Acon_ori_mat.setZero();

    m_acov_ori_mat.resize(2*2,m_map_ori.mine().points.cols());    m_acov_ori_mat.setZero();
    if (basis)
    {
        m_acon_ori_mat.resize(2*2,m_map_ori.mine().points.cols());    m_acon_ori_mat.setZero();
    }

    for (index_t k=0; k!= m_map_ori.mine().points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,2,2);
        acov = m_map_ori.mine().jacobian(k);

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
gsMatrix<T,2,2> gsMaterialMatrixBaseDim<dim,T>::_getAcov_def(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    return m_Acov_def_mat.reshapeCol(k,2,2);
}

template <short_t dim, class T>
gsMatrix<T,2,2> gsMaterialMatrixBaseDim<dim,T>::_getAcon_def(index_t k, T z) const
{
    GISMO_ENSURE(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    return m_Acon_def_mat.reshapeCol(k,2,2);
}

template <short_t dim, class T>
gsMatrix<T,2,2> gsMaterialMatrixBaseDim<dim,T>::_getBcov_def(index_t k, T z) const
{
    return _getBcov_def_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,2> gsMaterialMatrixBaseDim<dim,T>::_getncov_def(index_t k, T z) const
{
    return _getncov_def_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getGcov_def(index_t k, T z) const
{
    return _getGcov_def_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getGcon_def(index_t k, T z) const
{
    gsMatrix<T,3,3> Gcov_def = _getGcon_def(k,z);
    return Gcov_def.inverse();
}

template <short_t dim, class T>
gsMatrix<T,3,2> gsMaterialMatrixBaseDim<dim,T>::_getacov_def(index_t k, T z) const
{
    GISMO_ENSURE(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    return m_acov_def_mat.reshapeCol(k,3,2);
}

template <short_t dim, class T>
gsMatrix<T,3,2> gsMaterialMatrixBaseDim<dim,T>::_getacon_def(index_t k, T z) const
{
    GISMO_ENSURE(m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    return m_acon_def_mat.reshapeCol(k,3,2);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getgcov_def(index_t k, T z) const
{
    return _getgcov_def_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getgcon_def(index_t k, T z) const
{
    gsMatrix<T,3,3> gcov_def = _getgcov_def(k,z);
    gsMatrix<T,3,3> Gcon_def = _getGcon_def(k,z);
    gsMatrix<T,3,3> gcon_def;
    for (index_t c = 0; c!=3; c++)
    {
        gcon_def.col(c) =   Gcon_def(c,0) * gcov_def.col(0)
                          + Gcon_def(c,1) * gcov_def.col(1)
                          + Gcon_def(c,2) * gcov_def.col(2);
    }
    return gcon_def;
}

template <short_t dim, class T>
gsMatrix<T,2,2> gsMaterialMatrixBaseDim<dim,T>::_getAcov_ori(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    return m_Acov_ori_mat.reshapeCol(k,2,2);
}

template <short_t dim, class T>
gsMatrix<T,2,2> gsMaterialMatrixBaseDim<dim,T>::_getAcon_ori(index_t k, T z) const
{
    GISMO_ENSURE(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    return m_Acon_ori_mat.reshapeCol(k,2,2);
}

template <short_t dim, class T>
gsMatrix<T,2,2> gsMaterialMatrixBaseDim<dim,T>::_getBcov_ori(index_t k, T z) const
{
    return _getBcov_ori_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,2> gsMaterialMatrixBaseDim<dim,T>::_getncov_ori(index_t k, T z) const
{
    return _getncov_ori_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getGcov_ori(index_t k, T z) const
{
    return _getGcov_ori_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getGcon_ori(index_t k, T z) const
{
    gsMatrix<T,3,3> Gcov_ori = _getGcon_ori(k,z);
    return Gcov_ori.inverse();
}

template <short_t dim, class T>
gsMatrix<T,3,2> gsMaterialMatrixBaseDim<dim,T>::_getacov_ori(index_t k, T z) const
{
    GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    return m_acov_ori_mat.reshapeCol(k,3,2);
}

template <short_t dim, class T>
gsMatrix<T,3,2> gsMaterialMatrixBaseDim<dim,T>::_getacon_ori(index_t k, T z) const
{
    GISMO_ENSURE(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    return m_acon_ori_mat.reshapeCol(k,3,2);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getgcov_ori(index_t k, T z) const
{
    return _getgcov_ori_impl<dim>(k,z);
}

template <short_t dim, class T>
gsMatrix<T,3,3> gsMaterialMatrixBaseDim<dim,T>::_getgcon_ori(index_t k, T z) const
{
    gsMatrix<T,3,3> gcov_ori = _getgcov_ori(k,z);
    gsMatrix<T,3,3> Gcon_ori = _getGcon_ori(k,z);
    gsMatrix<T,3,3> gcon_ori;
    for (index_t c = 0; c!=3; c++)
    {
        gcon_ori.col(c) =   Gcon_ori(c,0) * gcov_ori.col(0)
                          + Gcon_ori(c,1) * gcov_ori.col(1)
                          + Gcon_ori(c,2) * gcov_ori.col(2);
    }
    return gcon_ori;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,2,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getBcov_def_impl(index_t k, T z) const
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,2,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getBcov_def_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    return m_Bcov_def_mat.reshapeCol(k,2,2);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,3,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getncov_def_impl(index_t k, T z) const
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,3,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getncov_def_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_ncov_def_mat.cols()!=0,"Is the basis initialized?");
    return m_ncov_def_mat.reshapeCol(k,3,2);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getGcov_def_impl(index_t k, T z) const
{
    // metrics
    gsMatrix<T> Acov_def = _getAcov_def(k,z);

    // Compute full metric
    gsMatrix<T,3,3> Gcov_def;
    Gcov_def.setZero();
    Gcov_def.block(0,0,2,2)= Acov_def;
    Gcov_def(2,2) = 1.0;
    return Gcov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getGcov_def_impl(index_t k, T z) const
{
    // metrics
    gsMatrix<T> Acov_def = _getAcov_def(k,z);
    gsMatrix<T> Bcov_def = _getBcov_def(k,z);
    gsMatrix<T> ncov_def = _getncov_def(k,z);

    // Compute full metric
    gsMatrix<T,3,3> Gcov_def;
    Gcov_def.setZero();
    Gcov_def.block(0,0,2,2)= Acov_def - 2.0 * z * m_Bcov_def + z*z * ncov_def.transpose()*ncov_def;
    Gcov_def(2,2) = 1.0;
    return Gcov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getgcov_def_impl(index_t k, T z) const
{
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    gsMatrix<T,3,3> gcov_def;
    gcov_def.setZero();
    gcov_def.block(0,0,2,2) = _getacov_def(k,z);
    gcov_def.col(2) = normal;
    return gcov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getgcov_def_impl(index_t k, T z) const
{
    gsMatrix<T,3,1> normal = m_map_def.mine().normal(k).normalized();
    gsMatrix<T,3,2> acov_def = _getacov_def(k,z);
    gsMatrix<T,3,2> ncov_def = _getncov_def(k,z);
    gsMatrix<T,3,3> gcov_def;
    gcov_def.setZero();
    gcov_def.leftCols(2) = acov_def + z * ncov_def;
    gcov_def.col(2) = normal;
    return gcov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,2,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getBcov_ori_impl(index_t k, T z) const
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,2,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getBcov_ori_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");
    return m_Bcov_ori_mat.reshapeCol(k,2,2);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,3,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getncov_ori_impl(index_t k, T z) const
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,3,2>>::type
gsMaterialMatrixBaseDim<dim,T>::_getncov_ori_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");
    return m_ncov_ori_mat.reshapeCol(k,3,2);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getGcov_ori_impl(index_t k, T z) const
{
    // metrics
    gsMatrix<T> Acov_ori = _getAcov_ori(k,z);

    // Compute full metric
    gsMatrix<T,3,3> Gcov_ori;
    Gcov_ori.setZero();
    Gcov_ori.block(0,0,2,2)= Acov_ori;
    Gcov_ori(2,2) = 1.0;
    return Gcov_ori;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getGcov_ori_impl(index_t k, T z) const
{
    // metrics
    gsMatrix<T> Acov_ori = _getAcov_ori(k,z);
    gsMatrix<T> Bcov_ori = _getBcov_ori(k,z);
    gsMatrix<T> ncov_ori = _getncov_ori(k,z);

    // Compute full metric
    gsMatrix<T,3,3> Gcov_ori;
    Gcov_ori.setZero();
    Gcov_ori.block(0,0,2,2)= Acov_ori - 2.0 * z * m_Bcov_ori + z*z * ncov_ori.transpose()*ncov_ori;
    Gcov_ori(2,2) = 1.0;
    return Gcov_ori;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getgcov_ori_impl(index_t k, T z) const
{
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    gsMatrix<T,3,3> gcov_ori;
    gcov_ori.setZero();
    gcov_ori.block(0,0,2,2) = _getacov_ori(k,z);
    gcov_ori.col(2) = normal;
    return gcov_ori;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, gsMatrix<T,3,3>>::type
gsMaterialMatrixBaseDim<dim,T>::_getgcov_ori_impl(index_t k, T z) const
{
    gsMatrix<T,3,1> normal = m_map_ori.mine().normal(k).normalized();
    gsMatrix<T,3,2> acov_ori = _getacov_ori(k,z);
    gsMatrix<T,3,2> ncov_ori = _getncov_ori(k,z);
    gsMatrix<T,3,3> gcov_ori;
    gcov_ori.setZero();
    gcov_ori.leftCols(2) = acov_ori + z * ncov_ori;
    gcov_ori.col(2) = normal;
    return gcov_ori;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetric(index_t k, T z, const gsMatrix<T> & C, bool basis) const
{
    this->_getMetricDeformed(C);
    this->_getMetricUndeformed(k,z,basis);

    T ratio = m_Gcov_def.determinant() / m_Gcov_ori.determinant();
    GISMO_ENSURE(ratio >= 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_Gcov_ori.determinant());
    m_J0_sq = ratio;
}

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetric(index_t k, T z, bool basis) const
{
    this->_getMetricDeformed(k,z,basis);
    this->_getMetricUndeformed(k,z,basis);

    T ratio = m_Gcov_def.determinant() / m_Gcov_ori.determinant();
    GISMO_ENSURE(ratio >= 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_Gcov_ori.determinant());
    m_J0_sq = ratio;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_getMetricDeformed(const gsMatrix<T> & C) const
{
    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def(0,0) = C(0,0);
    m_Gcov_def(1,1) = C(1,0);
    m_Gcov_def(0,1) = m_Gcov_def(1,0) = C(2,0);
    m_Gcov_def(2,2) =
    m_Gcon_def(2,2) = 1.0;
    m_Gcon_def.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2).inverse();
    // CLear other members on the deformed geometry
    m_Acov_def.setZero();
    m_Acon_def.setZero();
    m_Bcov_def.setZero();
    m_Bcon_def.setZero();
    m_acov_def.setZero();
    m_acon_def.setZero();
    m_gcov_def.setZero();
    m_gcon_def.setZero();
}

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
    if (basis)
        m_ncov_def = m_ncov_def_mat.reshapeCol(k,3,2);

    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def - 2.0 * z * m_Bcov_def + z*z * m_ncov_def.transpose()*m_ncov_def;
    m_Gcon_def(2,2) = m_Gcov_def(2,2) = 1.0;
    m_Gcon_def.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2).inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_acov_def = m_acov_def_mat.reshapeCol(k,3,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,3,2);
    // g
    gsMatrix<T,3,1> normal = m_map_def.mine().normal(k).normalized();
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
    m_Gcon_def(2,2) = m_Gcov_def(2,2) = 1.0;
    m_Gcon_def.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2).inverse();

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
    if (basis)
        m_ncov_ori = m_ncov_ori_mat.reshapeCol(k,3,2);

    // Compute full metric
    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori - 2.0 * z * m_Bcov_ori + z*z * m_ncov_ori.transpose()*m_ncov_ori;
    m_Gcon_ori(2,2) = m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori.block(0,0,2,2) = m_Gcov_ori.block(0,0,2,2).inverse();

    if (!basis) return;
    // Compute full basis
    // basis vectors
    m_acov_ori = m_acov_ori_mat.reshapeCol(k,3,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,3,2);
    // g
    gsMatrix<T,3,1> normal = m_map_ori.mine().normal(k).normalized();
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
    m_Gcon_ori(2,2) = m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori.block(0,0,2,2) = m_Gcov_ori.block(0,0,2,2).inverse();

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
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixBaseDim<dim,T>::_evalStretch(const gsMatrix<T> & C, const gsMatrix<T> & gcon_ori) const
{
    gsVector<T> stretches;
    gsMatrix<T> stretchvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    stretches.resize(3,1);    stretches.setZero();
    stretchvec.resize(3,3);   stretchvec.setZero();

    Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;


    GISMO_ENSURE(gcon_ori.cols()!=0,"Is the basis initialized?");

    gsMatrix<T,3,3> B;
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += C(k,l) * gcon_ori.col(k) * gcon_ori.col(l).transpose();

    // gsDebugVar();
    // gsDebugVar((math::almostEqual(((B.block(0,0,2,2).diagonal()).array()-1).sum(),std::numeric_limits<T>::epsilon())));

    gsMatrix<T,3,3> eigVectors;
    gsMatrix<T,3,1> eigValues;
    if (math::abs(((B.block(0,0,2,2).diagonal()).array()-1).sum()) < 10*std::numeric_limits<T>::epsilon())
    {
        eigValues<<0,1,1; // eigenvalues when B = [1 0 0; 0 1 0; 0 0 0]
        eigVectors<<0,0,1,
                    0,1,0,
                    1,0,0;
    }
    else
    {
        eigSolver.compute(B);
        eigValues = eigSolver.eigenvalues();
        eigVectors= eigSolver.eigenvectors();
    }


    stretchvec.leftCols(2) = eigVectors.rightCols(2);
    stretchvec.col(2) = gcon_ori.col(2); // replace with: stretchvec.col(0).template head<3>().cross(stretchvec.col(1).template head<3>())
    stretches.block(0,0,2,1) = eigValues.block(1,0,2,1); // the eigenvalues are a 3x1 matrix, so we need to use matrix block-operations

    // m_stretches.at(2) = 1/m_J0_sq;
    stretches.at(2) = C(2,2);

    for (index_t k=0; k!=3; k++)
        stretches.at(k) = math::sqrt(stretches.at(k));

    result.first = stretches;
    result.second = stretchvec;

    return result;
}


template <short_t dim, class T>
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixBaseDim<dim,T>::_evalPStress(const gsMatrix<T> & S) const
{
    gsVector<T> pstresses;
    gsMatrix<T> pstressvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    pstresses.resize(2,1);    pstresses.setZero();
    pstressvec.resize(3,2);   pstressvec.setZero();

    typename Eigen::SelfAdjointEigenSolver< typename gsMatrix<T>::Base >  eigSolver;

    gsMatrix<T> B(3,3);
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += S(k,l) * m_gcov_ori.col(k) * m_gcov_ori.col(l).transpose();

    eigSolver.compute(B);

    index_t zeroIdx = -1;
    T tol = 1e-14;
    T max = eigSolver.eigenvalues().array().abs().maxCoeff();
    max = (max==0) ? 1 : max;
    for (index_t k=0; k!=3; k++)
        zeroIdx = std::abs(eigSolver.eigenvalues()[k] ) / max < tol ? k : zeroIdx;

    GISMO_ASSERT(zeroIdx!=-1,"No zero found?");

    index_t count = 0;

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

template <short_t dim, class T>
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixBaseDim<dim,T>::_evalPStrain(const gsMatrix<T> & S) const
{
    gsVector<T> pstrains;
    gsMatrix<T> pstrainvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    pstrains.resize(3,1);    pstrains.setZero();
    pstrainvec.resize(3,3);   pstrainvec.setZero();

    typename Eigen::SelfAdjointEigenSolver< typename gsMatrix<T>::Base >  eigSolver;

    gsMatrix<T> B(3,3);
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += S(k,l) * m_gcon_ori.col(k) * m_gcon_ori.col(l).transpose();

    eigSolver.compute(B);

    index_t zeroIdx = -1;
    T tol = 1e-14;
    T max = eigSolver.eigenvalues().array().abs().maxCoeff();
    max = (max==0) ? 1 : max;
    for (index_t k=0; k!=3; k++)
        zeroIdx = std::abs(eigSolver.eigenvalues()[k] ) / max < tol ? k : zeroIdx;

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

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computeStretch(const gsMatrix<T> & C, const gsMatrix<T> & gcon_ori) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalStretch(C,gcon_ori);
    m_stretches = result.first;
    m_stretchvec = result.second;
}

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computePStress(const gsMatrix<T> & S) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalPStress(S);
    m_pstress = result.first;
    m_pstressvec = result.second;
}

template <short_t dim, class T>
void gsMaterialMatrixBaseDim<dim,T>::_computePStrain(const gsMatrix<T> & E) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = _evalPStrain(E);
    m_pstrain = result.first;
    m_pstrainvec = result.second;
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
