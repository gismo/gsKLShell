/** @file gsMaterialMatrixLinear.h

    @brief Provides material matrices for the thin shell class.

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
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>

namespace gismo
{

// Linear material models
template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const std::vector<gsFunction<T>*> &pars
                                        )
                                        :
                                        m_patches(&mp),
                                        m_thickness(&thickness),
                                        m_pars(pars)
{
    initialize();
}

// Linear material models
template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & mp_def,
                                        const gsFunction<T> & thickness,
                                        const std::vector<gsFunction<T>*> &pars
                                        )
                                        :
                                        m_patches(&mp),
                                        m_defpatches(&mp_def),
                                        m_thickness(&thickness),
                                        m_pars(pars)
{
    initialize();
}

// Linear material models
template <short_t dim, class T >
gsMaterialMatrixLinear<dim,T>::gsMaterialMatrixLinear(
                                    const gsFunctionSet<T> & mp,
                                    const gsFunctionSet<T> & mp_def,
                                    const gsFunction<T> & thickness,
                                    const std::vector<gsFunction<T>*> &pars,
                                    const gsFunction<T> & density
                                    )
                                    :
                                    m_patches(&mp),
                                    m_defpatches(&mp_def),
                                    m_thickness(&thickness),
                                    m_pars(pars),
                                    m_density(&density)
{
    initialize();
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
void gsMaterialMatrixLinear<dim,T>::defaultOptions()
{
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T >
short_t gsMaterialMatrixLinear<dim,T>::domainDim() const { return 2; }

template <short_t dim, class T >
short_t gsMaterialMatrixLinear<dim,T>::targetDim() const
{
    if (m_outputType==2)
        return 9;
    else if (m_outputType==1)
        return 3;
    else if (m_outputType==0)
        return 1;
    else if (m_outputType==9)
        return 3;
    else if (m_outputType==10)
        return 2;
    else if (m_outputType==11)
        return 9;
    else
    {
        GISMO_ERROR("This option is unknown");
        return 1;
    }
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::initialize()
{
    // Set default options
    this->defaultOptions();

    // set flags
    m_map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;

    // Initialize some parameters
    m_moment = 0;
    m_outputType = 2;
    m_output = 0; // initialize output type
    m_numPars = m_pars.size();
}


template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::computePoints(const index_t patch, const gsMatrix<T> & u, bool deformed) const
{
    gsMatrix<T> tmp;

    m_map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    this->computeMetricUndeformed();

    if (m_defpatches->nPieces()!=0)
    {
        m_map_def.flags = m_map.flags;
        m_map_def.points = u;
        static_cast<const gsFunction<T>&>(m_defpatches->piece(patch)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
        this->computeMetricDeformed();
    }

    m_thickness->eval_into(m_map.values[0], m_Tmat);

    m_parmat.resize(m_numPars,m_map.values[0].cols());
    m_parmat.setZero();

    for (size_t v=0; v!=m_pars.size(); v++)
    {
        m_pars[v]->eval_into(m_map.values[0], tmp);
        m_parmat.row(v) = tmp;
    }

    m_parvals.resize(m_numPars);
}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.flags = NEED_VALUE;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    result.resize(1, u.cols());
    m_thickness->eval_into(m_map.values[0], m_Tmat);
    m_density->eval_into(m_map.values[0], m_rhomat);
    for (index_t i = 0; i != u.cols(); ++i) // points
    {
        result(0,i) = m_Tmat(0,i)*m_rhomat(0,i);
    }

}

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map.flags = NEED_VALUE;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    m_thickness->eval_into(m_map.values[0], result);
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->computePoints(patch,u);
    gsMatrix<T> result(9, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
                // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
                this->getMetric(k,z(j,k)); // on point i, on height z(0,j)

                gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j*u.cols()+k,3,3);
                /*
                    C = C1111,  C1122,  C1112
                        symm,   C2222,  C2212
                        symm,   symm,   C1212
                */
                C(0,0)          = Cijkl(0,0,0,0); // C1111
                C(1,1)          = Cijkl(1,1,1,1); // C2222
                C(2,2)          = Cijkl(0,1,0,1); // C1212
                C(1,0) = C(0,1) = Cijkl(0,0,1,1); // C1122
                C(2,0) = C(0,2) = Cijkl(0,0,0,1); // C1112
                C(2,1) = C(1,2) = Cijkl(1,1,0,1); // C2212
        }
    }

    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->computePoints(patch,u);
    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
                // this->computeMetric(i,z.at(j),true,true);
                this->getMetric(k,z(j,k)); // on point i, on height z(0,j)

                result(0,j*u.cols()+k) = Sij(0,0);
                result(1,j*u.cols()+k) = Sij(1,1);
                result(2,j*u.cols()+k) = Sij(0,1);
        }
    }

    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u) const
{
    gsMatrix<T> z(1,1);
    z.setZero();
    return eval3D_matrix(patch,u,z);
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixLinear<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u) const
{
    gsMatrix<T> z(1,1);
    z.setZero();
    return eval3D_vector(patch,u,z);
}


/*
    Available class members:
        - m_parvals
        - m_metric
        - m_metric_def
        - m_J0
        - m_J
*/
template <short_t dim, class T >
T gsMaterialMatrixLinear<dim,T>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ENSURE( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);


    // --------------------------
    // Saint Venant Kirchhoff
    // --------------------------
    T lambda, mu, Cconstant;

    mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    GISMO_ENSURE((1.-2.*m_parvals.at(1)) != 0, "Division by zero in construction of SvK material parameters! (1.-2.*m_parvals.at(1)) = "<<(1.-2.*m_parvals.at(1))<<"; m_parvals.at(1) = "<<m_parvals.at(1));
    lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1))) ;
    Cconstant = 2*lambda*mu/(lambda+2*mu);

    return Cconstant*m_Acon_ori(i,j)*m_Acon_ori(k,l) + mu*(m_Acon_ori(i,k)*m_Acon_ori(j,l) + m_Acon_ori(i,l)*m_Acon_ori(j,k));
}

template <short_t dim, class T >
T gsMaterialMatrixLinear<dim,T>::Sij(const index_t i, const index_t j) const
{
    gsMatrix<T> stress;
    // --------------------------
    // Saint Venant Kirchhoff
    // --------------------------
    if (m_moment==0)
    {
        GISMO_ENSURE( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
        stress = 0.5*(m_Acov_def - m_Acov_ori);
    }
    else if (m_moment==2)
    {
        GISMO_ENSURE( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
        stress = (m_Bcov_ori - m_Bcov_def);
    }
    else
    {
        stress.resize(2,2);
        stress.setZero();
        // GISMO_ERROR("Warning: no material model known in simplification, m_moment="<<m_moment);
    }

    // ALTERNATIVE
    // stress = 0.5 * (m_Gcov_def_L - m_Gcov_ori_L);

    T result =  Cijkl(i,j,0,0) * stress(0,0) + Cijkl(i,j,0,1) * stress(0,1)
                + Cijkl(i,j,1,0) * stress(1,0) + Cijkl(i,j,1,1) * stress(1,1);

    return result;
}

// ---------------------------------------------------------------------------------------------------------------------------------
//                                          Metric Computations
// ---------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::computeMetricDeformed() const
{
    computeMetricDeformed_impl<dim>();
}

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixLinear<dim,T>::computeMetricDeformed_impl() const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.points.cols());    m_Acon_def_mat.setZero();
    m_Bcov_def_mat.resize(4,m_map_def.points.cols());    m_Bcov_def_mat.setZero();

    m_acov_def_mat.resize(2*3,m_map_def.points.cols());    m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*3,m_map_def.points.cols());    m_acon_def_mat.setZero();
    m_ncov_def_mat.resize(2*3,m_map_def.points.cols());    m_ncov_def_mat.setZero();

    for (index_t k=0; k!= m_map_def.points.cols(); k++)
    {
        m_acov_def_mat.reshapeCol(k,3,2)   = m_map_def.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,3,2);
        acov = m_map_def.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map_def.deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map_def.normal(k).normalized();

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_Bcov_def_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_def_mat.reshapeCol(k,2,2);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

        gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_def_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);
    }
}

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixLinear<dim,T>::computeMetricDeformed_impl() const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.points.cols());    m_Acon_def_mat.setZero();

    m_acov_def_mat.resize(2*2,m_map_def.points.cols());    m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*2,m_map_def.points.cols());    m_acon_def_mat.setZero();

    for (index_t k=0; k!= m_map_def.points.cols(); k++)
    {
        m_acov_def_mat.reshapeCol(k,2,2)   = m_map_def.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,2,2);
        acov = m_map_def.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,2,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::computeMetricUndeformed() const
{
    computeMetricUndeformed_impl<dim>();
}

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixLinear<dim,T>::computeMetricUndeformed_impl() const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map.points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map.points.cols());    m_Acon_ori_mat.setZero();
    m_Bcov_ori_mat.resize(4,m_map.points.cols());    m_Bcov_ori_mat.setZero();

    m_acov_ori_mat.resize(2*3,m_map.points.cols());    m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*3,m_map.points.cols());    m_acon_ori_mat.setZero();
    m_ncov_ori_mat.resize(2*3,m_map.points.cols());    m_ncov_ori_mat.setZero();

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        m_acov_ori_mat.reshapeCol(k,3,2)   = m_map.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,3,2);
        acov = m_map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_ori_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map.deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map.normal(k).normalized();

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_Bcov_ori_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_ori_mat.reshapeCol(k,2,2);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

        gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_ori_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);
    }
}

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixLinear<dim,T>::computeMetricUndeformed_impl() const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map.points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map.points.cols());    m_Acon_ori_mat.setZero();

    m_acov_ori_mat.resize(2*2,m_map.points.cols());    m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*2,m_map.points.cols());    m_acon_ori_mat.setZero();

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        m_acov_ori_mat.reshapeCol(k,3,2)   = m_map.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,3,2);
        acov = m_map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_ori_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::getMetric(index_t k, T z) const
{
    this->getMetricDeformed(k,z);
    this->getMetricUndeformed(k,z);

    T ratio = m_Gcov_def.determinant() / m_Gcov_ori.determinant();
    GISMO_ENSURE(ratio > 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_Gcov_ori.determinant());
    m_J0_sq = ratio;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::getMetricDeformed(index_t k, T z) const
{
    getMetricDeformed_impl<dim>(k,z);
}

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixLinear<dim,T>::getMetricDeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_ncov_def_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    m_Acon_def = m_Acon_def_mat.reshapeCol(k,2,2);
    m_Bcov_def = m_Bcov_def_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_def = m_acov_def_mat.reshapeCol(k,3,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,3,2);
    m_ncov_def = m_ncov_def_mat.reshapeCol(k,3,2);
    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def - 2.0 * z * m_Bcov_def + z*z * m_ncov_def.transpose()*m_ncov_def;
    m_Gcov_def(2,2) = 1.0;
    m_Gcon_def = m_Gcov_def.inverse();
    // Compute full basis
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

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixLinear<dim,T>::getMetricDeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_def_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    m_Acon_def = m_Acon_def_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_def = m_acov_def_mat.reshapeCol(k,3,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,3,2);
    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def;
    m_Gcov_def(2,2) = 1.0;
    m_Gcon_def = m_Gcov_def.inverse();
    // Compute full basis
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_gcov_def.leftCols(2) = m_acov_def;
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

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::getMetricUndeformed(index_t k, T z) const
{
    getMetricUndeformed_impl<dim>(k,z);
}

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixLinear<dim,T>::getMetricUndeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    m_Acon_ori = m_Acon_ori_mat.reshapeCol(k,2,2);
    m_Bcov_ori = m_Bcov_ori_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_ori = m_acov_ori_mat.reshapeCol(k,3,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,3,2);
    m_ncov_ori = m_ncov_ori_mat.reshapeCol(k,3,2);
    // Compute full metric
    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori - 2.0 * z * m_Bcov_ori + z*z * m_ncov_ori.transpose()*m_ncov_ori;
    m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori = m_Gcov_ori.inverse();
    // Compute full basis
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

template <short_t dim, class T >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixLinear<dim,T>::getMetricUndeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    m_Acon_ori = m_Acon_ori_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_ori = m_acov_ori_mat.reshapeCol(k,3,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,3,2);
    // Compute full metric
    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori;
    m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori = m_Gcov_ori.inverse();
    // Compute full basis
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_gcov_ori.leftCols(2) = m_acov_ori;
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

template <short_t dim, class T >
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrixLinear<dim,T>::evalStretch(const gsMatrix<T> & C) const
{
    gsVector<T> stretches;
    gsMatrix<T> stretchvec;
    std::pair<gsVector<T>,gsMatrix<T>> result;
    stretches.resize(3,1);    stretches.setZero();
    stretchvec.resize(3,3);   stretchvec.setZero();

    Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

    gsMatrix<T> B(3,3);
    B.setZero();
    for (index_t k = 0; k != 2; k++)
        for (index_t l = 0; l != 2; l++)
            B += C(k,l) * m_gcon_ori.col(k) * m_gcon_ori.col(l).transpose();

    eigSolver.compute(B);

    stretchvec.leftCols(2) = eigSolver.eigenvectors().rightCols(2);
    stretchvec.col(2) = m_gcon_ori.col(2);
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

template <short_t dim, class T >
void gsMaterialMatrixLinear<dim,T>::computeStretch(const gsMatrix<T> & C) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = evalStretch(C);
    m_stretches = result.first;
    m_stretchvec = result.second;
}


} // end namespace
