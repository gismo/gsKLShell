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

#include <gsKLShell/gsMaterialMatrix.h>

using namespace gismo;

template <short_t d, typename T>
class Cfun : public gsFunction<T>
{
public:
    Cfun(const gsMaterialMatrixBase<T> * materialMat, index_t patch, const gsVector<T> & u, const T & z)
    :
    m_materialMat(materialMat),
    m_patchID(patch),
    m_points(u),
    m_z(z)
    {

    }

    Cfun(const gsMaterialMatrixBase<T> * materialMat, index_t patch, const gsVector<T> & u)
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
    const gsMaterialMatrixBase<T> * m_materialMat;
    mutable index_t m_patchID;
    mutable gsVector<T> m_points;
    mutable T m_z;
};

namespace gismo
{

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness
                                        )
                                        :
                                        Base(&mp,nullptr,&thickness,nullptr)
{
    _initialize();
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const std::vector<gsFunction<T>*> &pars
                                        )
                                        :
                                        Base(&mp,nullptr,&thickness,nullptr)
{
    GISMO_ASSERT(pars.size()>=2,"Two or more material parameters should be assigned!");
    m_pars = pars;
    _initialize();
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                    const gsFunctionSet<T> & mp,
                                    const gsFunction<T> & thickness,
                                    const std::vector<gsFunction<T>*> &pars,
                                    const gsFunction<T> & Density
                                    )
                                    :
                                    Base(&mp,nullptr,&thickness,&Density)
{
    GISMO_ASSERT(pars.size()>=2,"Two or more material parameters should be assigned!");
    m_pars = pars;
    _initialize();
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::info() const
{
    gsInfo  <<"---------------------------------------------------------------------\n"
            <<"---------------------Hyperelastic Material Info----------------------\n"
            <<"---------------------------------------------------------------------\n\n";

    gsInfo  <<"Material model: \t";
    if (comp)
        gsInfo<<"Compressible ";
    else
        gsInfo<<"Incompressible ";

    if      (mat==Material::SvK)
        gsInfo<<"Saint-Venant Kirchhoff";
    else if (mat==Material::NH)
        gsInfo<<"Neo-Hookean";
    else if (mat==Material::MR)
        gsInfo<<"Mooney-Rivlin";
    else if (mat==Material::OG)
        gsInfo<<"Ogden";
    else if (mat==Material::NH_ext)
        gsInfo<<"Neo-Hookean Extended";
    else
        gsWarn<<"Not specified";
    gsInfo<<"\n";

    gsInfo  <<"Implementation: \t";
    if      (imp==Implementation::Analytical)
        gsInfo<<"Analytical";
    else if (imp==Implementation::Generalized)
        gsInfo<<"Generalized";
    else if (imp==Implementation::Spectral)
        gsInfo<<"Spectral";
    else
        gsWarn<<"Not specified";
    gsInfo<<" implementation\n";

    gsInfo  <<"---------------------------------------------------------------------\n\n";

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_defaultOptions()
{
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_initialize()
{
    // Set default options
    this->_defaultOptions();

    // set flags
    m_map_ori.mine().flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map_ori.mine().flags = NEED_VALUE;
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);

    result.resize(1, u.cols());
    m_thickness->eval_into(m_map_ori.mine().values[0], m_Tmat);
    m_density->eval_into(m_map_ori.mine().values[0], m_rhomat);
    for (index_t i = 0; i != u.cols(); ++i) // points
    {
        result(0,i) = m_Tmat(0,i)*m_rhomat(0,i);
    }

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);

    this->_computePoints(patch,u);

    _stretch_into_impl<comp>(u,result);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(3, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0); // on point i, with height 0.0

        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_J0_sq;

        res = this->_evalStretch(C,m_gcon_ori);
        result.col(i) = res.first;
    }
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(3, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,i);

        this->_getMetric(i,0.0); // on point i, with height 0.0

        // Define objects
        gsMatrix<T,3,3> c, cinv;
        T S33, C3333, dc33;
        // T S33_old;
        S33 = 0.0;
        dc33 = 0.0;
        C3333 = 1.0;

        index_t itmax = 100;
        T tol = 1e-10;

        // Initialize c
        c.setZero();
        c.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        c(2,2) = math::pow(m_J0_sq,-1.0); // c33
        // c(2,2) = 1.0; // c33
        cinv.block(0,0,2,2) = m_Gcon_def.block(0,0,2,2);
        cinv(2,2) = 1.0/c(2,2);

        m_J_sq = m_J0_sq * c(2,2);
        S33 = _Sij(2,2,c,cinv);
        // S33_old = (S33 == 0.0) ? 1.0 : S33;
        C3333   = _Cijkl3D(2,2,2,2,c,cinv);

        dc33 = -2. * S33 / C3333;
        for (index_t it = 0; it < itmax; it++)
        {
            c(2,2) += dc33;

            GISMO_ENSURE(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2) ;

            S33     = _Sij(2,2,c,cinv);
            C3333   = _Cijkl3D(2,2,2,2,c,cinv); //  or _Cijkl???

            dc33 = -2. * S33 / C3333;
            if (abs(dc33) < tol)
            {
                res = this->_evalStretch(c,m_gcon_ori);
                result.col(i) = res.first;
                break;
            }
            GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, abs(dc33) = "<<abs(dc33)<<" and tolerance = "<<tol<<"\n");
        }
    }
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);

    this->_computePoints(patch,u);

    _stretchDir_into_impl<comp>(u,result);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->_getMetric(i,0.0); // on point i, with height 0.0

        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_J0_sq;

        res = this->_evalStretch(C,m_gcon_ori);
        result.col(i) = res.second.reshape(9,1);
    }
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,i);

        this->_getMetric(i,0.0); // on point i, with height 0.0

        // Define objects
        gsMatrix<T,3,3> c, cinv;
        T S33, C3333, dc33;
        // T S33_old;
        S33 = 0.0;
        dc33 = 0.0;
        C3333 = 1.0;

        index_t itmax = 100;
        T tol = 1e-10;

        // Initialize c
        c.setZero();
        c.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        c(2,2) = math::pow(m_J0_sq,-1.0); // c33
        // c(2,2) = 1.0; // c33
        cinv.block(0,0,2,2) = m_Gcon_def.block(0,0,2,2);
        cinv(2,2) = 1.0/c(2,2);

        m_J_sq = m_J0_sq * c(2,2);
        S33 = _Sij(2,2,c,cinv);
        // S33_old = (S33 == 0.0) ? 1.0 : S33;
        C3333   = _Cijkl3D(2,2,2,2,c,cinv);

        dc33 = -2. * S33 / C3333;
        for (index_t it = 0; it < itmax; it++)
        {
            c(2,2) += dc33;

            GISMO_ENSURE(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2) ;

            S33     = _Sij(2,2,c,cinv);
            C3333   = _Cijkl3D(2,2,2,2,c,cinv); //  or _Cijkl???

            dc33 = -2. * S33 / C3333;
            if (abs(dc33) < tol)
            {
                res = this->_evalStretch(c,m_gcon_ori);
                result.col(i) = res.second.reshape(9,1);
                break;
            }
            GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, abs(dc33) = "<<abs(dc33)<<" and tolerance = "<<tol<<"\n");
        }
    }
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map_ori.mine().flags = NEED_VALUE;
    m_map_ori.mine().points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map_ori);
    m_thickness->eval_into(m_map_ori.mine().values[0], result);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::covtransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::contransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput out) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    return _eval3D_matrix_C_impl<mat,comp>(Cmat,patch,u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_matrix_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result = _eval3D_Incompressible_matrix_C(Cmat,patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_matrix_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result = _eval3D_Incompressible_matrix_C(Cmat,patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_matrix_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result = _eval3D_Compressible_matrix_C(Cmat,patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_dmatrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return _eval3D_dmatrix_impl<mat,comp>(patch,u,z,out);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!(!_comp && _mat==Material::NH), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_dmatrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return eval3D_dmatrix_num(patch,u,z,out);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if< (!_comp && _mat==Material::NH), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_dmatrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return eval3D_dmatrix_num(patch,u,z,out);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_dmatrix_num(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // return gsMatrix<T>::Zero(27,u.cols()*z.rows());

    this->_computePoints(patch,u);

    gsMatrix<T> result(27,u.cols()*z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {

        // // Evaluate material properties on the quadrature point
        // for (index_t v=0; v!=m_parmat.rows(); v++)
        //     m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_Tmat(0, k), false); // on point i, on height z(0,j)
            result.col(j*u.cols()+k) = dCijkl(patch,u.col(k),z(j,k));
        }
    }

    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dCijkl(const index_t patch, const gsVector<T> & u, const T z) const
{
    return dCijkl_impl<mat,comp>(patch,u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!(!_comp && _mat==Material::NH), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dCijkl_impl(const index_t patch, const gsVector<T> & u, const T z) const
{
    gsMatrix<T> result;
    Cfun<dim,T> Cfunc(this,patch,u,z);
    gsMatrix<T> Cvoight(3,1);

    Cvoight(0,0) = m_Gcov_def(0,0);
    Cvoight(1,0) = m_Gcov_def(1,1);
    Cvoight(2,0) = m_Gcov_def(0,1);

    Cfunc.deriv_into(Cvoight,result);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<(!_comp && _mat==Material::NH), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dCijkl_impl(const index_t patch, const gsVector<T> & u, const T z) const
{
    gsMatrix<T> dCdC(3,9);
    dCdC(0,0) = dCijkl_dCmn(0,0,0,0,  0,0); // C1111dC11
    dCdC(1,0) = dCijkl_dCmn(0,0,0,0,  1,1); // C1111dC22
    dCdC(2,0) = dCijkl_dCmn(0,0,0,0,  0,1); // C1111dC12

    dCdC(0,1) = dCijkl_dCmn(0,0,1,1,  0,0); // C2211dC11
    dCdC(1,1) = dCijkl_dCmn(0,0,1,1,  1,1); // C2211dC22
    dCdC(2,1) = dCijkl_dCmn(0,0,1,1,  0,1); // C2211dC12

    dCdC(0,2) = dCijkl_dCmn(0,0,0,1,  0,0); // C1211dC11
    dCdC(1,2) = dCijkl_dCmn(0,0,0,1,  1,1); // C1211dC22
    dCdC(2,2) = dCijkl_dCmn(0,0,0,1,  0,1); // C1211dC12

    dCdC(0,3) = dCijkl_dCmn(1,1,0,0,  0,0); // C1122dC11
    dCdC(1,3) = dCijkl_dCmn(1,1,0,0,  1,1); // C1122dC22
    dCdC(2,3) = dCijkl_dCmn(1,1,0,0,  0,1); // C1122dC12

    dCdC(0,4) = dCijkl_dCmn(1,1,1,1,  0,0); // C2222dC11
    dCdC(1,4) = dCijkl_dCmn(1,1,1,1,  1,1); // C2222dC22
    dCdC(2,4) = dCijkl_dCmn(1,1,1,1,  0,1); // C2222dC12

    dCdC(0,5) = dCijkl_dCmn(1,1,0,1,  0,0); // C1222dC11
    dCdC(1,5) = dCijkl_dCmn(1,1,0,1,  1,1); // C1222dC22
    dCdC(2,5) = dCijkl_dCmn(1,1,0,1,  0,1); // C1222dC12

    dCdC(0,6) = dCijkl_dCmn(0,1,0,0,  0,0); // C1112dC11
    dCdC(1,6) = dCijkl_dCmn(0,1,0,0,  1,1); // C1112dC22
    dCdC(2,6) = dCijkl_dCmn(0,1,0,0,  0,1); // C1112dC12

    dCdC(0,7) = dCijkl_dCmn(0,1,1,1,  0,0); // C2212dC11
    dCdC(1,7) = dCijkl_dCmn(0,1,1,1,  1,1); // C2212dC22
    dCdC(2,7) = dCijkl_dCmn(0,1,1,1,  0,1); // C2212dC12

    dCdC(0,8) = dCijkl_dCmn(0,1,0,1,  0,0); // C1212dC11
    dCdC(1,8) = dCijkl_dCmn(0,1,0,1,  1,1); // C1212dC22
    dCdC(2,8) = dCijkl_dCmn(0,1,0,1,  0,1); // C1212dC12
    return dCdC.reshape(27,1);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dCijkl_dCmn(const index_t i, const index_t j, const index_t k, const index_t l, const index_t m, const index_t n) const
{
    return dCijkl_dCmn_impl<mat,comp>(i,j,k,l,m,n);
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<(!_comp && _mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dCijkl_dCmn_impl(const index_t i, const index_t j, const index_t k, const index_t l, const index_t m, const index_t n) const
{
    // --------------------------
    // Neo-Hookean
    // --------------------------
    gsMatrix<T> Cinv = m_Gcon_def;
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    return
            -mu*1./m_J0_sq*Cinv(m,n)*(2.*Cinv(i,j)*Cinv(k,l) + Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k))
            -mu   *m_J0_sq*1./2. *(
                                        (Cinv(j,k)*Cinv(l,n) + Cinv(j,l)*Cinv(k,n) + 2*Cinv(j,n)*Cinv(k,l))*Cinv(i,m)
                                       +(Cinv(j,k)*Cinv(l,m) + Cinv(j,l)*Cinv(k,m) + 2*Cinv(j,m)*Cinv(k,l))*Cinv(i,n)
                                       +(Cinv(i,k)*Cinv(l,n) + Cinv(i,l)*Cinv(k,n))*Cinv(j,m)
                                       +(Cinv(i,k)*Cinv(l,m) + Cinv(i,l)*Cinv(k,m))*Cinv(j,n)
                                     +2*(Cinv(k,m)*Cinv(l,n) + Cinv(k,n)*Cinv(l,m))*Cinv(i,j)
                                    );

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!(!_comp && _mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dCijkl_dCmn_impl(const index_t i, const index_t j, const index_t k, const index_t l, const index_t m, const index_t n) const
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_dmatrix_an(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // return gsMatrix<T>::Zero(27,u.cols()*z.rows());

    this->_computePoints(patch,u);

    gsMatrix<T> result(27,u.cols()*z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {

        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_Tmat(0, k), false); // on point i, on height z(0,j)

            gsAsMatrix<T,Dynamic,Dynamic> dCdC = result.reshapeCol(k,3,9);
            dCdC(0,0) = dCijkl_dCmn(0,0,0,0,  0,0); // C1111dC11
            dCdC(1,0) = dCijkl_dCmn(0,0,0,0,  1,1); // C1111dC22
            dCdC(2,0) = dCijkl_dCmn(0,0,0,0,  0,1); // C1111dC12

            dCdC(0,1) = dCijkl_dCmn(0,0,1,1,  0,0); // C2211dC11
            dCdC(1,1) = dCijkl_dCmn(0,0,1,1,  1,1); // C2211dC22
            dCdC(2,1) = dCijkl_dCmn(0,0,1,1,  0,1); // C2211dC12

            dCdC(0,2) = dCijkl_dCmn(0,0,0,1,  0,0); // C1211dC11
            dCdC(1,2) = dCijkl_dCmn(0,0,0,1,  1,1); // C1211dC22
            dCdC(2,2) = dCijkl_dCmn(0,0,0,1,  0,1); // C1211dC12

            dCdC(0,3) = dCijkl_dCmn(1,1,0,0,  0,0); // C1122dC11
            dCdC(1,3) = dCijkl_dCmn(1,1,0,0,  1,1); // C1122dC22
            dCdC(2,3) = dCijkl_dCmn(1,1,0,0,  0,1); // C1122dC12

            dCdC(0,4) = dCijkl_dCmn(1,1,1,1,  0,0); // C2222dC11
            dCdC(1,4) = dCijkl_dCmn(1,1,1,1,  1,1); // C2222dC22
            dCdC(2,4) = dCijkl_dCmn(1,1,1,1,  0,1); // C2222dC12

            dCdC(0,5) = dCijkl_dCmn(1,1,0,1,  0,0); // C1222dC11
            dCdC(1,5) = dCijkl_dCmn(1,1,0,1,  1,1); // C1222dC22
            dCdC(2,5) = dCijkl_dCmn(1,1,0,1,  0,1); // C1222dC12

            dCdC(0,6) = dCijkl_dCmn(0,1,0,0,  0,0); // C1112dC11
            dCdC(1,6) = dCijkl_dCmn(0,1,0,0,  1,1); // C1112dC22
            dCdC(2,6) = dCijkl_dCmn(0,1,0,0,  0,1); // C1112dC12

            dCdC(0,7) = dCijkl_dCmn(0,1,1,1,  0,0); // C2212dC11
            dCdC(1,7) = dCijkl_dCmn(0,1,1,1,  1,1); // C2212dC22
            dCdC(2,7) = dCijkl_dCmn(0,1,1,1,  0,1); // C2212dC12

            dCdC(0,8) = dCijkl_dCmn(0,1,0,1,  0,0); // C1212dC11
            dCdC(1,8) = dCijkl_dCmn(0,1,0,1,  1,1); // C1212dC22
            dCdC(2,8) = dCijkl_dCmn(0,1,0,1,  0,1); // C1212dC12

            dCdC.row(2) *= 2;
        }
    }

    return result;
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    return _eval3D_matrix_impl<mat,comp>(patch,u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result = _eval3D_Incompressible_matrix(patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result = _eval3D_Incompressible_matrix(patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    this->_computePoints(patch,u);
    gsMatrix<T> result = _eval3D_Compressible_matrix(patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return this->eval3D_stress(patch,u,z,out);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{

    this->_computePoints(patch,u);

    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();
    gsMatrix<T,3,3> E;
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_Tmat(0, k)); // on point i, on height z(0,j)

            E = 0.5 * (m_Gcov_def - m_Gcon_ori);
            res = this->_evalPStrain(E);
            result.col(j * u.cols() + k) = res.first;
        }
    }

    return result;

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
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
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    this->_computePoints(patch,u);

    gsMatrix<T> result(3, u.cols() * z.rows());
    gsMatrix<T,2,2> E;
    result.setZero();
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_Tmat(0, k)); // on point i, on height z(0,j)

            E = 0.5 * (m_Gcov_def.block(0,0,2,2) - m_Gcov_ori.block(0,0,2,2));
            result(0,j*u.cols() + k) = E(0,0);
            result(1,j*u.cols() + k) = E(1,1);
            result(2,j*u.cols() + k) = E(0,1) + E(1,0);
        }
    }

    return result;

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput ) const
{
    this->_computePoints(patch,u);
    return this->_eval3D_stress_impl<mat,comp>(patch,u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_stress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = _eval3D_Incompressible_stress(patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_stress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = _eval3D_Incompressible_stress(patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_stress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = _eval3D_Compressible_stress(patch, u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_Compressible_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    gsMatrix<T> result(2, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
            this->_getMetric(k, z(j, k) * m_Tmat(0, k)); // on point i, on height z(0,j)

            // Define objects
            gsMatrix<T,3,3> c, cinv;
            T S33, C3333, dc33;
            // T S33_old;
            S33 = 0.0;
            dc33 = 0.0;
            C3333 = 1.0;

            index_t itmax = 100;
            T tol = 1e-10;

            // Initialize c
            c.setZero();
            c.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
            c(2,2) = math::pow(m_J0_sq,-1.0); // c33
            // c(2,2) = 1.0; // c33

            cinv.setZero();
            cinv.block(0,0,2,2) = m_Gcon_def.block(0,0,2,2);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2);
            S33 = _Sij(2,2,c,cinv);
            // S33_old = (S33 == 0.0) ? 1.0 : S33;
            C3333   = _Cijkl3D(2,2,2,2,c,cinv);

            dc33 = -2. * S33 / C3333;
            for (index_t it = 0; it < itmax; it++)
            {
                c(2,2) += dc33;

                GISMO_ENSURE(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
                cinv(2,2) = 1.0/c(2,2);

                m_J_sq = m_J0_sq * c(2,2) ;

                S33     = _Sij(2,2,c,cinv);
                C3333   = _Cijkl3D(2,2,2,2,c,cinv); //  or _Cijkl???

                dc33 = -2. * S33 / C3333;
                // if (abs(S33/S33_old) < tol)
                if (abs(dc33) < tol)
                {

                    result(0,j*u.cols()+k) = _Sij(0,0,c,cinv); // S11
                    result(1,j*u.cols()+k) = _Sij(1,1,c,cinv); // S22
                    result(2,j*u.cols()+k) = _Sij(0,1,c,cinv); // S12

                    break;
                }
                GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n");
            }
        }
    }
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_Incompressible_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    gsMatrix<T> result(3, u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            this->_getMetric(k, z(j, k) * m_Tmat(0, k)); // on point i, on height z(0,j)

            result(0, j * u.cols() + k) = _Sij(0, 0);
            result(1, j * u.cols() + k) = _Sij(1, 1);
            result(2, j * u.cols() + k) = _Sij(0, 1);
        }
    }
    return result;
}

//////// Material getters and setters
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::setYoungsModulus(const gsFunction<T> & YoungsModulus)
{
    if ((index_t)m_pars.size() < 1)
        m_pars.resize(1);
    m_pars[0] = const_cast<gsFunction<T> *>(&YoungsModulus);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsFunction<T> * gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getYoungsModulus()
{
    GISMO_ASSERT((index_t)m_pars.size()>=1,"Something went wrong: Size of parameter container should be >= 1, but is "<<m_pars.size());
    return m_pars[0];
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::setPoissonsRatio(const gsFunction<T> & PoissonsRatio)
{
    if ((index_t)m_pars.size() < 2)
        m_pars.resize(2);
    m_pars[1] = const_cast<gsFunction<T> *>(&PoissonsRatio);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsFunction<T> * gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getPoissonsRatio()
{
    GISMO_ASSERT((index_t)m_pars.size()>=2,"Something went wrong: Size of parameter container should be >= 2, but is "<<m_pars.size());
    return m_pars[1];
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_setRatio_impl(const gsFunction<T> & Ratio)
{
    if ((index_t)m_pars.size() < 3)
        m_pars.resize(3);
    m_pars[2] = const_cast<gsFunction<T> *>(&Ratio);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat!=Material::MR, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_setRatio_impl(const gsFunction<T> & Ratio)
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, gsFunction<T> * >::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_getRatio_impl()
{
    GISMO_ASSERT((index_t)m_pars.size()>2,"Something went wrong: Size of parameter container should be > 2, but is "<<m_pars.size());
    return m_pars[2];
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat!=Material::MR, gsFunction<T> * >::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_getRatio_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::OG, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_setMu_impl(const index_t & i, const gsFunction<T> & Mu_i)
{
    if ((index_t)m_pars.size() < 2+2*i+1)
        m_pars.resize(2+2*i+1);
    m_pars[2+2*i] = const_cast<gsFunction<T> *>(&Mu_i);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat!=Material::OG, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_setMu_impl(const index_t & i, const gsFunction<T> & Mu_i)
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::OG, gsFunction<T> * >::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_getMu_impl(const index_t & i)
{
    GISMO_ASSERT((index_t)m_pars.size()>2+2*i,"Something went wrong: Size of parameter container should be > 2, but is "<<m_pars.size());
    return m_pars[2+2*i];
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat!=Material::OG, gsFunction<T> * >::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_getMu_impl(const index_t & i)
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::OG, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_setAlpha_impl(const index_t & i, const gsFunction<T> & Alpha_i)
{
    if ((index_t)m_pars.size() < 2+2*i+2)
        m_pars.resize(2+2*i+2);
    m_pars[2+2*i+1] = const_cast<gsFunction<T> *>(&Alpha_i);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat!=Material::OG, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_setAlpha_impl(const index_t & i, const gsFunction<T> & Alpha_i)
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::OG, gsFunction<T> * >::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_getAlpha_impl(const index_t & i)
{
    GISMO_ASSERT((index_t)m_pars.size()>2+2*i+1,"Something went wrong: Size of parameter container should be > 2, but is "<<m_pars.size());
    return m_pars[2+2*i+1];
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat!=Material::OG, gsFunction<T> * >::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_getAlpha_impl(const index_t & i)
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    gsMatrix<T> Spmat = eval3D_pstress(patch,u,z);
    gsMatrix<T> Epmat = eval3D_pstrain(patch,u,z);
    gsVector<T> Sp, Ep;
    gsMatrix<T> result(1, u.cols() * z.rows());
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            // See Nakashino 2020
            // Smin = Sp[0], Smax = Sp[1], S33 = Sp[2]
            // Emin = Ep[0], Emax = Ep[1], E33 = Ep[2]
            Sp = Spmat.col(j*u.cols() + k);
            Ep = Epmat.col(j*u.cols() + k);
            if (Sp[0] >= 0) // taut
                result.col(j * u.cols() + k) << 1;
            else if (Sp[1] < 0) // slack
            // else if (Ep[1] < 0) // slack
                result.col(j * u.cols() + k) << -1;
            else // wrinkled
                result.col(j * u.cols() + k) << 0;
        }
    }
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_Incompressible_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    gsMatrix<T> result(9, u.cols());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        this->_getMetric(0, z * m_Tmat(0, 0), Cmat, false); // on point i, on height z(0,j)

        gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(k,3,3);
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

    return result;
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_Incompressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

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
                this->_getMetric(k,z(j,k) * m_Tmat(0,k) ); // on point i, on height z(0,j)

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

/*
    Available class members:
        - m_parvals
        - m_metric
        - m_metric_def
        - m_J0
        - m_J
*/
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ENSURE( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(!comp,"Material model is not incompressible?");

    return _Cijkl_impl<mat,imp>(i,j,k,l);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    // --------------------------
    // Neo-Hookean
    // --------------------------
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    return mu*1./m_J0_sq*(2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k));
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    // --------------------------
    // Mooney-Rivlin
    // Parameter 3 is the ratio between c1 and c2.; c1 = m_parvals.at(2)*c2
    // --------------------------
    GISMO_ENSURE(m_pars.size()==3,"Mooney-Rivlin model needs to be a 3 parameter model");
    T traceCt =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);

    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    T c2 = mu/(m_parvals.at(2) + 1);
    T c1 = m_parvals.at(2)*c2;

    T Gabcd = - 1./2. * ( m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k) );

    return (c1 + c2 * traceCt) *1./m_J0_sq*(2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k))// correct
            - 2. * c2 / m_J0_sq * ( m_Gcon_ori(i,j) * m_Gcon_def(k,l) + m_Gcon_def(i,j)*m_Gcon_ori(k,l)) // correct
            + 2. * c2 * m_J0_sq * ( Gabcd + m_Gcon_def(i,j)*m_Gcon_def(k,l) ); // Roohbakhshan
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp==Implementation::Spectral, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    // --------------------------
    // Stretch-based implementations
    // --------------------------
    T tmp = 0.0;
    gsMatrix<T> C(3,3);
    C.setZero();
    C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
    // C.block(0,0,2,2) = (m_gcov_def.transpose()*m_gcov_def).block(0,0,2,2);
    // gsDebugVar(m_gcov_def.transpose()*m_gcov_def);
    C(2,2) = 1./m_J0_sq;

    this->_computeStretch(C,m_gcon_ori);

    tmp = 0.0;
    for (index_t a = 0; a != 2; a++)
    {
        // C_iiii
        tmp +=  _Cabcd(a,a,a,a)*(
                ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                );

        for (index_t b = a+1; b != 2; b++)
        {
            // C_iijj = C_jjii
            tmp +=  _Cabcd(a,a,b,b)*(
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                        +
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                    );

            // C_ijij = Cjiji = Cijji = Cjiij
            tmp +=  _Cabcd(a,b,a,b)*(
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                        +
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                        +
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                        +
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                    );
        }
    }
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp==Implementation::Generalized, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    // --------------------------
    // General implementations
    // --------------------------
    return 4.0 * _d2Psi(i,j,k,l) + 4.0 * _d2Psi(2,2,2,2)*math::pow(m_J0_sq,-2.0)*m_Gcon_def(i,j)*m_Gcon_def(k,l)
            - 4.0/ m_J0_sq  * ( _d2Psi(2,2,i,j)*m_Gcon_def(k,l) + _d2Psi(2,2,k,l)*m_Gcon_def(i,j) )
            + 2.0 * _dPsi(2,2) / m_J0_sq * (2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k));
}

// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// template <enum Material _mat>
// // OTHER CASES!
// typename std::enable_if<(_mat >= 30) && !(_mat==0) && !(_mat==2) && !(_mat==3), T>::type
// gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
// {
//     GISMO_ERROR("Material model unknown (model = "<<_mat<<"). Use gsMaterialMatrix<dim,T,matId,comp,mat,imp>::info() to see the options.");
// }

        // else if (m_material==4)
        //         GISMO_ERROR("Material model 4 is not invariant-based! Use 14 instead...");
        // else if (m_material==5)
        //         GISMO_ERROR("Material model 5 is only for compressible materials...");


// Condensation of the 3D tensor for compressible materials
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE(c.cols()==c.rows(),"Matrix c must be square");
    GISMO_ENSURE(c.cols()==3,"Matrix c must be 3x3");
    GISMO_ENSURE(cinv.cols()==cinv.rows(),"Matrix cinv must be square");
    GISMO_ENSURE(cinv.cols()==3,"Matrix cinv must be 3x3");
    GISMO_ENSURE( ( (i <2) && (j <2) && (k <2) && (l <2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(comp,"Material model is not compressible?");

    return _Cijkl_impl<imp>(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Implementation _imp>
typename std::enable_if< _imp==Implementation::Spectral, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // static condensation is done before the projection
    this->_computeStretch(c,m_gcon_ori);
    return _Cijkl3D(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Implementation _imp>
typename std::enable_if<!(_imp==Implementation::Spectral), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    return _Cijkl3D(i,j,k,l,c,cinv) - ( _Cijkl3D(i,j,2,2,c,cinv) * _Cijkl3D(2,2,k,l,c,cinv) ) / _Cijkl3D(2,2,2,2,c,cinv);
}

// 3D tensor for compressible materials
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE( ( (i <3) && (j <3) && (k <3) && (l <3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    return _Cijkl3D_impl<mat,imp>(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::SvK && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ERROR("Compressible material matrix requested, but not needed. How?");
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Neo-Hookean
    // --------------------------

    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");

    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    return 1.0 / 9.0 * mu * math::pow( m_J_sq , -1.0/3.0 ) * ( 2.0 * I_1 * ( cinv(i,j)*cinv(k,l) - 3.0 * dCinv )
                        - 6.0 *( m_Gcon_ori(i,j)*cinv(k,l) + cinv(i,j)*m_Gcon_ori(k,l) ) )
            + K * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1)*dCinv );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::MR && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Mooney-Rivlin
    // --------------------------
    GISMO_ENSURE(m_pars.size()==3,"Mooney-Rivlin model needs to be a 3 parameter model");

    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");

    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    T I_2 = c(2,2) * traceCt + m_J0_sq;
    T d2I_2 = idelta(i,2)*idelta(j,2)*idelta(k,2)*idelta(l,2)*( m_J0_sq*( cinv(i,j)*cinv(k,l) + dCinv ) )
            + delta(i,2)*delta(j,2)*idelta(k,2)*idelta(l,2)*_dI_1(k,l)
            + idelta(i,2)*idelta(j,2)*delta(k,2)*delta(l,2)*_dI_1(i,j);
    T c2 = mu/(m_parvals.at(2) + 1);
    T c1 = m_parvals.at(2)*c2;

    return  1.0/9.0 * c1 * math::pow(m_J_sq, -1.0/3.0) *  ( 2.0*I_1*cinv(i,j)*cinv(k,l) - 6.0*I_1*dCinv
                                                            - 6.0*_dI_1(i,j)*cinv(k,l)     - 6.0*cinv(i,j)*_dI_1(k,l) ) // + 9*d2I_1 = 0
            + 1.0/9.0 * c2 * math::pow(m_J_sq, -2.0/3.0) *  ( 8.0*I_2*cinv(i,j)*cinv(k,l) - 12.0*I_2*dCinv
                                                                - 12.0*_dI_2(i,j,c,cinv)*cinv(k,l)- 12.0*cinv(i,j)*_dI_2(k,l,c,cinv)
                                                                + 18.0*d2I_2 )
            + K * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1)*dCinv );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH_ext && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Neo-Hookean 2
    // --------------------------
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");

    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    return - 2.0 * mu * dCinv + lambda * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1)*dCinv );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Spectral, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Stretch-based implementations
    // --------------------------
    T tmp;
    if ( (i==2) && (j==2) && (k==2) && (l==2) ) // if C3333 (no static condensation)
        tmp = _Cabcd(2,2,2,2);
    else
    {
        tmp = 0.0;
        T C = 0.0;
        T C2222 = _Cabcd(2,2,2,2);
        // T Cab22,C22ab;
        for (index_t a = 0; a != 2; a++)
        {
            // C_iiii
            // if (!((i==2 || j==2 || k==2 || l==2) && a!=2))
            // C = _Cabcd(a,a,a,a);
            C = _Cabcd(a,a,a,a) - math::pow(_Cabcd(2,2,a,a),2) / C2222;
            tmp +=  C*(
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                    );

            for (index_t b = a+1; b != 2; b++)
            {
                // C_iijj
                // C = _Cabcd(a,a,b,b);
                C = _Cabcd(a,a,b,b) - _Cabcd(a,a,2,2) * _Cabcd(2,2,b,b) / C2222;
                tmp +=  C*(
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                            +
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                        );

                // C_ijij = Cjiji = Cijji = Cjiij
                // C = _Cabcd(a,b,a,b);
                C = _Cabcd(a,b,a,b) - math::pow(_Cabcd(2,2,a,b),2) / C2222;
                tmp +=  C*(
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                            +
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                            +
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                            +
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                        );
            }
        }
    }

    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Generalized, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // General implementations
    // --------------------------
    return 4.0 * _d2Psi(i,j,k,l,c,cinv);
}

// Plan of approach:
// - Make function getMetricDeformed(k,z,basis,E), where k,z and basis are used to compute the undeformed metric
// - this new getMetricDeformed will compute the metric based on E, using m_Gcov_def etc =2E+m_Gcov_ori etc
// - Make an eval3D_matrix(u,z,...,E) which computes using E
// - Call that eval matrix inside the TFT given the u,z where the strain is also evaluated

// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_C(const gsMatrix<T> & Cmat) const
// {
//     gsMatrix<T> result(9, u.cols() * z.rows());
//     result.setZero();

//     for (index_t k=0; k!=u.cols(); k++)
//     {
//         // Evaluate material properties on the quadrature point
//         for (index_t v=0; v!=m_parmat.rows(); v++)
//             m_parvals.at(v) = m_parmat(v,k);

//         for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
//         {
//             this->_getMetric(k, z(j,k) * m_Tmat(0,k), Cmat);

//         }
//     }

// }

/*
    Available class members:
        - m_parvals.at(0)
        - m_parvals.at(1)
        - m_metric
        - m_metric_def
        - m_J0
        - m_J
        - m_Cinv
*/
// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij(const index_t i, const index_t j) const { _Sij(i,j,NULL,NULL); }

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij(const index_t i, const index_t j) const
{
    return _Sij_impl<mat,imp>(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::SvK && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j) const
{
    GISMO_ERROR("Not implemented");
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j) const
{
    // --------------------------
    // Neo-Hoookean
    // --------------------------
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    return mu * (m_Gcon_ori(i,j) - 1./m_J0_sq * m_Gcon_def(i,j) );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::MR && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j) const
{
    // --------------------------
    // Mooney-Rivlin
    // Parameter 3 is the ratio between c1 and c2.
    // --------------------------
    GISMO_ENSURE(m_pars.size()==3,"Mooney-Rivlin model needs to be a 3 parameter model");
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    T traceCt =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);

    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;

    return  c1 * ( m_Gcon_ori(i,j) - 1/m_J0_sq * m_Gcon_def(i,j) )
            + c2 / m_J0_sq * (m_Gcon_ori(i,j) - traceCt * m_Gcon_def(i,j) ) + c2 * m_J0_sq * m_Gcon_def(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Spectral, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j) const
{
    T tmp = 0.0;
    gsMatrix<T> C(3,3);
    C.setZero();
    C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
    C(2,2) = 1./m_J0_sq;

    this->_computeStretch(C,m_gcon_ori);

    for (index_t a = 0; a != 2; a++)
    {
        tmp += _Sa(a)*(
                    ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )
                    );
    }
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Generalized, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j) const
{
    // --------------------------
    // Generalized
    // --------------------------
    return 2.0 * _dPsi(i,j) - 2.0 * _dPsi(2,2) * math::pow(m_J0_sq,-1.0)*m_Gcon_def(i,j);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    return _Sij_impl<mat,imp>(i,j,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Neo-Hoookean
    // --------------------------
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);

    return  mu * math::pow( m_J_sq , -1.0/3.0 ) * ( m_Gcon_ori(i,j) - 1.0/3.0 * I_1 * cinv(i,j) )
            + K * 0.5 * ( m_J_sq - 1.0 ) * cinv(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::MR && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Mooney-Rivlin
    // Parameter 3 is the ratio between c1 and c2.
    // --------------------------
    GISMO_ENSURE(m_pars.size()==3,"Mooney-Rivlin model needs to be a 3 parameter model");
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    T I_2   = c(2,2) * traceCt + m_J0_sq;

    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;

    return  c1 * math::pow( m_J_sq , -1.0/3.0 ) * ( m_Gcon_ori(i,j) - 1.0/3.0 * I_1 * cinv(i,j) )
            + c2 * math::pow( m_J_sq , -2.0/3.0 ) * ( _dI_2(i,j,c,cinv)- 2.0/3.0 * I_2 * cinv(i,j) )
            + K * 0.5 * ( m_J_sq - 1.0 ) * cinv(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH_ext && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Neo-Hookean 2
    // --------------------------
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
    return mu * m_Gcon_ori(i,j) - mu * cinv(i,j) + lambda / 2.0 * ( m_J_sq - 1 ) * cinv(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Spectral, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T tmp = 0.0;
    this->_computeStretch(c,m_gcon_ori);
    for (index_t a = 0; a != 3; a++)
    {
        tmp += _Sa(a)*(
                    ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )
                    );
    }
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Generalized, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Generalized
    // --------------------------
    return 2.0 * _dPsi(i,j,c,cinv);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sii(const index_t i) const // principle stresses
{
    gsMatrix<T> C(3,3);
    C.setZero();
    C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
    C(2,2) = 1./m_J0_sq;
    return _Sii(i,C);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sii(const index_t i, const gsMatrix<T> & c) const
{
    this->_computeStretch(c,m_gcon_ori);
    return _Sa(i);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_Compressible_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    gsMatrix<T> result(9, u.cols());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        this->_getMetric(0, z * m_Tmat(0, 0), Cmat, false); // on point i, on height z(0,j)


        // Define objects
        gsMatrix<T,3,3> c, cinv;
        T S33, C3333, dc33;
        // T S33_old;
        S33 = 0.0;
        dc33 = 0.0;
        C3333 = 1.0;

        index_t itmax = 100;
        T tol = 1e-10;

        // Initialize c
        c.setZero();
        c.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        c(2,2) = math::pow(m_J0_sq,-1.0); // c33
        // c(2,2) = 1.0; // c33
        cinv.setZero();
        cinv.block(0,0,2,2) = m_Gcon_def.block(0,0,2,2);
        cinv(2,2) = 1.0/c(2,2);

        m_J_sq = m_J0_sq * c(2,2);
        S33 = _Sij(2,2,c,cinv);
        // S33_old = (S33 == 0.0) ? 1.0 : S33;
        C3333   = _Cijkl3D(2,2,2,2,c,cinv);

        dc33 = -2. * S33 / C3333;
        for (index_t it = 0; it < itmax; it++)
        {
            c(2,2) += dc33;

            GISMO_ENSURE(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2) ;

            S33     = _Sij(2,2,c,cinv);
            C3333   = _Cijkl3D(2,2,2,2,c,cinv); //  or _Cijkl???

            dc33 = -2. * S33 / C3333;
            // if (abs(S33/S33_old) < tol)
            if (abs(dc33) < tol)
            {
                gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(k,3,3);
                /*
                    C = C1111,  C1122,  C1112
                        symm,   C2222,  C2212
                        symm,   symm,   C1212
                */
                C(0,0)          = _Cijkl(0,0,0,0,c,cinv); // C1111
                C(1,1)          = _Cijkl(1,1,1,1,c,cinv); // C2222
                C(2,2)          = _Cijkl(0,1,0,1,c,cinv); // C1212
                C(1,0) = C(0,1) = _Cijkl(0,0,1,1,c,cinv); // C1122
                C(2,0) = C(0,2) = _Cijkl(0,0,0,1,c,cinv); // C1112
                C(2,1) = C(1,2) = _Cijkl(1,1,0,1,c,cinv); // C2212

                break;
            }
            GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n");
        }
    }
    return result;
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_eval3D_Compressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
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
            this->_getMetric(k, z(j, k) * m_Tmat(0, k)); // on point i, on height z(0,j)

            // Define objects
            gsMatrix<T,3,3> c, cinv;
            T S33, C3333, dc33;
            // T S33_old;
            S33 = 0.0;
            dc33 = 0.0;
            C3333 = 1.0;

            index_t itmax = 100;
            T tol = 1e-10;

            // Initialize c
            c.setZero();
            c.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
            c(2,2) = math::pow(m_J0_sq,-1.0); // c33
            // c(2,2) = 1.0; // c33
            cinv.setZero();
            cinv.block(0,0,2,2) = m_Gcon_def.block(0,0,2,2);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2);
            S33 = _Sij(2,2,c,cinv);
            // S33_old = (S33 == 0.0) ? 1.0 : S33;
            C3333   = _Cijkl3D(2,2,2,2,c,cinv);

            dc33 = -2. * S33 / C3333;
            for (index_t it = 0; it < itmax; it++)
            {
                c(2,2) += dc33;

                GISMO_ENSURE(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
                cinv(2,2) = 1.0/c(2,2);

                m_J_sq = m_J0_sq * c(2,2) ;

                S33     = _Sij(2,2,c,cinv);
                C3333   = _Cijkl3D(2,2,2,2,c,cinv); //  or _Cijkl???

                dc33 = -2. * S33 / C3333;
                // if (abs(S33/S33_old) < tol)
                if (abs(dc33) < tol)
                {
                    gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j*u.cols()+k,3,3);
                    /*
                        C = C1111,  C1122,  C1112
                            symm,   C2222,  C2212
                            symm,   symm,   C1212
                    */
                    C(0,0)          = _Cijkl(0,0,0,0,c,cinv); // C1111
                    C(1,1)          = _Cijkl(1,1,1,1,c,cinv); // C2222
                    C(2,2)          = _Cijkl(0,1,0,1,c,cinv); // C1212
                    C(1,0) = C(0,1) = _Cijkl(0,0,1,1,c,cinv); // C1122
                    C(2,0) = C(0,2) = _Cijkl(0,0,0,1,c,cinv); // C1112
                    C(2,1) = C(1,2) = _Cijkl(1,1,0,1,c,cinv); // C2212

                    break;
                }
                GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n");
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------------------------------------------------------------------
//                                          INCOMPRESSIBLE
// ---------------------------------------------------------------------------------------------------------------------------------
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi(const index_t i, const index_t j) const
{
    GISMO_ENSURE( ( (i < 3) && (j < 3) ) , "Index out of range. i="<<i<<", j="<<j);
    GISMO_ENSURE(!comp,"Material model is not incompressible?");
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return _dPsi_impl<mat>(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_impl(const index_t i, const index_t j) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    return 0.5 * mu * m_Gcon_ori(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_impl(const index_t i, const index_t j) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;
    T tmp;
    if ((i==2) && (j==2))
    {
        T traceCt =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                        m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                        m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                        m_Gcov_def(1,1)*m_Gcon_ori(1,1);
        tmp = c1/2.0 + c2 / 2.0 * traceCt;
    }
    else
        tmp =  c1 / 2. * m_Gcon_ori(i,j) + c2 / 2. * ( 1. / m_J0_sq * m_Gcon_ori(i,j) + m_J0_sq * m_Gcon_def(i,j) );

    return tmp;
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ENSURE( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(!comp,"Material model is not incompressible?");
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return _d2Psi_impl<mat>(i,j,k,l);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    return 0.0;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    T tmp;
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T c2 = mu/(m_parvals.at(2)+1);
    // T c1 = m_parvals.at(2)*c2;
    if      ( ((i==2) && (j==2)) && !((k==2) || (l==2)) ) // _dPsi/d22dkl
        tmp = c2 / 2.0 * m_Gcon_ori(k,l);
    else if ( !((i==2) && (j==2)) && ((k==2) || (l==2)) ) // _dPsi/dijd22
        tmp = c2 / 2.0 * m_Gcon_ori(i,j);
    else if ( ((i==2) && (j==2)) && ((k==2) || (l==2)) ) // _dPsi/d22d22
        tmp = 0.0;
    else
    {
        T Gabcd = - 1./2. * ( m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k) );
        tmp =  c2 / 2.0 * m_J0_sq * ( Gabcd + m_Gcon_def(i,j)*m_Gcon_def(k,l) );
    }
    return tmp;
}

// ---------------------------------------------------------------------------------------------------------------------------------
//                                          COMPRESSIBLE
// ---------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dI_1(const index_t i, const index_t j) const
{
    return m_Gcon_ori(i,j);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dI_2(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    return idelta(i,2)*idelta(j,2)*( c(2,2)*_dI_1(i,j) + m_J0_sq*cinv(i,j) ) + delta(i,2)*delta(j,2)*traceCt;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return _dPsi_impl<mat>(i,j,c,cinv);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    return mu/2.0 * math::pow(m_J_sq,-1./3.) * ( - 1.0/3.0 * I_1 * cinv(i,j) + _dI_1(i,j) ) + _dPsi_vol(i,j,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    T I_2 = c(2,2) * traceCt + m_J0_sq;
    T c2= mu/(m_parvals.at(2)+1);
    T c1= m_parvals.at(2)*c2;
    return  c1/2.0 * math::pow(m_J_sq,-1./3.) * ( - 1.0/3.0 * I_1 * cinv(i,j) + _dI_1(i,j) )
            + c2/2.0 * math::pow(m_J_sq,-2./3.) * ( - 2.0/3.0 * I_2 * cinv(i,j) + _dI_2(i,j,c,cinv) )
            + _dPsi_vol(i,j,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH_ext, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
    return mu / 2.0 * _dI_1(i,j) - mu / 2.0 * cinv(i,j) + lambda / 4.0 * ( m_J_sq - 1 ) * cinv(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_vol(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    return K * 0.25 * (m_J_sq - 1.0) * cinv(i,j);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(comp,"Material model is not compressible?");
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return _d2Psi_impl<mat>(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    return  1.0/9.0 * mu / 2.0 * math::pow(m_J_sq, -1.0/3.0) *  ( I_1*cinv(i,j)*cinv(k,l)
                                                            - 3.0*_dI_1(i,j)*cinv(k,l) - 3.0*cinv(i,j)*_dI_1(k,l)
                                                            - 3.0*I_1*dCinv  )
            + _d2Psi_vol(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    T I_2 = c(2,2) * traceCt + m_J0_sq;
    T d2I_2 = idelta(i,2)*idelta(j,2)*idelta(k,2)*idelta(l,2)*( m_J0_sq*( cinv(i,j)*cinv(k,l) + dCinv ) )
            + delta(i,2)*delta(j,2)*idelta(k,2)*idelta(l,2)*_dI_1(k,l)
            + idelta(i,2)*idelta(j,2)*delta(k,2)*delta(l,2)*_dI_1(i,j);
    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;
    // c1 = 0;
    return
          1.0/9.0 * c1 / 2.0 * math::pow(m_J_sq, -1.0/3.0) *  ( I_1*cinv(i,j)*cinv(k,l)
                                                            - 3.0*_dI_1(i,j)*cinv(k,l)       - 3.0*cinv(i,j)*_dI_1(k,l)
                                                            - 3.0*I_1*dCinv ) // + 9*d2I_1 = 0
        + 1.0/9.0 * c2 / 2.0 * math::pow(m_J_sq, -2.0/3.0) *  ( 4.0*I_2*cinv(i,j)*cinv(k,l) - 6.0*I_2*dCinv
                                                            - 6.0*_dI_2(i,j,c,cinv)*cinv(k,l)- 6.0*cinv(i,j)*_dI_2(k,l,c,cinv)
                                                            + 9.0*d2I_2 )
        + _d2Psi_vol(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH_ext, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
    return - mu / 2.0 * dCinv + lambda / 4.0 * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1.0)*dCinv );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_vol(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    return K * 0.25 * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1.0)*dCinv );
}

// ---------------------------------------------------------------------------------------------------------------------------------
//                                          STRETCHES
// ---------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da(const index_t a) const
{
    GISMO_ENSURE( a < 3 , "Index out of range. a="<<a);
    return _dPsi_da_impl<mat,comp>(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_impl(const index_t a) const
{
    GISMO_ENSURE( a < 3 , "Index out of range. a="<<a);
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T _dI_1a = 2*m_stretches(a);

    return  mu/2.0 * math::pow(m_J_sq,-1./3.) * ( -2./3. *  I_1 / m_stretches(a) + _dI_1a )
            + _dPsi_da_vol(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T _dI_1a = 2*m_stretches(a);
    return mu/2 * _dI_1a;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T _dI_1a = 2*m_stretches(a);
    T I_2   = math::pow(m_stretches(0),2.)*math::pow(m_stretches(1),2.)
            + math::pow(m_stretches(1),2.)*math::pow(m_stretches(2),2.)
            + math::pow(m_stretches(0),2.)*math::pow(m_stretches(2),2.);
    T _dI_2a  = 2*m_stretches(a)*( I_1 - math::pow(m_stretches(a),2.0) );

    T c2= mu/(m_parvals.at(2)+1);
    T c1= m_parvals.at(2)*c2;
    return c1/2.0 * math::pow(m_J_sq,-1./3.) * ( -2./3. *  I_1 / m_stretches(a) + _dI_1a )
         + c2/2.0 * math::pow(m_J_sq,-2./3.) * ( -4./3. *  I_2 / m_stretches(a) + _dI_2a )
         + _dPsi_da_vol(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T _dI_1a = 2*m_stretches(a);
    T _dI_2a  = 2*m_stretches(a)*( I_1 - math::pow(m_stretches(a),2.0) );

    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;
    return c1/2.0*_dI_1a + c2/2.0*_dI_2a;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::OG), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_impl(const index_t a) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    T tmp = 0.0;
    int n = (m_pars.size()-2)/2;
    T alpha_i, mu_i, Lambda;
    for (index_t k=0; k!=n; k++)
    {
        alpha_i = m_parvals.at(2*(k+1)+1);
        mu_i = m_parvals.at(2*(k+1));
        Lambda = math::pow(m_stretches(0),alpha_i) + math::pow(m_stretches(1),alpha_i) + math::pow(m_stretches(2),alpha_i);
        tmp += mu_i * math::pow(m_J_sq,-alpha_i/6.0) * ( math::pow(m_stretches(a),alpha_i-1) - 1./3. * 1./m_stretches(a) * Lambda );
    }
    return tmp + _dPsi_da_vol(a);
}
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::OG), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_impl(const index_t a) const
{
    T tmp = 0.0;
    int n = (m_pars.size()-2)/2;
    T alpha_i, mu_i;
    for (index_t k=0; k!=n; k++)
    {
        alpha_i = m_parvals.at(2*(k+1)+1);
        mu_i = m_parvals.at(2*(k+1));
        tmp += mu_i*math::pow(m_stretches(a),alpha_i-1);
    }
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::NH_ext), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T _dI_1a = 2*m_stretches(a);
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    //  choose compressibility function (and parameter)
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

    return mu/2.0 * _dI_1a - mu / m_stretches(a) + lambda / (m_stretches(a)*2) * (m_J_sq-1.0);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi_da_vol(const index_t a) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    T beta  = -2.0;
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    return K / (m_stretches(a)*beta) * (1.0 - math::pow(m_J_sq,-beta/2.0));
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab(const index_t a, const index_t b) const
{
    GISMO_ENSURE( ( (a < 3) && (b < 3) ) , "Index out of range. a="<<a<<", b="<<b);
    return _d2Psi_dab_impl<mat,comp>(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));

    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T _dI_1a = 2*m_stretches(a);
    T _dI_1b = 2*m_stretches(b);
    T d2I_1 = 2*delta(a,b);
    return  mu/2.0 * math::pow(m_J_sq,-1./3.) *   (
                                                    -2./3. * 1. / m_stretches(b) * ( -2./3. * I_1 / m_stretches(a) + _dI_1a )
                                                    -2./3. * 1. / m_stretches(a) * _dI_1b
                                                    +d2I_1
                                                    +2./3. * delta(a,b) * I_1 / (m_stretches(a)*m_stretches(a))
                                            )
            + _d2Psi_dab_vol(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T d2I_1 = 2*delta(a,b);
    return mu/2 * d2I_1;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));

    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T _dI_1a = 2*m_stretches(a);
    T _dI_1b = 2*m_stretches(b);
    T d2I_1 = 2*delta(a,b);

    T I_2   = math::pow(m_stretches(0),2.)*math::pow(m_stretches(1),2.)
            + math::pow(m_stretches(1),2.)*math::pow(m_stretches(2),2.)
            + math::pow(m_stretches(0),2.)*math::pow(m_stretches(2),2.);
    T _dI_2a = 2*m_stretches(a)*( I_1 - math::pow(m_stretches(a),2.0) );
    T _dI_2b = 2*m_stretches(b)*( I_1 - math::pow(m_stretches(b),2.0) );
    T d2I_2 = idelta(a,b)*4.0*m_stretches(a)*m_stretches(b) + delta(a,b)*2.0*(I_1 - m_stretches(a)*m_stretches(a));

    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;
    return
        c1/2.0 * math::pow(m_J_sq,-1./3.) *   (
                                                    -2./3. * 1. / m_stretches(b) * ( -2./3. * I_1 / m_stretches(a) + _dI_1a )
                                                    -2./3. * 1. / m_stretches(a) * _dI_1b
                                                    +d2I_1
                                                    +2./3. * delta(a,b) * I_1 / (m_stretches(a)*m_stretches(a))
                                            )
        + c2/2.0 * math::pow(m_J_sq,-2./3.) *   (
                                                    -4./3. * 1. / m_stretches(b) * ( -4./3. * I_2 / m_stretches(a) + _dI_2a )
                                                    -4./3. * 1. / m_stretches(a) * _dI_2b
                                                    +d2I_2
                                                    +4./3. * delta(a,b) * I_2 / (m_stretches(a)*m_stretches(a))
                                            )
        + _d2Psi_dab_vol(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_impl(const index_t a, const index_t b) const
{
    GISMO_ENSURE(m_pars.size()==3,"Mooney-Rivlin model needs to be a 3 parameter model");
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));

    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T d2I_1 = 2*delta(a,b);

    T d2I_2 = idelta(a,b)*4.0*m_stretches(a)*m_stretches(b) + delta(a,b)*2.0*(I_1 - m_stretches(a)*m_stretches(a));

    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;

    return c1/2.0 * d2I_1 + c2/2.0 * d2I_2;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::OG), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T tmp = 0.0;
    int n = (m_pars.size()-2)/2;
    T alpha_i, mu_i, Lambda;
    for (index_t k=0; k!=n; k++)
    {
        alpha_i = m_parvals.at(2*(k+1)+1);
        mu_i = m_parvals.at(2*(k+1));
        Lambda = math::pow(m_stretches(0),alpha_i) + math::pow(m_stretches(1),alpha_i) + math::pow(m_stretches(2),alpha_i);
        tmp += mu_i * math::pow(m_J_sq,-alpha_i/6.0) *
                (   - alpha_i/3. * ( math::pow(m_stretches(a),alpha_i-1.0) / m_stretches(b) + math::pow(m_stretches(b),alpha_i-1.0) / m_stretches(a)
                                    - 1./3. * 1. / (m_stretches(a)*m_stretches(b)) * Lambda )
                    + delta(a,b) * ( (alpha_i - 1.) * math::pow(m_stretches(a),alpha_i-2.0) + Lambda / 3. * math::pow(m_stretches(a),-2.0) )
                );
    }
    tmp += _d2Psi_dab_vol(a,b);
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::OG), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T tmp = 0.0;
    int n = (m_pars.size()-2)/2;
    T alpha_i, mu_i;
    for (index_t k=0; k!=n; k++)
    {
        alpha_i = m_parvals.at(2*(k+1)+1);
        mu_i = m_parvals.at(2*(k+1));
        tmp += mu_i*math::pow(m_stretches(a),alpha_i-2)*(alpha_i-1)*delta(a,b);
    }
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::NH_ext), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T d2I_1 = 2*delta(a,b);
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

    return mu/2.0 * d2I_1 + mu * delta(a,b) / ( m_stretches(a) * m_stretches(b) ) + lambda / (2*m_stretches(a)*m_stretches(b)) * ( 2*m_J_sq - delta(a,b) * (m_J_sq - 1.0) );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi_dab_vol(const index_t a, const index_t b) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    m_J_sq = math::pow(m_stretches(0)*m_stretches(1)*m_stretches(2),2.0);
    T beta  = -2.0;
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    return K / (beta*m_stretches(a)*m_stretches(b)) * ( beta*math::pow(m_J_sq,-beta/2.0) + delta(a,b) * (math::pow(m_J_sq,-beta/2.0) - 1.0) );

}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dJ_da(const index_t a) const
{
    return 1.0/m_stretches(a);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2J_dab(const index_t a, const index_t b) const
{
    if (a==b)
        return 0.0;
    else
        return 1.0  / ( m_stretches(a) * m_stretches(b) );
    // return ( 1.0 - delta(a,b) ) / ( m_stretches(a) * m_stretches(b) );
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_p() const
{
    return m_stretches(2) * _dPsi_da(2);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dp_da(const index_t a) const
{
    if (a==2)
        return m_stretches(2) * _d2Psi_dab(2,a) + _dPsi_da(2);
    else
        return m_stretches(2) * _d2Psi_dab(2,a);

    // return m_stretches(2) * _d2Psi_dab(2,a) + delta(a,2) * _dPsi_da(2);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sa(const index_t a) const
{
   return _Sa_impl<comp>(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sa_impl(const index_t a) const
{
    return 1.0/m_stretches(a) * _dPsi_da(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Sa_impl(const index_t a) const
{
    return 1.0/m_stretches(a) * (_dPsi_da(a) - _p() * _dJ_da(a) );
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dSa_db(const index_t a, const index_t b) const
{
    return _dSa_db_impl<comp>(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dSa_db_impl(const index_t a, const index_t b) const
{
    T tmp = 1.0/m_stretches(a) * _d2Psi_dab(a,b);
    if (a==b)
        tmp += - 1.0 / math::pow(m_stretches(a),2) * _dPsi_da(a);
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dSa_db_impl(const index_t a, const index_t b) const
{
    T tmp = 1.0/m_stretches(a) * ( _d2Psi_dab(a,b) - _dp_da(a)*_dJ_da(b) - _dp_da(b)*_dJ_da(a) - _p() * _d2J_dab(a,b) );
    if (a==b)
        tmp += - 1.0 / math::pow(m_stretches(a),2) * (_dPsi_da(a) - _p() * _dJ_da(a));
    return tmp;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cabcd(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    return _Cabcd_impl<comp>(a,b,c,d);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    // Compute part with stress tensor involved.
    T frac = 0.0;
    T tmp = 0.0;

    if (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-14)
    {
        // gsDebug<<"Stretches are equal; (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) = "<<abs((m_stretches(a) - m_stretches(b)) / m_stretches(a))<<"\n";
        frac = 1.0 / (2.0 * m_stretches(a) ) * ( _dSa_db(b,b) - _dSa_db(a,b));
    }
    else
        frac = ( _Sa(b)-_Sa(a) ) / (math::pow(m_stretches(b),2) - math::pow(m_stretches(a),2));

    GISMO_ENSURE( ( (a < 3) && (b < 3) && (c < 3) && (d < 3) ) , "Index out of range. a="<<a<<", b="<<b<<", c="<<c<<", d="<<d);
    if ( ( (a==b) && (c==d)) )
        tmp = 1/m_stretches(c) * _dSa_db(a,c);
    else if (( (a==d) && (b==c) && (a!=b) ) || ( ( (a==c) && (b==d) && (a!=b)) ))
        tmp = frac;
    // return 1/m_stretches(c) * _dSa_db(a,c) * delta(a,b) * delta(c,d) + frac * (delta(a,c)*delta(b,d) + delta(a,d)*delta(b,c)) * (1-delta(a,b));

    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    // Compute part with stress tensor involved.
    T frac = 0.0;
    T tmp = 0.0;

    if (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-14)
    {
        // gsDebug<<"Stretches are equal; (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) = "<<abs((m_stretches(a) - m_stretches(b)) / m_stretches(a))<<"\n";
        frac = 1.0 / (2.0 * m_stretches(a) ) * ( _dSa_db(b,b) - _dSa_db(a,b));
    }
    else
        frac = ( _Sa(b)-_Sa(a) ) / (math::pow(m_stretches(b),2) - math::pow(m_stretches(a),2));

    GISMO_ENSURE( ( (a < 2) && (b < 2) && (c < 2) && (d < 2) ) , "Index out of range. a="<<a<<", b="<<b<<", c="<<c<<", d="<<d);
    if ( ( (a==b) && (c==d)) )
        tmp = 1/m_stretches(c) * _dSa_db(a,c) + 1/(math::pow(m_stretches(a),2) * math::pow(m_stretches(c),2)) * ( math::pow(m_stretches(2),2) * _d2Psi_dab(2,2) + 2*_dPsi_da(2)*m_stretches(2) );
    else if (( (a==d) && (b==c) && (a!=b) ) || ( ( (a==c) && (b==d) && (a!=b)) ))
        tmp = frac;
    // return 1/m_stretches(c) * _dSa_db(a,c) * delta(a,b) * delta(c,d) + frac * (delta(a,c)*delta(b,d) + delta(a,d)*delta(b,c)) * (1-delta(a,b))
                // + delta(a,b)*delta(c,d)*1/(math::pow(m_stretches(a),2) * math::pow(m_stretches(c),2)) * ( math::pow(m_stretches(2),2) * _d2Psi_dab(2,2) + 2*_dPsi_da(2)*m_stretches(2) );

    return tmp;
}

//--------------------------------------------------------------------------------------------------------------------------------------

// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_dPsi(const index_t a) const
// {
//     T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));

//     GISMO_ENSURE( ( (a < 3) ) , "Index out of range. a="<<a);
//     GISMO_ENSURE(!comp,"Material model is not incompressible?");

//     if (m_material==9)
//     {
//         return mu * m_stretches(a,0);
//     }
//     else if (m_material==3)
//     {
//         // return 2.0*m_parvals.at(0)*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
//     }
//     else if (m_material==4)
//     {
//         gsMatrix<T> C(3,3);
//         C.setZero();
//         C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
//         C(2,2) = math::pow(m_J0,-2.0);
//         this->_computeStretch(C,m_gcon_ori);

//         return mu*m_stretches.at(a);
//     }
//     else
//         GISMO_ERROR("Material model not implemented.");
// }

// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::_d2Psi(const index_t a, const index_t b) const
// {
//     T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));

//     GISMO_ENSURE( ( (a < 3) && (b < 3) ) , "Index out of range. a="<<a<<", b="<<b);
//     GISMO_ENSURE(!comp,"Material model is not incompressible?");

//     if (m_material==9)
//     {
//         return ( a==b ? mu : 0.0 );
//     }
//     else if (m_material==3)
//     {
//         // return 2.0*m_parvals.at(0)*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
//     }
//     else if (m_material==4)
//     {
//         if (a==b)
//         {
//             return mu;
//         }
//         else
//             return 0.0;
//     }
//     else
//         GISMO_ERROR("Material model not implemented.");
// }

} // end namespace
