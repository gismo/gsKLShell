/** @file gsMaterialMatrix.h

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

#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixEval.h>

namespace gismo
{

// ?
int delta(const int a, const int b)
{
    return (a==b) ? 1 : 0;
}

int idelta(const int a, const int b)
{
    return (a!=b) ? 1 : 0;
}

// Linear material models
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const std::vector<gsFunction<T>*> &pars
                                        )
                                        :
                                        m_patches(&mp),
                                        m_thickness(&thickness),
                                        m_pars(pars),
                                        m_piece(nullptr)
{
    initialize();
}

// Linear material models
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                        const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & mp_def,
                                        const gsFunction<T> & thickness,
                                        const std::vector<gsFunction<T>*> &pars
                                        )
                                        :
                                        m_patches(&mp),
                                        m_defpatches(&mp_def),
                                        m_thickness(&thickness),
                                        m_pars(pars),
                                        m_piece(nullptr)
{
    initialize();
}

// Linear material models
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
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
                                    m_density(&density),
                                    m_piece(nullptr)
{
    initialize();
}

// Composite material model
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                        const gsFunctionSet<T>            & mp,
                                        const std::vector<T>                thickness,
                                        const std::vector<std::pair<T,T>> & YoungsModuli,
                                        const std::vector<T>              & ShearModuli,
                                        const std::vector<std::pair<T,T>> & PoissonRatios,
                                        const std::vector<T>                phis
                                        )
                                        :
                                        m_patches(&mp),
                                        m_YoungsModuli(YoungsModuli),
                                        m_ShearModuli(ShearModuli),
                                        m_PoissonRatios(PoissonRatios),
                                        m_thickValues(thickness),
                                        m_phis(phis),
                                        m_piece(nullptr)
{
    initialize();
    m_numPars=2;
    m_options.setInt("MaterialLaw",material_law::SvK_Orthotropic);
}
// Composite material model
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                        const gsFunctionSet<T>            & mp,
                                        const std::vector<T>                thickness,
                                        const std::vector<std::pair<T,T>> & YoungsModuli,
                                        const std::vector<T>              & ShearModuli,
                                        const std::vector<std::pair<T,T>> & PoissonRatios,
                                        const std::vector<T>                phis,
                                        const std::vector<T>                densities
                                        )
                                        :
                                        m_patches(&mp),
                                        m_YoungsModuli(YoungsModuli),
                                        m_ShearModuli(ShearModuli),
                                        m_PoissonRatios(PoissonRatios),
                                        m_thickValues(thickness),
                                        m_phis(phis),
                                        m_densities(densities),
                                        m_piece(nullptr)
{
    initialize();
    m_numPars=2;
    m_options.setInt("MaterialLaw",material_law::SvK_Orthotropic);
}

// Composite material model
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(
                                        const gsFunctionSet<T>            & mp,
                                        const gsFunctionSet<T>            & mp_def,
                                        const std::vector<T>                thickness,
                                        const std::vector<std::pair<T,T>> & YoungsModuli,
                                        const std::vector<T>              & ShearModuli,
                                        const std::vector<std::pair<T,T>> & PoissonRatios,
                                        const std::vector<T>                phis
                                        )
                                        :
                                        m_patches(&mp),
                                        m_defpatches(&mp_def),
                                        m_YoungsModuli(YoungsModuli),
                                        m_ShearModuli(ShearModuli),
                                        m_PoissonRatios(PoissonRatios),
                                        m_thickValues(thickness),
                                        m_phis(phis),
                                        m_piece(nullptr)
{
    initialize();
    m_numPars=2;
    m_options.setInt("MaterialLaw",material_law::SvK_Orthotropic);
}
// Composite material model
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(  const gsFunctionSet<T>            & mp,
                                        const gsFunctionSet<T>            & mp_def,
                                        const std::vector<T>                thickness,
                                        const std::vector<std::pair<T,T>> & YoungsModuli,
                                        const std::vector<T>              & ShearModuli,
                                        const std::vector<std::pair<T,T>> & PoissonRatios,
                                        const std::vector<T>                phis,
                                        const std::vector<T>                densities
                                        )
                                        :
                                        m_patches(&mp),
                                        m_defpatches(&mp_def),
                                        m_YoungsModuli(YoungsModuli),
                                        m_ShearModuli(ShearModuli),
                                        m_PoissonRatios(PoissonRatios),
                                        m_thickValues(thickness),
                                        m_phis(phis),
                                        m_densities(densities),
                                        m_piece(nullptr)
{
    initialize();
    m_numPars=2;
    m_options.setInt("MaterialLaw",material_law::SvK_Orthotropic);
}

// // Nonlinear material models
// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// gsMaterialMatrix<dim,T,matId,comp,mat,imp>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
//                                         const gsFunctionSet<T> & mp_def,
//                                         const gsFunction<T> & thickness,
//                                         const gsFunction<T> & par1,
//                                         const gsFunction<T> & par2
//                                         )
//                                         :
//                                         m_patches(&mp),
//                                         m_defpatches(&mp_def),
//                                         m_thickness(&thickness),
//                                         m_par1(&par1),
//                                         m_par2(&par2),
//                                         m_piece(nullptr)
// {
//     // m_material = 0;
//     m_moment = 0;
//     m_map.flags = m_map_def.flags =  NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;

//     // normal and deriv2 only needed when ComputeFull
// }

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::info() const
{
    gsInfo  <<"---------------------------------------------------------------------\n"
            <<"------------------------Material matrix class------------------------\n"
            <<"---------------------------------------------------------------------\n\n";

    gsInfo  <<"Material models:\n";
    gsInfo  <<"Specification of material models is done using parameters I and M \n"
            <<"combined and parameter C separately. E.g. for a Stretch-based \n"
            <<"implementation of the NH model, use 12 as material model. \n"
            <<"The conventions are:\n";

    gsInfo  <<"* Implementations (I):\n"
            <<"\t <none> \t Explicit implementation\n"
            <<"\t 1      \t Stretch-based implementation\n"
            <<"\t 2      \t General (using derivatives of Psi) implementation\n"
            <<"* Material models (M):\n"
            <<"\t 0      \t SvK Isotropic (Linear)\n"
            <<"\t 1      \t SvK Orthotropic (Linear)\n"
            <<"\t 2      \t NH Isotropic (Nonlinear)\n"
            <<"\t 3      \t MR Isotropic (Nonlinear)\n"
            <<"* Compressibility (C):\n"
            <<"\t 0      \t Incompressible\n"
            <<"\t 1      \t Compressible\n";
    gsInfo  <<"---------------------------------------------------------------------\n\n";

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::defaultOptions()
{
    m_options.addInt("MaterialLaw","Material law: \n"
                    "0 : St. Venant-Kirchhof Isotropic,\n"
                    "1 : St. Venant-Kirchhof Orthotropic,\n"
                    "2 : Neo-Hookean,\n"
                    "3 : Mooney-Rivlin,\n"
                    ,material_law::SvK_Isotropic);

    m_options.addInt("Compressibility","Specifies whether the material is modelled compressibile or incompressible",compressibility::incompressible);

    m_options.addInt("IntegrationMethod","Specifies thicknessintegration method; 0: Drectly Decoupled (DD); 1: Analytically Projected (AP); 2: Numerically Projected (NP)",integration::NP);

    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
short_t gsMaterialMatrix<dim,T,matId,comp,mat,imp>::domainDim() const { return 2; }

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
short_t gsMaterialMatrix<dim,T,matId,comp,mat,imp>::targetDim() const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getOptions() const
{

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::initialize()
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
    m_pIndex = 0;
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computePoints(const gsMatrix<T> & u, bool deformed) const
{
    gsMatrix<T> tmp;

    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    this->computeMetricUndeformed();

    if (m_defpatches->nPieces()!=0)
    {
        m_map_def.flags = m_map.flags;
        m_map_def.points = u;
        static_cast<const gsFunction<T>&>(m_defpatches->piece(m_pIndex)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
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

    computePoints_impl<mat>(u,deformed);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::OG, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computePoints_impl(const gsMatrix<T>& u, bool deformed) const
{
    T prod, sum, mu;
    for (index_t c=0; c!=m_parmat.cols(); c++)
    {
        for (index_t r=0; r!=m_parmat.rows(); r++)
            m_parvals.at(r) = m_parmat(r,c);

        mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
        GISMO_ENSURE((m_numPars-2 )% 2 ==0, "Ogden material models must have an even number of parameters (tuples of alpha_i and mu_i). m_numPars = "<< m_numPars);
        int n = (m_numPars-2)/2;
        sum = 0.0;
        for (index_t k=0; k!=n; k++)
        {
            prod = m_parvals.at(2*(k+1))*m_parvals.at(2*(k+1)+1);
            GISMO_ENSURE(prod > 0.0,"Product of coefficients must be positive for all indices");
            sum += prod;
        }
        GISMO_ENSURE((sum-2.*mu)/sum<1e-10,"Sum of products must be equal to 2*mu! sum = "<<sum<<"; 2*mu = "<<2.*mu);
    }
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat!=Material::OG, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computePoints_impl(const gsMatrix<T>& u, bool deformed) const
{

}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->getOptions();

    // This is for density and thickness output
    if (m_outputType==0)
        this->density_into(u,result);
    else if (m_outputType==1)
            this->eval_into_NP(u,result);
    else if (m_outputType==2)
            this->eval_into_NP(u,result);
    else if (m_outputType==9)
        this->stretch_into(u,result); // midplane stretches
    else if (m_outputType==10)
        this->eval_into_NP(u,result);
    else if (m_outputType==11)
        this->stretchDir_into(u,result); // stretch directions
    else
        GISMO_ERROR("Output type unknown");
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::density_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.flags = NEED_VALUE;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    result.resize(1, u.cols());
    m_thickness->eval_into(m_map.values[0], m_Tmat);
    m_density->eval_into(m_map.values[0], m_rhomat);
    for (index_t i = 0; i != u.cols(); ++i) // points
    {
        result(0,i) = m_Tmat(0,i)*m_rhomat(0,i);
    }

}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretch_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    this->computePoints(u);

    stretch_into_impl<comp>(u,result);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(3, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->getMetric(i,0.0); // on point i, with height 0.0

        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_J0_sq;

        res = evalStretch(C);
        result.col(i) = res.first;
    }
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(3, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,i);

        this->getMetric(i,0.0); // on point i, with height 0.0

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
        S33 = Sij(2,2,c,cinv);
        // S33_old = (S33 == 0.0) ? 1.0 : S33;
        C3333   = Cijkl3D(2,2,2,2,c,cinv);

        dc33 = -2. * S33 / C3333;
        for (index_t it = 0; it < itmax; it++)
        {
            c(2,2) += dc33;

            GISMO_ENSURE(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2) ;

            S33     = Sij(2,2,c,cinv);
            C3333   = Cijkl3D(2,2,2,2,c,cinv); //  or Cijkl???

            dc33 = -2. * S33 / C3333;
            if (abs(dc33) < tol)
            {
                res = evalStretch(c);
                result.col(i) = res.first;
                break;
            }
            GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, abs(dc33) = "<<abs(dc33)<<" and tolerance = "<<tol<<"\n");
        }
    }
}


template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretchDir_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    this->computePoints(u);

    stretchDir_into_impl<comp>(u,result);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        this->getMetric(i,0.0); // on point i, with height 0.0

        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = 1./m_J0_sq;

        res = evalStretch(C);
        result.col(i) = res.second.reshape(9,1);
        break;
    }
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(9, u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;
    for (index_t i=0; i!= u.cols(); i++)
    {
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,i);

        this->getMetric(i,0.0); // on point i, with height 0.0

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
        S33 = Sij(2,2,c,cinv);
        // S33_old = (S33 == 0.0) ? 1.0 : S33;
        C3333   = Cijkl3D(2,2,2,2,c,cinv);

        dc33 = -2. * S33 / C3333;
        for (index_t it = 0; it < itmax; it++)
        {
            c(2,2) += dc33;

            GISMO_ENSURE(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2) ;

            S33     = Sij(2,2,c,cinv);
            C3333   = Cijkl3D(2,2,2,2,c,cinv); //  or Cijkl???

            dc33 = -2. * S33 / C3333;
            if (abs(dc33) < tol)
            {
                res = evalStretch(c);
                result.col(i) = res.second.reshape(9,1);
                gsDebugVar(res.second.reshape(9,1));
                break;
            }
            GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, abs(dc33) = "<<abs(dc33)<<" and tolerance = "<<tol<<"\n");
        }
    }
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_into_NP(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    // Define the moment to take
    if      (m_outputType==1) // output is a vector
        if      (m_output==0) // N vector
            m_moment=0;
        else if (m_output==1) // M vector
            m_moment=1;
        else
            GISMO_ERROR("Something went wrong, m_output = "<<m_output);
    else if (m_outputType==2) // output is a matrix
        if      (m_output==0) // A matrix
            m_moment=0;
        else if (m_output==1) // B matrix
            m_moment=1;
        else if (m_output==2) // C matrix
            m_moment=1;
        else if (m_output==3) // D matrix
            m_moment=2;
        else
            GISMO_ERROR("Something went wrong, m_output = "<<m_output);
    else if (m_outputType==10) // output is a vector
        if      (m_output==0) // N vector
            m_moment=0;
        else if (m_output==1) // M vector
            m_moment=1;
        else
            GISMO_ERROR("Something went wrong, m_output = "<<m_output);
    else
        GISMO_ERROR("Something went wrong");

    // if (m_material==0)
    // {
    //     result = multiplyZ(u);
    // }
    // else
    result = integrateZ(u);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::thickness_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map.flags = NEED_VALUE;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    m_thickness->eval_into(m_map.values[0], result);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D(const index_t i, const gsMatrix<T> & z) const
{
    return eval3D_impl<mat,comp>(i,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_impl(const index_t i, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible(i, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_impl(const index_t i, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible(i, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_impl(const index_t i, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Compressible(i, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_matrix(const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    return eval3D_matrix_impl<mat,comp>(u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_matrix_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible_matrix(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_matrix_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible_matrix(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_matrix_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Compressible_matrix(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_vector(const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    return eval3D_vector_impl<mat,comp>(u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_vector_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible_vector(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_vector_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible_vector(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_vector_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Compressible_vector(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_pstress(const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    return eval3D_pstress_impl<mat,comp>(u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_pstress_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible_pstress(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_pstress_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Incompressible_pstress(u, z);
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && !(_mat==Material::SvK), gsMatrix<T>>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_pstress_impl(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    gsMatrix<T> result = eval_Compressible_pstress(u, z);
    return result;
}

// else if (m_material==-1) // test for integration
// {
//     result.resize(this->targetDim(), z.cols());
//     result.setZero();
//     for (index_t k=0; k!=z.cols(); k++)
//         for (index_t r = 0; r!=this->targetDim(); r++)
//         {
//             result(r,k) = math::pow(z(0,k),r);
//         }

// }

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D(const index_t i) const
{
    gsMatrix<T> z(1,1);
    z.setZero();
    return eval3D(i,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_matrix(const gsMatrix<T> & u) const
{
    gsMatrix<T> z(1,1);
    z.setZero();
    return eval3D_matrix(u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_vector(const gsMatrix<T> & u) const
{
    gsMatrix<T> z(1,1);
    z.setZero();
    return eval3D_vector(u,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval3D_pstress(const gsMatrix<T> & u) const
{
    gsMatrix<T> z(1,1);
    z.setZero();
    return eval3D_pstress(u,z);
}

// multiplies the material matrix by the thickness on all points
// NOTE: this function is a little outdated but it works in its current configuration (note also the implementation of the SvK stress in the Sij function).
// This function needs to be implemented in analytically projected integrals (see Roohbakashan & Sauer)
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::multiplyZ(const gsMatrix<T>& u) const
{
    this->computePoints(u);

    if ((m_output==1) && (m_outputType==1))
        m_moment = 2; // NEEDED SINCE m_moment=2 IS FOR THE OUTPUT OF THE M TENSOR, WHICH IN FACT HAS MOMENT 2. THIS IS BY CHOICE OF THE COMPUTATION OF THE STRAINS IN THE Sij() FUNCTION
    // Input: points in R2
    // Ouput: results in targetDim
    gsMatrix<T> result(this->targetDim(),u.cols());

    if (m_moment % 2 != 0)  //then the moment is odd
        result.setZero();
    else                    // then the moment is even
    {
        T fac;
        for (index_t i = 0; i != u.cols(); ++i) // points
        {
            m_evalPoints = this->eval3D(i);
            fac = 2.0/(m_moment+1) * math::pow( m_Tmat(0,i) / 2.0 , m_moment + 1);
            result.col(i) = m_evalPoints * fac;
        }

    }
    return result;
}

// integrates the material matrix over its last coordinate (thickness)
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::integrateZ(const gsMatrix<T>& u) const
{
    // Input: points in R2
    // Ouput: results in targetDim

    this->computePoints(u);

    if ((m_output==1) && (m_outputType==1))
        m_moment = 2; // NEEDED SINCE m_moment=2 IS FOR THE OUTPUT OF THE M TENSOR, WHICH IN FACT HAS MOMENT 2. THIS IS BY CHOICE OF THE COMPUTATION OF THE STRAINS IN THE Sij() FUNCTION

    // Perform integration
    m_numGauss = m_options.getInt("NumGauss");

    gsMatrix<T> result(9,1);
    result.resize(this->targetDim(),u.cols());
    result.setZero();

    // m_points.conservativeResize(m_points.rows()+1,Eigen::NoChange);

    T res;
    for (index_t j = 0; j != u.cols(); ++j) // for all points
    {
        gsGaussRule<T> gauss(m_numGauss);
        gsMatrix<T> pts(1,m_numGauss);
        gsMatrix<T> evalPoints;
        // m_points3D.resize(1,m_numGauss);

        gsMatrix<T> quNodes(1,m_numGauss);
        gsVector<T> quWeights(m_numGauss);

        // Compute values of in-plane points
        // m_points3D.topRows(u.rows()) = u.col(j).replicate(1,m_numGauss);
        // m_points3D.row(0) = u.col(j).replicate(1,m_numGauss);

        // set new integration point
        m_tHalf = m_Tmat(0,j)/2.0;
        gauss.mapTo(-m_tHalf,m_tHalf,quNodes,quWeights);
        pts.row(0) = quNodes;

        evalPoints = this->eval3D(j, pts);
        for (index_t i=0; i!=this->targetDim(); ++i) // components
        {
            res = 0.0;
            for (index_t k = 0; k != m_numGauss; ++k) // compute integral
                res += quWeights.at(k) * math::pow(quNodes(0,k),m_moment) * evalPoints(i,k);
            result(i,j) = res;
        }

    }
    return result;
}
//
// integrates the material matrix over its last coordinate (thickness)
// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::integrateZ(const gsMatrix<T>& u) const
// {
//     // Input: points in R2
//     // Ouput: results in targetDim

//     this->computePoints(u);

//     if ((m_output==1) && (m_outputType==1))
//         m_moment = 2; // NEEDED SINCE m_moment=2 IS FOR THE OUTPUT OF THE M TENSOR, WHICH IN FACT HAS MOMENT 2. THIS IS BY CHOICE OF THE COMPUTATION OF THE STRAINS IN THE Sij() FUNCTION

//     // Perform integration
//     m_numGauss = m_options.getInt("NumGauss");

//     gsMatrix<T> result(9,1);
//     result.resize(this->targetDim(),u.cols());
//     result.setZero();

//     m_gauss = gsGaussRule<T>(m_numGauss);
//     gsMatrix<T> z(m_numGauss,u.cols());
//     gsMatrix<T> w(m_numGauss,u.cols());
//     gsMatrix<T> vals;
//     // m_points3D.resize(1,m_numGauss);

//     // pre-compute Z
//     for (index_t k = 0; k != u.cols(); ++k) // for all points
//     {
//         gsMatrix<T> quNodes(1,m_numGauss);
//         gsVector<T> quWeights(m_numGauss);
//         // set new integration point
//         m_tHalf = m_Tmat(0,k)/2.0;
//         m_gauss.mapTo(-m_tHalf,m_tHalf,quNodes,quWeights);
//         w.col(k)=quWeights;
//         z.col(k)=quNodes.transpose();
//     }

//     vals = this->eval3D(u,z);

//     T res;
//     for (index_t k = 0; k != u.cols(); ++k) // for all points
//     {
//         for (index_t i=0; i!=this->targetDim(); ++i) // components
//         {
//             res = 0.0;
//             for (index_t j = 0; j != m_numGauss; ++j) // compute integral
//                 res += w(j,k) * math::pow(z(j,k),m_moment) * vals(i,j*u.cols() + k);
//             result(i,k) = res;
//         }

//     }
//     return result;
// }

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Incompressible(const index_t i, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    gsMatrix<T> result(this->targetDim(), z.cols());
    result.setZero();

    for( index_t j=0; j < z.cols(); ++j ) // through-thickness points
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,i);

        if (m_outputType==2)
        {
            // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
            this->getMetric(i,z.at(j)); // on point i, on height z(0,j)

            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j,3,3);
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
        else if (m_outputType==1)
        {
            // this->computeMetric(i,z.at(j),true,true);
            this->getMetric(i,z.at(j)); // on point i, on height z(0,j)

            result(0,j) = Sij(0,0);
            result(1,j) = Sij(1,1);
            result(2,j) = Sij(0,1);
        }
        else if (m_outputType==10)
        {

            // GISMO_ENSURE(m_material >= 10 && m_material < 20, "Only available for stretch-based materials.");

            // this->computeMetric(i,z.at(j),true,true);
            this->getMetric(i,z.at(j)); // on point i, on height z(0,j)

            result(0,j) = Sii(0);
            result(1,j) = Sii(1);
        }
        else
            GISMO_ERROR("no vector or matrix produced");
    }

    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Incompressible(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]


    this->computePoints(u);
    gsMatrix<T> result(this->targetDim(), u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            if (m_outputType==2)
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
            else if (m_outputType==1)
            {
                // this->computeMetric(i,z.at(j),true,true);
                this->getMetric(k,z(j,k)); // on point i, on height z(0,j)

                result(0,j*u.cols()+k) = Sij(0,0);
                result(1,j*u.cols()+k) = Sij(1,1);
                result(2,j*u.cols()+k) = Sij(0,1);
            }
            else if (m_outputType==10)
            {

                // GISMO_ENSURE(m_material >= 10 && m_material < 20, "Only available for stretch-based materials.");

                // this->computeMetric(i,z.at(j),true,true);
                this->getMetric(k,z(j,k)); // on point i, on height z(0,j)

                result(0,j*u.cols()+k) = Sii(0);
                result(1,j*u.cols()+k) = Sii(1);
            }
            else
                GISMO_ERROR("no vector or matrix produced");
        }
    }

    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Incompressible_matrix(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->computePoints(u);
    gsMatrix<T> result(this->targetDim(), u.cols() * z.rows());
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Incompressible_vector(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->computePoints(u);
    gsMatrix<T> result(this->targetDim(), u.cols() * z.rows());
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Incompressible_pstress(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->computePoints(u);
    gsMatrix<T> result(this->targetDim(), u.cols() * z.rows());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
                this->getMetric(k,z(j,k)); // on point i, on height z(0,j)

                result(0,j*u.cols()+k) = Sii(0);
                result(1,j*u.cols()+k) = Sii(1);
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
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ENSURE( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(!comp,"Material model is not incompressible?");

    return Cijkl_impl<mat,imp>(i,j,k,l);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    // --------------------------
    // Mooney-Rivlin
    // Parameter 3 is the ratio between c1 and c2.; c1 = m_parvals.at(2)*c2
    // --------------------------
    GISMO_ENSURE(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
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

    computeStretch(C);

    tmp = 0.0;
    for (index_t a = 0; a != 2; a++)
    {
        // C_iiii
        tmp +=  Cabcd(a,a,a,a)*(
                ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                );

        for (index_t b = a+1; b != 2; b++)
        {
            // C_iijj = C_jjii
            tmp +=  Cabcd(a,a,b,b)*(
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                        +
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                    );

            // C_ijij = Cjiji = Cijji = Cjiij
            tmp +=  Cabcd(a,b,a,b)*(
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    // --------------------------
    // General implementations
    // --------------------------
    return 4.0 * d2Psi(i,j,k,l) + 4.0 * d2Psi(2,2,2,2)*math::pow(m_J0_sq,-2.0)*m_Gcon_def(i,j)*m_Gcon_def(k,l)
            - 4.0/ m_J0_sq  * ( d2Psi(2,2,i,j)*m_Gcon_def(k,l) + d2Psi(2,2,k,l)*m_Gcon_def(i,j) )
            + 2.0 * dPsi(2,2) / m_J0_sq * (2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k));
}

// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// template <enum Material _mat>
// // OTHER CASES!
// typename std::enable_if<(_mat >= 30) && !(_mat==0) && !(_mat==2) && !(_mat==3), T>::type
// gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
// {
//     GISMO_ERROR("Material model unknown (model = "<<_mat<<"). Use gsMaterialMatrix<dim,T,matId,comp,mat,imp>::info() to see the options.");
// }

        // else if (m_material==4)
        //         GISMO_ERROR("Material model 4 is not invariant-based! Use 14 instead...");
        // else if (m_material==5)
        //         GISMO_ERROR("Material model 5 is only for compressible materials...");


// Condensation of the 3D tensor for compressible materials
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE(c.cols()==c.rows(),"Matrix c must be square");
    GISMO_ENSURE(c.cols()==3,"Matrix c must be 3x3");
    GISMO_ENSURE(cinv.cols()==cinv.rows(),"Matrix cinv must be square");
    GISMO_ENSURE(cinv.cols()==3,"Matrix cinv must be 3x3");
    GISMO_ENSURE( ( (i <2) && (j <2) && (k <2) && (l <2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(comp,"Material model is not compressible?");

    return Cijkl_impl<imp>(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Implementation _imp>
typename std::enable_if< _imp==Implementation::Spectral, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // static condensation is done before the projection
    computeStretch(c);
    return Cijkl3D(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Implementation _imp>
typename std::enable_if<!(_imp==Implementation::Spectral), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    return Cijkl3D(i,j,k,l,c,cinv) - ( Cijkl3D(i,j,2,2,c,cinv) * Cijkl3D(2,2,k,l,c,cinv) ) / Cijkl3D(2,2,2,2,c,cinv);
}

// 3D tensor for compressible materials
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE( ( (i <3) && (j <3) && (k <3) && (l <3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    return Cijkl3D_impl<mat,imp>(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::SvK && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ERROR("Compressible material matrix requested, but not needed. How?");
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Mooney-Rivlin
    // --------------------------
    GISMO_ENSURE(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");

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
            + delta(i,2)*delta(j,2)*idelta(k,2)*idelta(l,2)*dI_1(k,l)
            + idelta(i,2)*idelta(j,2)*delta(k,2)*delta(l,2)*dI_1(i,j);
    T c2 = mu/(m_parvals.at(2) + 1);
    T c1 = m_parvals.at(2)*c2;

    return  1.0/9.0 * c1 * math::pow(m_J_sq, -1.0/3.0) *  ( 2.0*I_1*cinv(i,j)*cinv(k,l) - 6.0*I_1*dCinv
                                                            - 6.0*dI_1(i,j)*cinv(k,l)     - 6.0*cinv(i,j)*dI_1(k,l) ) // + 9*d2I_1 = 0
            + 1.0/9.0 * c2 * math::pow(m_J_sq, -2.0/3.0) *  ( 8.0*I_2*cinv(i,j)*cinv(k,l) - 12.0*I_2*dCinv
                                                                - 12.0*dI_2(i,j,c,cinv)*cinv(k,l)- 12.0*cinv(i,j)*dI_2(k,l,c,cinv)
                                                                + 18.0*d2I_2 )
            + K * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1)*dCinv );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH_ext && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Stretch-based implementations
    // --------------------------
    T tmp;
    if ( (i==2) && (j==2) && (k==2) && (l==2) ) // if C3333 (no static condensation)
        tmp = Cabcd(2,2,2,2);
    else
    {
        tmp = 0.0;
        T C = 0.0;
        T C2222 = Cabcd(2,2,2,2);
        // T Cab22,C22ab;
        for (index_t a = 0; a != 2; a++)
        {
            // C_iiii
            // if (!((i==2 || j==2 || k==2 || l==2) && a!=2))
            // C = Cabcd(a,a,a,a);
            C = Cabcd(a,a,a,a) - math::pow(Cabcd(2,2,a,a),2) / C2222;
            tmp +=  C*(
                        ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                        ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                    );

            for (index_t b = a+1; b != 2; b++)
            {
                // C_iijj
                // C = Cabcd(a,a,b,b);
                C = Cabcd(a,a,b,b) - Cabcd(a,a,2,2) * Cabcd(2,2,b,b) / C2222;
                tmp +=  C*(
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(b)) )
                            +
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(b)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                            ( m_gcon_ori.col(k).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(a)) )
                        );

                // C_ijij = Cjiji = Cijji = Cjiij
                // C = Cabcd(a,b,a,b);
                C = Cabcd(a,b,a,b) - math::pow(Cabcd(2,2,a,b),2) / C2222;
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // General implementations
    // --------------------------
    return 4.0 * d2Psi(i,j,k,l,c,cinv);
}

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
// T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij(const index_t i, const index_t j) const { Sij(i,j,NULL,NULL); }

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij(const index_t i, const index_t j) const
{
    return Sij_impl<mat,imp>(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::SvK && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j) const
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

    return    Cijkl(i,j,0,0) * stress(0,0) + Cijkl(i,j,0,1) * stress(0,1)
            + Cijkl(i,j,1,0) * stress(1,0) + Cijkl(i,j,1,1) * stress(1,1);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j) const
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j) const
{
    // --------------------------
    // Mooney-Rivlin
    // Parameter 3 is the ratio between c1 and c2.
    // --------------------------
    GISMO_ENSURE(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j) const
{
    T tmp = 0.0;
    gsMatrix<T> C(3,3);
    C.setZero();
    C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
    C(2,2) = 1./m_J0_sq;

    computeStretch(C);

    for (index_t a = 0; a != 2; a++)
    {
        tmp += Sa(a)*(
                    ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )
                    );
    }
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Generalized, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j) const
{
    // --------------------------
    // Generalized
    // --------------------------
    return 2.0 * dPsi(i,j) - 2.0 * dPsi(2,2) * math::pow(m_J0_sq,-1.0)*m_Gcon_def(i,j);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    return Sij_impl<mat,imp>(i,j,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Mooney-Rivlin
    // Parameter 3 is the ratio between c1 and c2.
    // --------------------------
    GISMO_ENSURE(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
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
            + c2 * math::pow( m_J_sq , -2.0/3.0 ) * ( dI_2(i,j,c,cinv)- 2.0/3.0 * I_2 * cinv(i,j) )
            + K * 0.5 * ( m_J_sq - 1.0 ) * cinv(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_mat==Material::NH_ext && _imp == Implementation::Analytical, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T tmp = 0.0;
    computeStretch(c);
    for (index_t a = 0; a != 3; a++)
    {
        tmp += Sa(a)*(
                    ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )
                    );
    }
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, enum Implementation _imp>
typename std::enable_if<_imp == Implementation::Generalized, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    // --------------------------
    // Generalized
    // --------------------------
    return 2.0 * dPsi(i,j,c,cinv);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sii(const index_t i) const // principle stresses
{
    gsMatrix<T> C(3,3);
    C.setZero();
    C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
    C(2,2) = 1./m_J0_sq;
    return Sii(i,C);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sii(const index_t i, const gsMatrix<T> & c) const
{
    computeStretch(c);
    return Sa(i);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Compressible(const index_t i, const gsMatrix<T>& z) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    gsMatrix<T> result(this->targetDim(), z.cols());
    result.setZero();

    for( index_t j=0; j < z.cols(); ++j ) // through thickness
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,i);

        // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
        this->getMetric(i,z.at(j)); // on point i, on height z(0,j)

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
        S33 = Sij(2,2,c,cinv);
        // S33_old = (S33 == 0.0) ? 1.0 : S33;
        C3333   = Cijkl3D(2,2,2,2,c,cinv);

        dc33 = -2. * S33 / C3333;
        for (index_t it = 0; it < itmax; it++)
        {
            c(2,2) += dc33;

            GISMO_ENSURE(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
            cinv(2,2) = 1.0/c(2,2);

            m_J_sq = m_J0_sq * c(2,2) ;

            S33     = Sij(2,2,c,cinv);
            C3333   = Cijkl3D(2,2,2,2,c,cinv); //  or Cijkl???

            dc33 = -2. * S33 / C3333;
            // if (abs(S33/S33_old) < tol)
            if (abs(dc33) < tol)
            {
                // gsInfo<<"Converged in "<<it<<" iterations, abs(S33) = "<<abs(dc33)<<" and tolerance = "<<tol<<"\n";
                if (m_outputType==2)
                    {
                        gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j,3,3);
                        /*
                            C = C1111,  C1122,  C1112
                                symm,   C2222,  C2212
                                symm,   symm,   C1212
                        */
                        C(0,0)          = Cijkl(0,0,0,0,c,cinv); // C1111
                        C(1,1)          = Cijkl(1,1,1,1,c,cinv); // C2222
                        C(2,2)          = Cijkl(0,1,0,1,c,cinv); // C1212
                        C(1,0) = C(0,1) = Cijkl(0,0,1,1,c,cinv); // C1122
                        C(2,0) = C(0,2) = Cijkl(0,0,0,1,c,cinv); // C1112
                        C(2,1) = C(1,2) = Cijkl(1,1,0,1,c,cinv); // C2212
                    }
                    else if (m_outputType==1)
                    {
                        result(0,j) = Sij(0,0,c,cinv); // S11
                        result(1,j) = Sij(1,1,c,cinv); // S22
                        result(2,j) = Sij(0,1,c,cinv); // S12
                    }
                    else if (m_outputType==10)
                    {
                        GISMO_ENSURE(imp==Implementation::Spectral, "Only available for stretch-based materials.");

                        result(0,j) = Sii(0,c); // S11
                        result(1,j) = Sii(1,c); // S22
                    }
                else
                    GISMO_ERROR("no vector or matrix produced");

                    break;
            }
            GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n");
        }
    }
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Compressible_matrix(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    this->computePoints(u);
    gsMatrix<T> result(9, u.cols() * z.cols());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.cols(); ++j ) // through-thickness points
        {
            // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
            this->getMetric(k,z.at(j)); // on point i, on height z(0,j)

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
            S33 = Sij(2,2,c,cinv);
            // S33_old = (S33 == 0.0) ? 1.0 : S33;
            C3333   = Cijkl3D(2,2,2,2,c,cinv);

            dc33 = -2. * S33 / C3333;
            for (index_t it = 0; it < itmax; it++)
            {
                c(2,2) += dc33;

                GISMO_ENSURE(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
                cinv(2,2) = 1.0/c(2,2);

                m_J_sq = m_J0_sq * c(2,2) ;

                S33     = Sij(2,2,c,cinv);
                C3333   = Cijkl3D(2,2,2,2,c,cinv); //  or Cijkl???

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
                    C(0,0)          = Cijkl(0,0,0,0,c,cinv); // C1111
                    C(1,1)          = Cijkl(1,1,1,1,c,cinv); // C2222
                    C(2,2)          = Cijkl(0,1,0,1,c,cinv); // C1212
                    C(1,0) = C(0,1) = Cijkl(0,0,1,1,c,cinv); // C1122
                    C(2,0) = C(0,2) = Cijkl(0,0,0,1,c,cinv); // C1112
                    C(2,1) = C(1,2) = Cijkl(1,1,0,1,c,cinv); // C2212

                    break;
                }
                GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n");
            }
        }
    }
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Compressible_vector(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    this->computePoints(u);
    gsMatrix<T> result(3, u.cols() * z.cols());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.cols(); ++j ) // through-thickness points
        {
            // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
            this->getMetric(k,z.at(j)); // on point i, on height z(0,j)

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
            S33 = Sij(2,2,c,cinv);
            // S33_old = (S33 == 0.0) ? 1.0 : S33;
            C3333   = Cijkl3D(2,2,2,2,c,cinv);

            dc33 = -2. * S33 / C3333;
            for (index_t it = 0; it < itmax; it++)
            {
                c(2,2) += dc33;

                GISMO_ENSURE(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
                cinv(2,2) = 1.0/c(2,2);

                m_J_sq = m_J0_sq * c(2,2) ;

                S33     = Sij(2,2,c,cinv);
                C3333   = Cijkl3D(2,2,2,2,c,cinv); //  or Cijkl???

                dc33 = -2. * S33 / C3333;
                // if (abs(S33/S33_old) < tol)
                if (abs(dc33) < tol)
                {

                    result(0,j*u.cols()+k) = Sij(0,0,c,cinv); // S11
                    result(1,j*u.cols()+k) = Sij(1,1,c,cinv); // S22
                    result(2,j*u.cols()+k) = Sij(0,1,c,cinv); // S12

                    break;
                }
                GISMO_ENSURE(it != itmax-1,"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n");
            }
        }
    }
    return result;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
gsMatrix<T> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::eval_Compressible_pstress(const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    this->computePoints(u);
    gsMatrix<T> result(2, u.cols() * z.cols());
    result.setZero();

    for (index_t k=0; k!=u.cols(); k++)
    {
        // Evaluate material properties on the quadrature point
        for (index_t v=0; v!=m_parmat.rows(); v++)
            m_parvals.at(v) = m_parmat(v,k);

        for( index_t j=0; j < z.cols(); ++j ) // through-thickness points
        {
            // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
            this->getMetric(k,z.at(j)); // on point i, on height z(0,j)

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
            S33 = Sij(2,2,c,cinv);
            // S33_old = (S33 == 0.0) ? 1.0 : S33;
            C3333   = Cijkl3D(2,2,2,2,c,cinv);

            dc33 = -2. * S33 / C3333;
            for (index_t it = 0; it < itmax; it++)
            {
                c(2,2) += dc33;

                GISMO_ENSURE(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
                cinv(2,2) = 1.0/c(2,2);

                m_J_sq = m_J0_sq * c(2,2) ;

                S33     = Sij(2,2,c,cinv);
                C3333   = Cijkl3D(2,2,2,2,c,cinv); //  or Cijkl???

                dc33 = -2. * S33 / C3333;
                // if (abs(S33/S33_old) < tol)
                if (abs(dc33) < tol)
                {
                    GISMO_ENSURE(imp==Implementation::Spectral, "Only available for stretch-based materials.");

                    result(0,j*u.cols()+k) = Sii(0,c); // S11
                    result(1,j*u.cols()+k) = Sii(1,c); // S22

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
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi(const index_t i, const index_t j) const
{
    GISMO_ENSURE( ( (i < 3) && (j < 3) ) , "Index out of range. i="<<i<<", j="<<j);
    GISMO_ENSURE(!comp,"Material model is not incompressible?");
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return dPsi_impl<mat>(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_impl(const index_t i, const index_t j) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    return 0.5 * mu * m_Gcon_ori(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_impl(const index_t i, const index_t j) const
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
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ENSURE( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(!comp,"Material model is not incompressible?");
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return d2Psi_impl<mat>(i,j,k,l);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    return 0.0;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    T tmp;
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T c2 = mu/(m_parvals.at(2)+1);
    // T c1 = m_parvals.at(2)*c2;
    if      ( ((i==2) && (j==2)) && !((k==2) || (l==2)) ) // dPsi/d22dkl
        tmp = c2 / 2.0 * m_Gcon_ori(k,l);
    else if ( !((i==2) && (j==2)) && ((k==2) || (l==2)) ) // dPsi/dijd22
        tmp = c2 / 2.0 * m_Gcon_ori(i,j);
    else if ( ((i==2) && (j==2)) && ((k==2) || (l==2)) ) // dPsi/d22d22
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
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dI_1(const index_t i, const index_t j) const
{
    return m_Gcon_ori(i,j);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dI_2(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    return idelta(i,2)*idelta(j,2)*( c(2,2)*dI_1(i,j) + m_J0_sq*cinv(i,j) ) + delta(i,2)*delta(j,2)*traceCt;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return dPsi_impl<mat>(i,j,c,cinv);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    return mu/2.0 * math::pow(m_J_sq,-1./3.) * ( - 1.0/3.0 * I_1 * cinv(i,j) + dI_1(i,j) ) + dPsi_vol(i,j,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
    return  c1/2.0 * math::pow(m_J_sq,-1./3.) * ( - 1.0/3.0 * I_1 * cinv(i,j) + dI_1(i,j) )
            + c2/2.0 * math::pow(m_J_sq,-2./3.) * ( - 2.0/3.0 * I_2 * cinv(i,j) + dI_2(i,j,c,cinv) )
            + dPsi_vol(i,j,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH_ext, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
    return mu / 2.0 * dI_1(i,j) - mu / 2.0 * cinv(i,j) + lambda / 4.0 * ( m_J_sq - 1 ) * cinv(i,j);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_vol(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    return K * 0.25 * (m_J_sq - 1.0) * cinv(i,j);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ENSURE( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ENSURE(comp,"Material model is not compressible?");
    GISMO_ENSURE(imp==Implementation::Generalized,"Not generalized implementation");
    return d2Psi_impl<mat>(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
                                                            - 3.0*dI_1(i,j)*cinv(k,l) - 3.0*cinv(i,j)*dI_1(k,l)
                                                            - 3.0*I_1*dCinv  )
            + d2Psi_vol(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::MR, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
            + delta(i,2)*delta(j,2)*idelta(k,2)*idelta(l,2)*dI_1(k,l)
            + idelta(i,2)*idelta(j,2)*delta(k,2)*delta(l,2)*dI_1(i,j);
    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;
    // c1 = 0;
    return
          1.0/9.0 * c1 / 2.0 * math::pow(m_J_sq, -1.0/3.0) *  ( I_1*cinv(i,j)*cinv(k,l)
                                                            - 3.0*dI_1(i,j)*cinv(k,l)       - 3.0*cinv(i,j)*dI_1(k,l)
                                                            - 3.0*I_1*dCinv ) // + 9*d2I_1 = 0
        + 1.0/9.0 * c2 / 2.0 * math::pow(m_J_sq, -2.0/3.0) *  ( 4.0*I_2*cinv(i,j)*cinv(k,l) - 6.0*I_2*dCinv
                                                            - 6.0*dI_2(i,j,c,cinv)*cinv(k,l)- 6.0*cinv(i,j)*dI_2(k,l,c,cinv)
                                                            + 9.0*d2I_2 )
        + d2Psi_vol(i,j,k,l,c,cinv);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat>
typename std::enable_if<_mat==Material::NH_ext, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
    return - mu / 2.0 * dCinv + lambda / 4.0 * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1.0)*dCinv );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_vol(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
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
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da(const index_t a) const
{
    GISMO_ENSURE( a < 3 , "Index out of range. a="<<a);
    return dPsi_da_impl<mat,comp>(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_impl(const index_t a) const
{
    GISMO_ENSURE( a < 3 , "Index out of range. a="<<a);
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T dI_1a = 2*m_stretches(a);

    return  mu/2.0 * math::pow(m_J_sq,-1./3.) * ( -2./3. *  I_1 / m_stretches(a) + dI_1a )
            + dPsi_da_vol(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T dI_1a = 2*m_stretches(a);
    return mu/2 * dI_1a;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T dI_1a = 2*m_stretches(a);
    T I_2   = math::pow(m_stretches(0),2.)*math::pow(m_stretches(1),2.)
            + math::pow(m_stretches(1),2.)*math::pow(m_stretches(2),2.)
            + math::pow(m_stretches(0),2.)*math::pow(m_stretches(2),2.);
    T dI_2a  = 2*m_stretches(a)*( I_1 - math::pow(m_stretches(a),2.0) );

    T c2= mu/(m_parvals.at(2)+1);
    T c1= m_parvals.at(2)*c2;
    return c1/2.0 * math::pow(m_J_sq,-1./3.) * ( -2./3. *  I_1 / m_stretches(a) + dI_1a )
         + c2/2.0 * math::pow(m_J_sq,-2./3.) * ( -4./3. *  I_2 / m_stretches(a) + dI_2a )
         + dPsi_da_vol(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T dI_1a = 2*m_stretches(a);
    T dI_2a  = 2*m_stretches(a)*( I_1 - math::pow(m_stretches(a),2.0) );

    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;
    return c1/2.0*dI_1a + c2/2.0*dI_2a;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::OG), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_impl(const index_t a) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    T tmp = 0.0;
    int n = (m_numPars-2)/2;
    T alpha_i, mu_i, Lambda;
    for (index_t k=0; k!=n; k++)
    {
        alpha_i = m_parvals.at(2*(k+1)+1);
        mu_i = m_parvals.at(2*(k+1));
        Lambda = math::pow(m_stretches(0),alpha_i) + math::pow(m_stretches(1),alpha_i) + math::pow(m_stretches(2),alpha_i);
        tmp += mu_i * math::pow(m_J_sq,-alpha_i/6.0) * ( math::pow(m_stretches(a),alpha_i-1) - 1./3. * 1./m_stretches(a) * Lambda );
    }
    return tmp + dPsi_da_vol(a);
}
template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::OG), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_impl(const index_t a) const
{
    T tmp = 0.0;
    int n = (m_numPars-2)/2;
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_impl(const index_t a) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T dI_1a = 2*m_stretches(a);
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    //  choose compressibility function (and parameter)
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

    return mu/2.0 * dI_1a - mu / m_stretches(a) + lambda / (m_stretches(a)*2) * (m_J_sq-1.0);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi_da_vol(const index_t a) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    T beta  = -2.0;
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    return K / (m_stretches(a)*beta) * (1.0 - math::pow(m_J_sq,-beta/2.0));
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab(const index_t a, const index_t b) const
{
    GISMO_ENSURE( ( (a < 3) && (b < 3) ) , "Index out of range. a="<<a<<", b="<<b);
    return d2Psi_dab_impl<mat,comp>(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));

    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T dI_1a = 2*m_stretches(a);
    T dI_1b = 2*m_stretches(b);
    T d2I_1 = 2*delta(a,b);
    return  mu/2.0 * math::pow(m_J_sq,-1./3.) *   (
                                                    -2./3. * 1. / m_stretches(b) * ( -2./3. * I_1 / m_stretches(a) + dI_1a )
                                                    -2./3. * 1. / m_stretches(a) * dI_1b
                                                    +d2I_1
                                                    +2./3. * delta(a,b) * I_1 / (m_stretches(a)*m_stretches(a))
                                            )
            + d2Psi_dab_vol(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::NH), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T d2I_1 = 2*delta(a,b);
    return mu/2 * d2I_1;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));

    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T dI_1a = 2*m_stretches(a);
    T dI_1b = 2*m_stretches(b);
    T d2I_1 = 2*delta(a,b);

    T I_2   = math::pow(m_stretches(0),2.)*math::pow(m_stretches(1),2.)
            + math::pow(m_stretches(1),2.)*math::pow(m_stretches(2),2.)
            + math::pow(m_stretches(0),2.)*math::pow(m_stretches(2),2.);
    T dI_2a = 2*m_stretches(a)*( I_1 - math::pow(m_stretches(a),2.0) );
    T dI_2b = 2*m_stretches(b)*( I_1 - math::pow(m_stretches(b),2.0) );
    T d2I_2 = idelta(a,b)*4.0*m_stretches(a)*m_stretches(b) + delta(a,b)*2.0*(I_1 - m_stretches(a)*m_stretches(a));

    T c2 = mu/(m_parvals.at(2)+1);
    T c1 = m_parvals.at(2)*c2;
    return
        c1/2.0 * math::pow(m_J_sq,-1./3.) *   (
                                                    -2./3. * 1. / m_stretches(b) * ( -2./3. * I_1 / m_stretches(a) + dI_1a )
                                                    -2./3. * 1. / m_stretches(a) * dI_1b
                                                    +d2I_1
                                                    +2./3. * delta(a,b) * I_1 / (m_stretches(a)*m_stretches(a))
                                            )
        + c2/2.0 * math::pow(m_J_sq,-2./3.) *   (
                                                    -4./3. * 1. / m_stretches(b) * ( -4./3. * I_2 / m_stretches(a) + dI_2a )
                                                    -4./3. * 1. / m_stretches(a) * dI_2b
                                                    +d2I_2
                                                    +4./3. * delta(a,b) * I_2 / (m_stretches(a)*m_stretches(a))
                                            )
        + d2Psi_dab_vol(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::MR), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_impl(const index_t a, const index_t b) const
{
    GISMO_ENSURE(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T tmp = 0.0;
    int n = (m_numPars-2)/2;
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
    tmp += d2Psi_dab_vol(a,b);
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <enum Material _mat, bool _comp>
typename std::enable_if<!_comp && (_mat==Material::OG), T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T tmp = 0.0;
    int n = (m_numPars-2)/2;
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
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_impl(const index_t a, const index_t b) const
{
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T d2I_1 = 2*delta(a,b);
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

    return mu/2.0 * d2I_1 + mu * delta(a,b) / ( m_stretches(a) * m_stretches(b) ) + lambda / (2*m_stretches(a)*m_stretches(b)) * ( 2*m_J_sq - delta(a,b) * (m_J_sq - 1.0) );
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi_dab_vol(const index_t a, const index_t b) const
{
    GISMO_ENSURE(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
    m_J_sq = math::pow(m_stretches(0)*m_stretches(1)*m_stretches(2),2.0);
    T beta  = -2.0;
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    return K / (beta*m_stretches(a)*m_stretches(b)) * ( beta*math::pow(m_J_sq,-beta/2.0) + delta(a,b) * (math::pow(m_J_sq,-beta/2.0) - 1.0) );

}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dJ_da(const index_t a) const
{
    return 1.0/m_stretches(a);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2J_dab(const index_t a, const index_t b) const
{
    if (a==b)
        return 0.0;
    else
        return 1.0  / ( m_stretches(a) * m_stretches(b) );
    // return ( 1.0 - delta(a,b) ) / ( m_stretches(a) * m_stretches(b) );
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::p() const
{
    return m_stretches(2) * dPsi_da(2);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dp_da(const index_t a) const
{
    if (a==2)
        return m_stretches(2) * d2Psi_dab(2,a) + dPsi_da(2);
    else
        return m_stretches(2) * d2Psi_dab(2,a);

    // return m_stretches(2) * d2Psi_dab(2,a) + delta(a,2) * dPsi_da(2);
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sa(const index_t a) const
{
   return Sa_impl<comp>(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sa_impl(const index_t a) const
{
    return 1.0/m_stretches(a) * dPsi_da(a);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Sa_impl(const index_t a) const
{
    return 1.0/m_stretches(a) * (dPsi_da(a) - p() * dJ_da(a) );
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dSa_db(const index_t a, const index_t b) const
{
    return dSa_db_impl<comp>(a,b);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dSa_db_impl(const index_t a, const index_t b) const
{
    T tmp = 1.0/m_stretches(a) * d2Psi_dab(a,b);
    if (a==b)
        tmp += - 1.0 / math::pow(m_stretches(a),2) * dPsi_da(a);
    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dSa_db_impl(const index_t a, const index_t b) const
{
    T tmp = 1.0/m_stretches(a) * ( d2Psi_dab(a,b) - dp_da(a)*dJ_da(b) - dp_da(b)*dJ_da(a) - p() * d2J_dab(a,b) );
    if (a==b)
        tmp += - 1.0 / math::pow(m_stretches(a),2) * (dPsi_da(a) - p() * dJ_da(a));
    return tmp;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cabcd(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    return Cabcd_impl<comp>(a,b,c,d);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    // Compute part with stress tensor involved.
    T frac = 0.0;
    T tmp = 0.0;

    if (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-14)
    {
        // gsDebug<<"Stretches are equal; (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) = "<<abs((m_stretches(a) - m_stretches(b)) / m_stretches(a))<<"\n";
        frac = 1.0 / (2.0 * m_stretches(a) ) * ( dSa_db(b,b) - dSa_db(a,b));
    }
    else
        frac = ( Sa(b)-Sa(a) ) / (math::pow(m_stretches(b),2) - math::pow(m_stretches(a),2));

    GISMO_ENSURE( ( (a < 3) && (b < 3) && (c < 3) && (d < 3) ) , "Index out of range. a="<<a<<", b="<<b<<", c="<<c<<", d="<<d);
    if ( ( (a==b) && (c==d)) )
        tmp = 1/m_stretches(c) * dSa_db(a,c);
    else if (( (a==d) && (b==c) && (a!=b) ) || ( ( (a==c) && (b==d) && (a!=b)) ))
        tmp = frac;
    // return 1/m_stretches(c) * dSa_db(a,c) * delta(a,b) * delta(c,d) + frac * (delta(a,c)*delta(b,d) + delta(a,d)*delta(b,c)) * (1-delta(a,b));

    return tmp;
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <bool _comp>
typename std::enable_if<!_comp, T>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    // Compute part with stress tensor involved.
    T frac = 0.0;
    T tmp = 0.0;

    if (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-14)
    {
        // gsDebug<<"Stretches are equal; (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) = "<<abs((m_stretches(a) - m_stretches(b)) / m_stretches(a))<<"\n";
        frac = 1.0 / (2.0 * m_stretches(a) ) * ( dSa_db(b,b) - dSa_db(a,b));
    }
    else
        frac = ( Sa(b)-Sa(a) ) / (math::pow(m_stretches(b),2) - math::pow(m_stretches(a),2));

    GISMO_ENSURE( ( (a < 2) && (b < 2) && (c < 2) && (d < 2) ) , "Index out of range. a="<<a<<", b="<<b<<", c="<<c<<", d="<<d);
    if ( ( (a==b) && (c==d)) )
        tmp = 1/m_stretches(c) * dSa_db(a,c) + 1/(math::pow(m_stretches(a),2) * math::pow(m_stretches(c),2)) * ( math::pow(m_stretches(2),2) * d2Psi_dab(2,2) + 2*dPsi_da(2)*m_stretches(2) );
    else if (( (a==d) && (b==c) && (a!=b) ) || ( ( (a==c) && (b==d) && (a!=b)) ))
        tmp = frac;
    // return 1/m_stretches(c) * dSa_db(a,c) * delta(a,b) * delta(c,d) + frac * (delta(a,c)*delta(b,d) + delta(a,d)*delta(b,c)) * (1-delta(a,b))
                // + delta(a,b)*delta(c,d)*1/(math::pow(m_stretches(a),2) * math::pow(m_stretches(c),2)) * ( math::pow(m_stretches(2),2) * d2Psi_dab(2,2) + 2*dPsi_da(2)*m_stretches(2) );

    return tmp;
}

//--------------------------------------------------------------------------------------------------------------------------------------

// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::dPsi(const index_t a) const
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
//         computeStretch(C);

//         return mu*m_stretches.at(a);
//     }
//     else
//         GISMO_ERROR("Material model not implemented.");
// }

// template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
// T gsMaterialMatrix<dim,T,matId,comp,mat,imp>::d2Psi(const index_t a, const index_t b) const
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
// ---------------------------------------------------------------------------------------------------------------------------------
//                                          Metric Computations
// ---------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computeMetricDeformed() const
{
    computeMetricDeformed_impl<dim>();
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computeMetricDeformed_impl() const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computeMetricDeformed_impl() const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computeMetricUndeformed() const
{
    computeMetricUndeformed_impl<dim>();
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computeMetricUndeformed_impl() const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computeMetricUndeformed_impl() const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getMetric(index_t k, T z) const
{
    this->getMetricDeformed(k,z);
    this->getMetricUndeformed(k,z);

    T ratio = m_Gcov_def.determinant() / m_Gcov_ori.determinant();
    GISMO_ENSURE(ratio > 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_Gcov_ori.determinant());
    m_J0_sq = ratio;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getMetricDeformed(index_t k, T z) const
{
    getMetricDeformed_impl<dim>(k,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getMetricDeformed_impl(index_t k, T z) const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getMetricDeformed_impl(index_t k, T z) const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getMetricUndeformed(index_t k, T z) const
{
    getMetricUndeformed_impl<dim>(k,z);
}

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getMetricUndeformed_impl(index_t k, T z) const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrix<dim,T,matId,comp,mat,imp>::getMetricUndeformed_impl(index_t k, T z) const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrix<dim,T,matId,comp,mat,imp>::evalStretch(const gsMatrix<T> & C) const
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

template <short_t dim, class T, index_t matId, bool comp, enum Material mat, enum Implementation imp >
void gsMaterialMatrix<dim,T,matId,comp,mat,imp>::computeStretch(const gsMatrix<T> & C) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = evalStretch(C);
    m_stretches = result.first;
    m_stretchvec = result.second;
}


} // end namespace
