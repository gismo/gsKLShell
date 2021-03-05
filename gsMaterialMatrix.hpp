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

namespace gismo
{

int delta(const int a, const int b)
{
    return (a==b) ? 1 : 0;
}

int idelta(const int a, const int b)
{
    return (a!=b) ? 1 : 0;
}

// Linear material models
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
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
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
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
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
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
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T>            & mp,
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
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T>            & mp,
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
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T>            & mp,
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
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T>            & mp,
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
// template<class T>
// gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
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

template <class T>
void gsMaterialMatrix<T>::info() const
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

template <class T>
void gsMaterialMatrix<T>::defaultOptions()
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

template <class T>
short_t gsMaterialMatrix<T>::domainDim() const { return 2; }

template <class T>
short_t gsMaterialMatrix<T>::targetDim() const
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

template <class T>
void gsMaterialMatrix<T>::getOptions() const
{
    m_material = m_options.getInt("MaterialLaw");
    m_compressible = m_options.getInt("Compressibility");
    m_integration = m_options.getInt("IntegrationMethod");
}

template <class T>
void gsMaterialMatrix<T>::initialize()
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


template <class T>
void gsMaterialMatrix<T>::computePoints(const gsMatrix<T> & u, bool deformed) const
{
    gsMatrix<T> tmp;

    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    this->computeMetricUndeformed();

    if (m_defpatches->nPieces()!=0)
    {
        m_map_def.flags = m_map.flags;
        m_map_def.points = u;
        static_cast<const gsFunction<T>&>(m_defpatches->piece(0)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
        this->computeMetricDeformed();
    }

    // Compute the thickness
    if (m_material!=1)
    {
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
    if (m_material==14) // check ogden conditions!
    {
        T prod, sum, mu;
        for (index_t c=0; c!=m_parmat.cols(); c++)
        {
            for (index_t r=0; r!=m_parmat.rows(); r++)
                m_parvals.at(r) = m_parmat(r,c);

            mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
            GISMO_ASSERT((m_numPars-2 )% 2 ==0, "Ogden material models must have an even number of parameters (tuples of alpha_i and mu_i). m_numPars = "<< m_numPars);
            int n = (m_numPars-2)/2;
            sum = 0.0;
            for (index_t k=0; k!=n; k++)
            {
                prod = m_parvals.at(2*(k+1))*m_parvals.at(2*(k+1)+1);
                GISMO_ASSERT(prod > 0.0,"Product of coefficients must be positive for all indices");
                sum += prod;
            }
            GISMO_ASSERT((sum-2.*mu)/sum<1e-10,"Sum of products must be equal to 2*mu! sum = "<<sum<<"; 2*mu = "<<2.*mu);
        }
    }
}

template <class T>
void gsMaterialMatrix<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->getOptions();
    if ((m_outputType==1) || (m_outputType==2)) // for matrix and vector
    {
        if (m_integration==0) // DD method
            this->eval_into_NP(u,result);
            // this->eval_into_DD(u,result);
        else if (m_integration==1) // AP method
            this->eval_into_NP(u,result);
            // this->eval_into_AP(u,result);
        else if (m_integration==2) // NP method
            this->eval_into_NP(u,result);
    }
    // This is for density and thickness output
    else if (m_outputType==0)
        this->eval_into_dens(u,result);
    else if (m_outputType==9)
        this->eval_into_stretch(u,result); // midplane stretches
    else if (m_outputType==10)
        this->eval_into_NP(u,result);
    else if (m_outputType==11)
        this->eval_into_stretchdir(u,result); // stretch directions
    else
        GISMO_ERROR("Output type unknown");
}

template <class T>
void gsMaterialMatrix<T>::eval_into_dens(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.flags = NEED_VALUE;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    if (m_material==1)
    {
        result.resize(this->targetDim(), 1);
        GISMO_ASSERT(m_thickValues.size()==m_densities.size(),"Size of vectors of thickness and densities is not equal: " << m_thickValues.size()<<" & "<<m_phis.size());
        for (size_t i = 0; i != m_thickValues.size(); ++i) // loop over laminates
        {
            m_t = m_thickValues[i];
            m_rho = m_densities[i];
            result(0,0) = m_t*m_rho;
        }
        result.replicate(0,u.cols());
    }
    else
    {
        result.resize(this->targetDim(), u.cols());
        m_thickness->eval_into(m_map.values[0], m_Tmat);
        m_density->eval_into(m_map.values[0], m_rhomat);
        for (index_t i = 0; i != u.cols(); ++i) // points
        {
            result(0,i) = m_Tmat(0,i)*m_rhomat(0,i);
        }
    }
}

template <class T>
void gsMaterialMatrix<T>::eval_into_stretch(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    this->computePoints(u);
    result.resize(this->targetDim(), u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;

    if (!m_compressible)
    {
        for (index_t i=0; i!= u.cols(); i++)
        {
            this->getMetric(i,0.0); // on point i, with height 0.0

            gsMatrix<T> C(3,3);
            C.setZero();
            C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
            // C.block(0,0,2,2) = (m_gcov_def.transpose()*m_gcov_def).block(0,0,2,2);
            // gsDebugVar(m_gcov_def.transpose()*m_gcov_def);
            C(2,2) = 1./m_J0_sq;

            res = evalStretch(C);
            result.col(i) = res.first;
        }
    }
    else
    {
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

                GISMO_ASSERT(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
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
                GISMO_ASSERT(it != itmax-1,"Error: Method did not converge, abs(dc33) = "<<abs(dc33)<<" and tolerance = "<<tol<<"\n");
            }
        }
    }
}

template <class T>
void gsMaterialMatrix<T>::eval_into_stretchdir(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    this->computePoints(u);
    result.resize(this->targetDim(), u.cols());
    std::pair<gsVector<T>,gsMatrix<T>> res;

    if (!m_compressible)
    {
        for (index_t i=0; i!= u.cols(); i++)
        {
            this->getMetric(i,0.0); // on point i, with height 0.0

            gsMatrix<T> C(3,3);
            C.setZero();
            C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
            // C.block(0,0,2,2) = (m_gcov_def.transpose()*m_gcov_def).block(0,0,2,2);
            // gsDebugVar(m_gcov_def.transpose()*m_gcov_def);
            C(2,2) = 1./m_J0_sq;

            res = evalStretch(C);
            result.col(i) = res.second.reshape(9,1);
            break;
        }
    }
    else
    {
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

                GISMO_ASSERT(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
                cinv(2,2) = 1.0/c(2,2);

                m_J_sq = m_J0_sq * c(2,2) ;

                S33     = Sij(2,2,c,cinv);
                C3333   = Cijkl3D(2,2,2,2,c,cinv); //  or Cijkl???

                dc33 = -2. * S33 / C3333;
                if (abs(dc33) < tol)
                {
                    res = evalStretch(c);
                    result.col(i) = res.second.reshape(9,1);
                    break;
                }
                GISMO_ASSERT(it != itmax-1,"Error: Method did not converge, abs(dc33) = "<<abs(dc33)<<" and tolerance = "<<tol<<"\n");
            }
        }
    }
}

template <class T>
void gsMaterialMatrix<T>::eval_into_NP(const gsMatrix<T>& u, gsMatrix<T>& result) const
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

    this->computePoints(u);
    if (m_material==0)
    {
        result = multiplyZ(u);
    }
    else if (m_material==1)
        result = eval_Composite(u,m_moment); // thickness integration is done inside
    else
        result = integrateZ(u);
}

/*
    To do:
    - make incompressible and compressible possible
    - for all incompressible material models, we want to change only the formulations of S and C. So they should be functors or someting
*/

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval3D(const index_t i, const gsMatrix<T> & z) const
{
    gsMatrix<T> result;

    if (m_material == 0 ) // svk
        result = eval_Incompressible(i, z);
    else if (m_material == 1) // composite
        gsInfo<<"No eval3D function for composites available. Something went wrong...";
    else if ((m_material >= 2) && (m_compressible)) // NHK
        result = eval_Compressible(i, z);
    else if ((m_material >= 2) && (!m_compressible)) // NHK
        result = eval_Incompressible(i, z);

    // else if ((m_material == 4) && (m_compressible)) // NHK
    //     result = eval_Compressible(i, z);
    // else if ((m_material == 4) && (!m_compressible)) // NHK
    //     result = eval_Incompressible(i, z);


    else if (m_material==-1) // test for integration
    {
        result.resize(this->targetDim(), z.cols());
        result.setZero();
        for (index_t k=0; k!=z.cols(); k++)
            for (index_t r = 0; r!=this->targetDim(); r++)
            {
                result(r,k) = math::pow(z(0,k),r);
            }

    }
    else
        GISMO_ERROR("no function available.");
    // else if ((m_material == 3) && (m_compressible)) // NHK
    //     result = eval_Compressible(i, z);
    // else if ((m_material == 3) && (!m_compressible)) // NHK
    //     result = eval_Incompressible(i, z);
    return result;
}

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval3D(const index_t i) const
{
    gsMatrix<T> z(1,1);
    z.setZero();
    return eval3D(i,z);
}

// multiplies the material matrix by the thickness on all points
// NOTE: this function is a little outdated but it works in its current configuration (note also the implementation of the SvK stress in the Sij function).
// This function needs to be implemented in analytically projected integrals (see Roohbakashan & Sauer)
template <class T>
gsMatrix<T> gsMaterialMatrix<T>::multiplyZ(const gsMatrix<T>& u) const
{
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
template <class T>
gsMatrix<T> gsMaterialMatrix<T>::integrateZ(const gsMatrix<T>& u) const
{
    // Input: points in R2
    // Ouput: results in targetDim

    // Perform integration
    m_numGauss = m_options.getInt("NumGauss");

    gsMatrix<T> result(9,1);
    result.resize(this->targetDim(),u.cols());
    result.setZero();

    m_gauss = gsGaussRule<T>(m_numGauss);
    gsMatrix<T> pts(1,m_numGauss);
    gsMatrix<T> evalPoints;
    // m_points3D.resize(1,m_numGauss);

    gsMatrix<T> quNodes(1,m_numGauss);
    gsVector<T> quWeights(m_numGauss);
    // m_points.conservativeResize(m_points.rows()+1,Eigen::NoChange);

    T res;
    for (index_t j = 0; j != u.cols(); ++j) // for all points
    {
        // Compute values of in-plane points
        // m_points3D.topRows(u.rows()) = u.col(j).replicate(1,m_numGauss);
        // m_points3D.row(0) = u.col(j).replicate(1,m_numGauss);

        // set new integration point
        m_tHalf = m_Tmat(0,j)/2.0;
        m_gauss.mapTo(-m_tHalf,m_tHalf,quNodes,quWeights);
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

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval_Incompressible(const index_t i, const gsMatrix<T>& z) const
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

            // GISMO_ASSERT(m_material >= 10 && m_material < 20, "Only available for stretch-based materials.");

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

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval_Composite(const gsMatrix<T>& u, const index_t moment) const
{
    // Initialize and result
    gsMatrix<T> result(this->targetDim(), u.cols());
    gsMatrix<T,3,3> Dmat, Tmat, Cmat;

    Cmat.setZero();
    Dmat.setZero();
    Tmat.setZero();
    // static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);
    GISMO_ASSERT(m_YoungsModuli.size()==m_PoissonRatios.size(),"Size of vectors of Youngs Moduli and Poisson Ratios is not equal: " << m_YoungsModuli.size()<<" & "<<m_PoissonRatios.size());
    GISMO_ASSERT(m_YoungsModuli.size()==m_ShearModuli.size(),"Size of vectors of Youngs Moduli and Shear Moduli is not equal: " << m_YoungsModuli.size()<<" & "<<m_ShearModuli.size());
    GISMO_ASSERT(m_thickValues.size()==m_phis.size(),"Size of vectors of thickness and angles is not equal: " << m_thickValues.size()<<" & "<<m_phis.size());
    GISMO_ASSERT(m_YoungsModuli.size()==m_thickValues.size(),"Size of vectors of material properties and laminate properties is not equal: " << m_YoungsModuli.size()<<" & "<<m_thickValues.size());
    GISMO_ASSERT(m_YoungsModuli.size()!=0,"No laminates defined");

    // Compute total thickness (sum of entries)
    m_t_tot = std::accumulate(m_thickValues.begin(), m_thickValues.end(), 0.0);

    // compute mid-plane height of total plate
    m_z_mid = m_t_tot / 2.0;

    // now we use t_temp to add the thickness of all plies iteratively
    m_t_temp = 0.0;

    for (size_t i = 0; i != m_phis.size(); ++i) // loop over laminates
    {
        // Compute all quantities
        m_E1 = m_YoungsModuli[i].first;
        m_E2 = m_YoungsModuli[i].second;
        m_G12 = m_ShearModuli[i];
        m_nu12 = m_PoissonRatios[i].first;
        m_nu21 = m_PoissonRatios[i].second;
        m_t = m_thickValues[i];
        m_phi = m_phis[i];

        GISMO_ASSERT(m_nu21*m_E1 == m_nu12*m_E2, "No symmetry in material properties for ply "<<i<<". nu12*E2!=nu21*E1:\n"<<
                "\tnu12 = "<<m_nu12<<"\t E2 = "<<m_E2<<"\t nu12*E2 = "<<m_nu12*m_E2<<"\n"
              <<"\tnu21 = "<<m_nu21<<"\t E1 = "<<m_E1<<"\t nu21*E1 = "<<m_nu21*m_E1);

        // Fill material matrix
        Dmat(0,0) = m_E1 / (1-m_nu12*m_nu21);
        Dmat(1,1) = m_E2 / (1-m_nu12*m_nu21);
        Dmat(2,2) = m_G12;
        Dmat(0,1) = m_nu21*m_E1 / (1-m_nu12*m_nu21);
        Dmat(1,0) = m_nu12*m_E2 / (1-m_nu12*m_nu21);
        Dmat(2,0) = Dmat(0,2) = Dmat(2,1) = Dmat(1,2) = 0.0;

        // Make transformation matrix
        Tmat(0,0) = Tmat(1,1) = math::pow(math::cos(m_phi),2);
        Tmat(0,1) = Tmat(1,0) = math::pow(math::sin(m_phi),2);
        Tmat(2,0) = Tmat(0,2) = Tmat(2,1) = Tmat(1,2) = math::sin(m_phi) * math::cos(m_phi);
        Tmat(2,0) *= -2.0;
        Tmat(2,1) *= 2.0;
        Tmat(1,2) *= -1.0;
        Tmat(2,2) = math::pow(math::cos(m_phi),2) - math::pow(math::sin(m_phi),2);

        // Compute laminate stiffness matrix
        Dmat = Tmat.transpose() * Dmat * Tmat;

        m_z = math::abs(m_z_mid - (m_t/2.0 + m_t_temp) ); // distance from mid-plane of plate

        switch (moment)
        {
            case 0:
                Cmat += Dmat * m_t; // A
                break;
            case 1:
                Cmat += Dmat * m_t*m_z; // B
                break;
            case 2:
                Cmat += Dmat * ( m_t*math::pow(m_z,2) + math::pow(m_t,3)/12.0 ); // D
                break;
        }
        // Note: we put everything in the first column, to replicate further after the loop;
        // since all cols are equal up to now
        // (ASSUMPTION THAT THE MATERIAL PARAMETERS ARE CONSTANT)

        m_t_temp += m_t;
    }

    GISMO_ASSERT(m_t_tot==m_t_temp,"Total thickness after loop is wrong. t_temp = "<<m_t_temp<<" and sum(thickness) = "<<m_t_tot);

    // m_map.points = m_map_def.points = u.topRows(2);
    // static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    // static_cast<const gsFunction<T>&>(m_defpatches->piece(m_pIndex)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    // TRANSFORMATION OF THE MATRIX FROM LOCAL CARTESIAN TO CURVILINEAR
    // Compute covariant bases in deformed and undeformed configuration

    m_covBasis.resize(3,3);
    m_conBasis.resize(3,3);
    for (index_t k = 0; k != u.cols(); ++k)
    {
        m_covBasis.leftCols(2) = m_map.jacobian(k);
        m_covBasis.col(2)      = m_map.normal(k).normalized();
        m_covMetric = m_covBasis.transpose() * m_covBasis;

        m_conMetric = m_covMetric.template block<3,3>(0,0).inverse();

        m_conBasis.col(0) = m_conMetric(0,0)*m_covBasis.col(0)+m_conMetric(0,1)*m_covBasis.col(1)+m_conMetric(0,2)*m_covBasis.col(2);
        m_conBasis.col(1) = m_conMetric(1,0)*m_covBasis.col(0)+m_conMetric(1,1)*m_covBasis.col(1)+m_conMetric(1,2)*m_covBasis.col(2);
        // conBasis.col(2) = m_conMetric(2,0)*m_covBasis.col(0)+m_conMetric(2,1)*m_covBasis.col(1)+m_conMetric(2,2)*m_covBasis.col(2);

        m_e1 = m_covBasis.col(0); m_e1.normalize();
        m_e2 = m_conBasis.col(1); m_e2.normalize();
        // e3 = normal;

        m_ac1 = m_conBasis.col(0);
        m_ac2 = m_conBasis.col(1);

        m_a1 = m_covBasis.col(0);
        m_a2 = m_covBasis.col(1);

        Tmat.setZero();
        // Covariant to local cartesian
        Tmat(0,0) = (m_e1.dot(m_a1))*(m_a1.dot(m_e1));
        Tmat(0,1) = (m_e1.dot(m_a2))*(m_a2.dot(m_e2));
        Tmat(0,2) = 2*(m_e1.dot(m_a1))*(m_a2.dot(m_e1));
        // Row 1
        Tmat(1,0) = (m_e2.dot(m_a1))*(m_a1.dot(m_e2));
        Tmat(1,1) = (m_e2.dot(m_a2))*(m_a2.dot(m_e2));
        Tmat(1,2) = 2*(m_e2.dot(m_a1))*(m_a2.dot(m_e2));
        // Row 2
        Tmat(2,0) = (m_e1.dot(m_a1))*(m_a1.dot(m_e2));
        Tmat(2,1) = (m_e1.dot(m_a2))*(m_a2.dot(m_e2));
        Tmat(2,2) = (m_e1.dot(m_a1))*(m_a2.dot(m_e2)) + (m_e1.dot(m_a2))*(m_a1.dot(m_e2));

        Tmat = Tmat.template block<3,3>(0,0).inverse(); // !!!!

        Cmat = Tmat * Cmat;

        Tmat.setZero();
        // Contravariant to local cartesian
        Tmat(0,0) = (m_e1.dot(m_ac1))*(m_ac1.dot(m_e1));
        Tmat(0,1) = (m_e1.dot(m_ac2))*(m_ac2.dot(m_e2));
        Tmat(0,2) = 2*(m_e1.dot(m_ac1))*(m_ac2.dot(m_e1));
        // Row 1
        Tmat(1,0) = (m_e2.dot(m_ac1))*(m_ac1.dot(m_e2));
        Tmat(1,1) = (m_e2.dot(m_ac2))*(m_ac2.dot(m_e2));
        Tmat(1,2) = 2*(m_e2.dot(m_ac1))*(m_ac2.dot(m_e2));
        // Row 2
        Tmat(2,0) = (m_e1.dot(m_ac1))*(m_ac1.dot(m_e2));
        Tmat(2,1) = (m_e1.dot(m_ac2))*(m_ac2.dot(m_e2));
        Tmat(2,2) = (m_e1.dot(m_ac1))*(m_ac2.dot(m_e2)) + (m_e1.dot(m_ac2))*(m_ac1.dot(m_e2));

        Cmat = Cmat * Tmat;

        if (m_outputType==2)
            result.reshapeCol(k,3,3) = Cmat;
        else if (m_outputType==1)
        {
            gsMatrix<T> Eij(2,2);
            // computeMetric(k,0.0,true,true); // height is set to 0, but we do not need m_G

            if (m_moment==0)
                Eij = 0.5*(m_Acov_def - m_Acov_ori);
            else if (m_moment==2)
                Eij = (m_Bcov_def - m_Bcov_ori);
            else
            {
                gsDebug<<"Warning: no material model known in simplification";
            }
            Eij(0,1) += Eij(1,0);
            std::swap(Eij(1,0), Eij(1,1));
            Eij.resize(4,1);
            Eij.conservativeResize(3,1);
            result.col(k) = Cmat * Eij;
        }
        else
            GISMO_ERROR("no vector or matrix produced");
    }

    return result;
};

/*
    Available class members:
        - m_parvals
        - m_metric
        - m_metric_def
        - m_J0
        - m_J
*/
template<class T>
T gsMaterialMatrix<T>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ASSERT( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    T tmp = 0.0;
    // --------------------------
    // Explicit implementations
    // --------------------------
        if (m_material==0) // svk
        {
            // --------------------------
            // Saint Venant Kirchhoff
            // --------------------------
            T lambda, mu, Cconstant;

            mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
            GISMO_ASSERT((1.-2.*m_parvals.at(1)) != 0, "Division by zero in construction of SvK material parameters! (1.-2.*m_parvals.at(1)) = "<<(1.-2.*m_parvals.at(1))<<"; m_parvals.at(1) = "<<m_parvals.at(1));
            lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1))) ;
            Cconstant = 2*lambda*mu/(lambda+2*mu);

            tmp = Cconstant*m_Acon_ori(i,j)*m_Acon_ori(k,l) + mu*(m_Acon_ori(i,k)*m_Acon_ori(j,l) + m_Acon_ori(i,l)*m_Acon_ori(j,k));
        }
        else if (m_material==1)
            GISMO_ERROR("Compressible material matrix  requested, but not needed. How?");
        else if (m_material==2)
        {
            // --------------------------
            // Neo-Hookean
            // --------------------------
            T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
            tmp = mu*1./m_J0_sq*(2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k));
        }

        else if (m_material==3)
        {
            // --------------------------
            // Mooney-Rivlin
            // Parameter 3 is the ratio between c1 and c2.; c1 = m_parvals.at(2)*c2
            // --------------------------
            GISMO_ASSERT(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
            T traceCt =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                        m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                        m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                        m_Gcov_def(1,1)*m_Gcon_ori(1,1);

            T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));
            T c2 = mu/(m_parvals.at(2) + 1);
            T c1 = m_parvals.at(2)*c2;

            T Gabcd = - 1./2. * ( m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k) );

            tmp = (c1 + c2 * traceCt) *1./m_J0_sq*(2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k))// correct
                    - 2. * c2 / m_J0_sq * ( m_Gcon_ori(i,j) * m_Gcon_def(k,l) + m_Gcon_def(i,j)*m_Gcon_ori(k,l)) // correct
                    + 2. * c2 * m_J0_sq * ( Gabcd + m_Gcon_def(i,j)*m_Gcon_def(k,l) ); // Roohbakhshan
        }
        else if (m_material==4)
                GISMO_ERROR("Material model 4 is not invariant-based! Use 14 instead...");
        else if (m_material==5)
                GISMO_ERROR("Material model 5 is only for compressible materials...");

    // --------------------------
    // Stretch-based implementations
    // --------------------------
        else if ((m_material >= 10) && (m_material < 20))
        {
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
        }
    // --------------------------
    // General implementations
    // --------------------------
        else if ((m_material >= 20) && (m_material < 30))
        {
            tmp = 4.0 * d2Psi(i,j,k,l) + 4.0 * d2Psi(2,2,2,2)*math::pow(m_J0_sq,-2.0)*m_Gcon_def(i,j)*m_Gcon_def(k,l)
                    - 4.0/ m_J0_sq  * ( d2Psi(2,2,i,j)*m_Gcon_def(k,l) + d2Psi(2,2,k,l)*m_Gcon_def(i,j) )
                    + 2.0 * dPsi(2,2) / m_J0_sq * (2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k));
        }
    // --------------------------
    // Error
    // --------------------------
        else
            GISMO_ERROR("Material model unknown (model = "<<m_material<<"). Use gsMaterialMatrix<T>::info() to see the options.");
    return tmp;
}

// Condensation of the 3D tensor for compressible materials
template<class T>
T gsMaterialMatrix<T>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ASSERT(c.cols()==c.rows(),"Matrix c must be square");
    GISMO_ASSERT(c.cols()==3,"Matrix c must be 3x3");
    GISMO_ASSERT(cinv.cols()==cinv.rows(),"Matrix cinv must be square");
    GISMO_ASSERT(cinv.cols()==3,"Matrix cinv must be 3x3");
    GISMO_ASSERT( ( (i <2) && (j <2) && (k <2) && (l <2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(m_compressible,"Material model is not compressible?");

    T tmp = 0.0;
    if ((m_material >= 10) && (m_material < 20))
    {
        // static condensation is done before the projection
        computeStretch(c);
        tmp = Cijkl3D(i,j,k,l,c,cinv);
    }
    else
        tmp = Cijkl3D(i,j,k,l,c,cinv) - ( Cijkl3D(i,j,2,2,c,cinv) * Cijkl3D(2,2,k,l,c,cinv) ) / Cijkl3D(2,2,2,2,c,cinv);

    return tmp;
}

// 3D tensor for compressible materials
template<class T>
T gsMaterialMatrix<T>::Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ASSERT( ( (i <3) && (j <3) && (k <3) && (l <3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);

    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    GISMO_ASSERT(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");

    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    // T d2I_1 = 0;

    T I_2 = c(2,2) * traceCt + m_J0_sq;
    T d2I_2 = idelta(i,2)*idelta(j,2)*idelta(k,2)*idelta(l,2)*( m_J0_sq*( cinv(i,j)*cinv(k,l) + dCinv ) )
            + delta(i,2)*delta(j,2)*idelta(k,2)*idelta(l,2)*dI_1(k,l)
            + idelta(i,2)*idelta(j,2)*delta(k,2)*delta(l,2)*dI_1(i,j);
            // + delta(i,2)*delta(j,2)*delta(k,2)*delta(l,2)*0;
    // --------------------------
    // Explicit implementations
    // --------------------------
        if (m_material==0 || m_material==1) // svk
            gsWarn<<"Compressible material matrix requested, but not needed. How?";
        else if (m_material==2)
            tmp = 1.0 / 9.0 * mu * math::pow( m_J_sq , -1.0/3.0 ) * ( 2.0 * I_1 * ( cinv(i,j)*cinv(k,l) - 3.0 * dCinv )
                                - 6.0 *( m_Gcon_ori(i,j)*cinv(k,l) + cinv(i,j)*m_Gcon_ori(k,l) ) )
                    + K * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1)*dCinv );
        else if (m_material==3)
        {
            GISMO_ASSERT(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
            T c2 = mu/(m_parvals.at(2) + 1);
            T c1 = m_parvals.at(2)*c2;
            tmp =     1.0/9.0 * c1 * math::pow(m_J_sq, -1.0/3.0) *  ( 2.0*I_1*cinv(i,j)*cinv(k,l) - 6.0*I_1*dCinv
                                                                    - 6.0*dI_1(i,j)*cinv(k,l)     - 6.0*cinv(i,j)*dI_1(k,l) ) // + 9*d2I_1 = 0
                    + 1.0/9.0 * c2 * math::pow(m_J_sq, -2.0/3.0) *  ( 8.0*I_2*cinv(i,j)*cinv(k,l) - 12.0*I_2*dCinv
                                                                        - 12.0*dI_2(i,j,c,cinv)*cinv(k,l)- 12.0*cinv(i,j)*dI_2(k,l,c,cinv)
                                                                        + 18.0*d2I_2 )
                    + K * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1)*dCinv );
        }
        else if (m_material==4)
                GISMO_ERROR("Material model 4 is not invariant-based! Use 14 instead...");
        else if (m_material==5)
            tmp =  - 2.0 * mu * dCinv + lambda * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1)*dCinv );

    // --------------------------
    // Stretch-based implementations
    // --------------------------
        else if ((m_material >= 10) && (m_material < 20))
        {

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
        }
    // --------------------------
    // General implementations
    // --------------------------
        else if ((m_material >= 20) && (m_material < 30))
            tmp = 4.0 * d2Psi(i,j,k,l,c,cinv);

    // --------------------------
    // Error
    // --------------------------
        else
            GISMO_ERROR("Material model unknown (model = "<<m_material<<"). Use gsMaterialMatrix<T>::info() to see the options.");
    return tmp;
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
// template<class T>
// T gsMaterialMatrix<T>::Sij(const index_t i, const index_t j) const { Sij(i,j,NULL,NULL); }

template<class T>
T gsMaterialMatrix<T>::Sij(const index_t i, const index_t j) const
{
    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2.*(1. + m_parvals.at(1)));

    // --------------------------
    // Explicit implementations
    // --------------------------
        if (m_material==0)
        {
            gsMatrix<T> stress;
            // --------------------------
            // Saint Venant Kirchhoff
            // --------------------------
            if (m_moment==0)
            {
                GISMO_ASSERT( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
                stress = 0.5*(m_Acov_def - m_Acov_ori);
            }
            else if (m_moment==2)
            {
                GISMO_ASSERT( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
                stress = (m_Bcov_ori - m_Bcov_def);
            }
            else
            {
                GISMO_ERROR("Warning: no material model known in simplification, m_moment="<<m_moment);
            }

            // ALTERNATIVE
            // stress = 0.5 * (m_Gcov_def_L - m_Gcov_ori_L);

            tmp =     Cijkl(i,j,0,0) * stress(0,0) + Cijkl(i,j,0,1) * stress(0,1)
                    + Cijkl(i,j,1,0) * stress(1,0) + Cijkl(i,j,1,1) * stress(1,1);
        }
        else if (m_material==1)
        {
            // --------------------------
            // Composite
            // --------------------------
            GISMO_ERROR("Incompressible material stress tensor requested, but not needed. How?");
        }
        else if (m_material==2)
        {
            // --------------------------
            // Neo-Hoookean
            // --------------------------
            // gsDebugVar(mu);
            // gsDebugVar(mu * (m_Gcon_ori(i,j) - 1./math::pow(m_J0,2.) * m_Gcon_def(i,j) ));
            tmp = mu * (m_Gcon_ori(i,j) - 1./m_J0_sq * m_Gcon_def(i,j) );
        }
        else if (m_material==3)
        {
            // --------------------------
            // Mooney-Rivlin
            // Parameter 3 is the ratio between c1 and c2.
            // --------------------------
            GISMO_ASSERT(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
            T traceCt =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                        m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                        m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                        m_Gcov_def(1,1)*m_Gcon_ori(1,1);

            T c2 = mu/(m_parvals.at(2)+1);
            T c1 = m_parvals.at(2)*c2;

            tmp = c1 * ( m_Gcon_ori(i,j) - 1/m_J0_sq * m_Gcon_def(i,j) )
                    + c2 / m_J0_sq * (m_Gcon_ori(i,j) - traceCt * m_Gcon_def(i,j) ) + c2 * m_J0_sq * m_Gcon_def(i,j);
        }
    // --------------------------
    // Stretch-based implementations
    // --------------------------
        else if ( (m_material>=10) && (m_material<20) )
        {
            gsMatrix<T> C(3,3);
            C.setZero();
            C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
            C(2,2) = 1./m_J0_sq;

            computeStretch(C);

            tmp = 0.0;
            for (index_t a = 0; a != 2; a++)
            {
                tmp += Sa(a)*(
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )
                            );
            }
        }
    // --------------------------
    // General implementations
    // --------------------------
        else if ( (m_material>=20) && (m_material<30) )
        {
            // --------------------------
            // Generalized
            // --------------------------
            tmp = 2.0 * dPsi(i,j) - 2.0 * dPsi(2,2) * math::pow(m_J0_sq,-1.0)*m_Gcon_def(i,j);
        }
    // --------------------------
    // Error
    // --------------------------
        else
            GISMO_ERROR("Material model unknown. Use gsMaterialMatrix<T>::info() to see the options.");
    return tmp;
}

template<class T>
T gsMaterialMatrix<T>::Sij(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ASSERT(c.cols()==c.rows(),"C must be square");
    GISMO_ASSERT(c.cols()==3,"C must be 3x3");
    GISMO_ASSERT(cinv.cols()==cinv.rows(),"Cinv must be square");
    GISMO_ASSERT(cinv.cols()==3,"Cinv must be 3x3");

    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    GISMO_ASSERT(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");

    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    T I_2   = c(2,2) * traceCt + m_J0_sq;

    if (m_material==0 || m_material==1)
        GISMO_ERROR("Incompressible material stress tensor requested, but not needed. How?");
    else if (m_material==2)
        tmp =  mu * math::pow( m_J_sq , -1.0/3.0 ) * ( m_Gcon_ori(i,j) - 1.0/3.0 * I_1 * cinv(i,j) )
                + K * 0.5 * ( m_J_sq - 1.0 ) * cinv(i,j);
    else if (m_material==3)
    {
        GISMO_ASSERT(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
        T c2 = mu/(m_parvals.at(2)+1);
        T c1 = m_parvals.at(2)*c2;

        tmp =     c1 * math::pow( m_J_sq , -1.0/3.0 ) * ( m_Gcon_ori(i,j) - 1.0/3.0 * I_1 * cinv(i,j) )
                + c2 * math::pow( m_J_sq , -2.0/3.0 ) * ( dI_2(i,j,c,cinv)- 2.0/3.0 * I_2 * cinv(i,j) )
                + K * 0.5 * ( m_J_sq - 1.0 ) * cinv(i,j);
    }
    else if (m_material==5)
        tmp = mu * m_Gcon_ori(i,j) - mu * cinv(i,j) + lambda / 2.0 * ( m_J_sq - 1 ) * cinv(i,j);

    // --------------------------
    // Stretch-based implementations
    // --------------------------
        else if ( (m_material>=10) && (m_material<20) )
        {
            computeStretch(c);

            tmp = 0.0;

            for (index_t a = 0; a != 3; a++)
            {
                tmp += Sa(a)*(
                            ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(a)) )
                            );
            }
        }

    // --------------------------
    // General implementations
    // --------------------------
        else if ((m_material >= 20) && (m_material < 30))
        {
            tmp = 2.0 * dPsi(i,j,c,cinv);
        }
    return tmp;
}

template<class T>
T gsMaterialMatrix<T>::Sii(const index_t i) const // principle stresses
{
    gsMatrix<T> C(3,3);
    C.setZero();
    C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
    C(2,2) = 1./m_J0_sq;
    return Sii(i,C);
}

template<class T>
T gsMaterialMatrix<T>::Sii(const index_t i, const gsMatrix<T> & c) const
{
    computeStretch(c);
    return Sa(i);
}

template<class T>
gsMatrix<T> gsMaterialMatrix<T>::eval_Compressible(const index_t i, const gsMatrix<T>& z) const
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
        cinv.setZero();
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

            GISMO_ASSERT(c(2,2)>= 0,"ERROR in iteration "<<it<<"; c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
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
                        GISMO_ASSERT(m_material >= 10 && m_material < 20, "Only available for stretch-based materials.");

                        result(0,j) = Sii(0,c); // S11
                        result(1,j) = Sii(1,c); // S22
                    }
                else
                    GISMO_ERROR("no vector or matrix produced");

                    break;
            }
            GISMO_ASSERT(it != itmax-1,"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n");
        }
    }
    return result;
}

// ---------------------------------------------------------------------------------------------------------------------------------
//                                          INCOMPRESSIBLE
// ---------------------------------------------------------------------------------------------------------------------------------
template<class T>
T gsMaterialMatrix<T>::dPsi(const index_t i, const index_t j) const
{
    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));

    GISMO_ASSERT( ( (i < 3) && (j < 3) ) , "Index out of range. i="<<i<<", j="<<j);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material == 22)
        tmp = 0.5 * mu * m_Gcon_ori(i,j);
    else if (m_material==23)
    {
        T c2 = mu/(m_parvals.at(2)+1);
        T c1 = m_parvals.at(2)*c2;
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
    }
    else if (m_material==24)
            GISMO_ERROR("Material model 24 is not invariant-based! Use 14 instead...");
    else if (m_material==25)
                GISMO_ERROR("Material model 25 is only for compressible materials...");
    else
        GISMO_ERROR("Material model not implemented.");

    return tmp;
}

template<class T>
T gsMaterialMatrix<T>::d2Psi(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));

    GISMO_ASSERT( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material == 22)
        tmp = 0.0;
    else if (m_material==23)
    {
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
    }
    else if (m_material==24)
            GISMO_ERROR("Material model 24 is not invariant-based! Use 14 instead...");
    else if (m_material==25)
            GISMO_ERROR("Material model 25 is only for compressible materials...");
    else
        GISMO_ERROR("Material model not implemented. material = "<<m_material);

    return tmp;
}
// ---------------------------------------------------------------------------------------------------------------------------------
//                                          COMPRESSIBLE
// ---------------------------------------------------------------------------------------------------------------------------------
template<class T>
T gsMaterialMatrix<T>::dI_1(const index_t i, const index_t j) const
{
    return m_Gcon_ori(i,j);
}

template<class T>
T gsMaterialMatrix<T>::dI_2(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    return idelta(i,2)*idelta(j,2)*( c(2,2)*dI_1(i,j) + m_J0_sq*cinv(i,j) ) + delta(i,2)*delta(j,2)*traceCt;
}

template<class T>
T gsMaterialMatrix<T>::dPsi(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));
    GISMO_ASSERT(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");

    //  choose compressibility function (and parameter)
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
    T dpsi_vol = K * 0.25 * (m_J_sq - 1.0) * cinv(i,j);

    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);
    T I_2 = c(2,2) * traceCt + m_J0_sq;

    GISMO_ASSERT( ( (i < 3) && (j < 3) ) , "Index out of range. i="<<i<<", j="<<j);
    GISMO_ASSERT(m_compressible,"Material model is not compressible?");

    if (m_material == 22)
        tmp = mu/2.0 * math::pow(m_J_sq,-1./3.) * ( - 1.0/3.0 * I_1 * cinv(i,j) + dI_1(i,j) );
    else if (m_material==23)
    {
        T c2= mu/(m_parvals.at(2)+1);
        T c1= m_parvals.at(2)*c2;
        // c1 = 0;
        tmp = c1/2.0 * math::pow(m_J_sq,-1./3.) * ( - 1.0/3.0 * I_1 * cinv(i,j) + dI_1(i,j) )
            + c2/2.0 * math::pow(m_J_sq,-2./3.) * ( - 2.0/3.0 * I_2 * cinv(i,j) + dI_2(i,j,c,cinv) );
    }
    else if (m_material==25)
        tmp = mu / 2.0 * dI_1(i,j) - mu / 2.0 * cinv(i,j) + lambda / 4.0 * ( m_J_sq - 1 ) * cinv(i,j) - dpsi_vol; // subtract dpsi_vol since it will be added later
    else
        GISMO_ERROR("Material model not implemented (Cijkl).");

    return tmp + dpsi_vol;
}

template<class T>
T gsMaterialMatrix<T>::d2Psi(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ASSERT( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(m_compressible,"Material model is not compressible?");

    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    GISMO_ASSERT(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");

    T dCinv = - 1./2.*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) );

    //  choose compressibility function (and parameter)
    T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
    T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

    T d2psi_vol = K * 0.25 * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1.0)*dCinv );

    T traceCt = m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1);
    T I_1   = traceCt + c(2,2);

    T I_2 = c(2,2) * traceCt + m_J0_sq;
    T d2I_2 = idelta(i,2)*idelta(j,2)*idelta(k,2)*idelta(l,2)*( m_J0_sq*( cinv(i,j)*cinv(k,l) + dCinv ) )
            + delta(i,2)*delta(j,2)*idelta(k,2)*idelta(l,2)*dI_1(k,l)
            + idelta(i,2)*idelta(j,2)*delta(k,2)*delta(l,2)*dI_1(i,j);
            // + delta(i,2)*delta(j,2)*delta(k,2)*delta(l,2)*0;

    if (m_material == 22)
        tmp = 1.0/9.0 * mu / 2.0 * math::pow(m_J_sq, -1.0/3.0) *  ( I_1*cinv(i,j)*cinv(k,l)
                                                                - 3.0*dI_1(i,j)*cinv(k,l)        - 3.0*cinv(i,j)*dI_1(k,l)
                                                                - 3.0*I_1*dCinv  ); // + 9*d2I_1 = 0
    else if (m_material==23)
    {
        T c2 = mu/(m_parvals.at(2)+1);
        T c1 = m_parvals.at(2)*c2;
        // c1 = 0;
        tmp =
              1.0/9.0 * c1 / 2.0 * math::pow(m_J_sq, -1.0/3.0) *  ( I_1*cinv(i,j)*cinv(k,l)
                                                                - 3.0*dI_1(i,j)*cinv(k,l)       - 3.0*cinv(i,j)*dI_1(k,l)
                                                                - 3.0*I_1*dCinv ) // + 9*d2I_1 = 0
            + 1.0/9.0 * c2 / 2.0 * math::pow(m_J_sq, -2.0/3.0) *  ( 4.0*I_2*cinv(i,j)*cinv(k,l) - 6.0*I_2*dCinv
                                                                - 6.0*dI_2(i,j,c,cinv)*cinv(k,l)- 6.0*cinv(i,j)*dI_2(k,l,c,cinv)
                                                                + 9.0*d2I_2 );
    }
    else if (m_material==25)
        tmp = - mu / 2.0 * dCinv + lambda / 4.0 * ( m_J_sq*cinv(i,j)*cinv(k,l) + (m_J_sq-1.0)*dCinv ) - d2psi_vol;  // subtract d2psi_vol since it will be added later

    else
        GISMO_ERROR("Material model not implemented.");

    return tmp + d2psi_vol;
}
// ---------------------------------------------------------------------------------------------------------------------------------
//                                          STRETCHES
// ---------------------------------------------------------------------------------------------------------------------------------

template<class T>
T gsMaterialMatrix<T>::dPsi_da(const index_t a) const
{
    GISMO_ASSERT( a < 3 , "Index out of range. a="<<a);
    T tmp = 0.0;
    T mu = m_parvals.at(0) / (2 * (1 + m_parvals.at(1)));
    T I_1   = m_stretches(0)*m_stretches(0) + m_stretches(1)*m_stretches(1) + m_stretches(2)*m_stretches(2);
    T dI_1a = 2*m_stretches(a);
    T I_2   = math::pow(m_stretches(0),2.)*math::pow(m_stretches(1),2.)
            + math::pow(m_stretches(1),2.)*math::pow(m_stretches(2),2.)
            + math::pow(m_stretches(0),2.)*math::pow(m_stretches(2),2.);
    T dI_2a  = 2*m_stretches(a)*( I_1 - math::pow(m_stretches(a),2.0) );

    if (!m_compressible) // incompressible
    {
            if (m_material==12 || m_material==2 || m_material==22)// second and third case for plotting principal stress
            {
                tmp  = mu/2 * dI_1a;
            }
            else if (m_material==13 || m_material==3 || m_material==23)// second and third case for plotting principal stress
            {
                T c2 = mu/(m_parvals.at(2)+1);
                T c1 = m_parvals.at(2)*c2;
                tmp  = c1/2.0*dI_1a + c2/2.0*dI_2a;
            }
            else if (m_material==14)
            {
                int n = (m_numPars-2)/2;
                T alpha_i, mu_i;
                for (index_t k=0; k!=n; k++)
                {
                    alpha_i = m_parvals.at(2*(k+1)+1);
                    mu_i = m_parvals.at(2*(k+1));
                    tmp += mu_i*math::pow(m_stretches(a),alpha_i-1);
                }
            }
            else if (m_material==15)
                GISMO_ERROR("Material model 15 is only for compressible materials...");
            else
                GISMO_ERROR("Material model "<<m_material<<" not available...");
    }
    else // compressible
    {
            GISMO_ASSERT(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
            T beta  = -2.0;
            //  choose compressibility function (and parameter)
            T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
            T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));

            T dpsi_vol = K / (m_stretches(a)*beta) * (1.0 - math::pow(m_J_sq,-beta/2.0));

            if (m_material==12)
                tmp = mu/2.0 * math::pow(m_J_sq,-1./3.) * ( -2./3. *  I_1 / m_stretches(a) + dI_1a );
            else if (m_material==13)
            {
                T c2= mu/(m_parvals.at(2)+1);
                T c1= m_parvals.at(2)*c2;
                tmp = c1/2.0 * math::pow(m_J_sq,-1./3.) * ( -2./3. *  I_1 / m_stretches(a) + dI_1a )
                    + c2/2.0 * math::pow(m_J_sq,-2./3.) * ( -4./3. *  I_2 / m_stretches(a) + dI_2a );
            }
            else if (m_material==14)
            {
                int n = (m_numPars-2)/2;
                T alpha_i, mu_i, Lambda;
                for (index_t k=0; k!=n; k++)
                {
                    alpha_i = m_parvals.at(2*(k+1)+1);
                    mu_i = m_parvals.at(2*(k+1));
                    Lambda = math::pow(m_stretches(0),alpha_i) + math::pow(m_stretches(1),alpha_i) + math::pow(m_stretches(2),alpha_i);
                    tmp += mu_i * math::pow(m_J_sq,-alpha_i/6.0) * ( math::pow(m_stretches(a),alpha_i-1) - 1./3. * 1./m_stretches(a) * Lambda );
                }
            }
            else if (m_material==15)
                tmp = mu/2.0 * dI_1a - mu / m_stretches(a) + lambda / (m_stretches(a)*2) * (m_J_sq-1.0) - dpsi_vol;

            else
                GISMO_ERROR("not available...");

            tmp += dpsi_vol;
    }


    return tmp;
}

template<class T>
T gsMaterialMatrix<T>::d2Psi_dab(const index_t a, const index_t b) const
{
    GISMO_ASSERT( ( (a < 3) && (b < 3) ) , "Index out of range. a="<<a<<", b="<<b);
    T tmp = 0.0;
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

    // T d2I_2;
    // if (a!=b)
    //     d2I_2 = 4.0*m_stretches(a)*m_stretches(b);
    // else
    //     d2I_2 = 2.0*(I_1 - m_stretches(a)*m_stretches(a));
    T d2I_2 = idelta(a,b)*4.0*m_stretches(a)*m_stretches(b) + delta(a,b)*2.0*(I_1 - m_stretches(a)*m_stretches(a));

    if (!m_compressible)
    {
        if (m_material==12)
        {
            tmp  = mu/2 * d2I_1;
        }
        else if (m_material==13)
        {
            GISMO_ASSERT(m_numPars==3,"Mooney-Rivlin model needs to be a 3 parameter model");
            T c2 = mu/(m_parvals.at(2)+1);
            T c1 = m_parvals.at(2)*c2;
            tmp  = c1/2.0 * d2I_1 + c2/2.0 * d2I_2;
        }
        else if (m_material==14)
        {
            int n = (m_numPars-2)/2;
            T alpha_i, mu_i;
            for (index_t k=0; k!=n; k++)
            {
                alpha_i = m_parvals.at(2*(k+1)+1);
                mu_i = m_parvals.at(2*(k+1));
                tmp += mu_i*math::pow(m_stretches(a),alpha_i-2)*(alpha_i-1)*delta(a,b);
            }
        }
        else if (m_material==15)
            GISMO_ERROR("Material model 15 is only for compressible materials...");
        else
            GISMO_ERROR("not available...");
    }
    else
    {
        GISMO_ASSERT(3 - 6 * m_parvals.at(1) != 0, "Bulk modulus is infinity for compressible material model. Try to use incompressible models.");
        m_J_sq = math::pow(m_stretches(0)*m_stretches(1)*m_stretches(2),2.0);

        T beta  = -2.0;
        T K  = m_parvals.at(0) / ( 3 - 6 * m_parvals.at(1));
        T lambda = m_parvals.at(0) * m_parvals.at(1) / ( (1. + m_parvals.at(1))*(1.-2.*m_parvals.at(1)));
        T d2psi_vol = K / (beta*m_stretches(a)*m_stretches(b)) * ( beta*math::pow(m_J_sq,-beta/2.0) + delta(a,b) * (math::pow(m_J_sq,-beta/2.0) - 1.0) );

        if (m_material==12)
        {
            tmp = mu/2.0 * math::pow(m_J_sq,-1./3.) *   (
                                                            -2./3. * 1. / m_stretches(b) * ( -2./3. * I_1 / m_stretches(a) + dI_1a )
                                                            -2./3. * 1. / m_stretches(a) * dI_1b
                                                            +d2I_1
                                                            +2./3. * delta(a,b) * I_1 / (m_stretches(a)*m_stretches(a))
                                                    );
        }
        else if (m_material==13)
        {
            T c2 = mu/(m_parvals.at(2)+1);
            T c1 = m_parvals.at(2)*c2;
            tmp =
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
                                                    );
        }
        else if (m_material==14)
        {
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
        }
        else if (m_material==15)
            tmp = mu/2.0 * d2I_1 + mu * delta(a,b) / ( m_stretches(a) * m_stretches(b) ) + lambda / (2*m_stretches(a)*m_stretches(b)) * ( 2*m_J_sq - delta(a,b) * (m_J_sq - 1.0) ) - d2psi_vol;
        else
            GISMO_ERROR("not available...");

        tmp += d2psi_vol;
    }
    return tmp;
}

template<class T>
T gsMaterialMatrix<T>::dJ_da(const index_t a) const
{
    return 1.0/m_stretches(a);
}

template<class T>
T gsMaterialMatrix<T>::d2J_dab(const index_t a, const index_t b) const
{
    if (a==b)
        return 0.0;
    else
        return 1.0  / ( m_stretches(a) * m_stretches(b) );
    // return ( 1.0 - delta(a,b) ) / ( m_stretches(a) * m_stretches(b) );
}

template<class T>
T gsMaterialMatrix<T>::p() const
{
    return m_stretches(2) * dPsi_da(2);
}

template<class T>
T gsMaterialMatrix<T>::dp_da(const index_t a) const
{
    if (a==2)
        return m_stretches(2) * d2Psi_dab(2,a) + dPsi_da(2);
    else
        return m_stretches(2) * d2Psi_dab(2,a);

    // return m_stretches(2) * d2Psi_dab(2,a) + delta(a,2) * dPsi_da(2);
}

template<class T>
T gsMaterialMatrix<T>::Sa(const index_t a) const
{
    // gsDebugVar(1.0 / m_stretches(a));
    // gsDebugVar(dPsi_da(a));
    // gsDebugVar(p());
    // gsDebugVar(dJ_da(a));
    // T result = 1.0/m_stretches(a) * (dPsi_da(a) - p() * dJ_da(a) );
    // T result = 1.0/m_stretches(a) * dPsi_da(a);
    if (!m_compressible)
        return 1.0/m_stretches(a) * (dPsi_da(a) - p() * dJ_da(a) );
    else
        return 1.0/m_stretches(a) * dPsi_da(a);
}

template<class T>
T gsMaterialMatrix<T>::dSa_db(const index_t a, const index_t b) const
{
    T tmp = 0.0;
    if (!m_compressible)
    {
        tmp = 1.0/m_stretches(a) * ( d2Psi_dab(a,b) - dp_da(a)*dJ_da(b) - dp_da(b)*dJ_da(a) - p() * d2J_dab(a,b) );
        if (a==b)
            tmp += - 1.0 / math::pow(m_stretches(a),2) * (dPsi_da(a) - p() * dJ_da(a));
        return tmp;

        // return 1.0/m_stretches(a) * ( d2Psi_dab(a,b) - dp_da(a)*dJ_da(b) - dp_da(b)*dJ_da(a) - p() * d2J_dab(a,b) ) - delta(a,b) / math::pow(m_stretches(a),2) * (dPsi_da(a) - p() * dJ_da(a));
    }
    else
    {
        tmp = 1.0/m_stretches(a) * d2Psi_dab(a,b);
        if (a==b)
            tmp += - 1.0 / math::pow(m_stretches(a),2) * dPsi_da(a);
        return tmp;

        // return 1.0/m_stretches(a) * d2Psi_dab(a,b) - delta(a,b) / math::pow(m_stretches(a),2) * dPsi_da(a);
    }
}

template<class T>
T gsMaterialMatrix<T>::Cabcd(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    // Compute part with stress tensor involved.
    T frac = 0.0;
    T tmp = 0.0;
    // if ( (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-12) && (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) > 1e-14) )
    //     gsInfo<<"Warning: difference in stretches is close to machine precision?\n";

    if (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-14)
    {
        // gsDebug<<"Stretches are equal; (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) = "<<abs((m_stretches(a) - m_stretches(b)) / m_stretches(a))<<"\n";
        frac = 1.0 / (2.0 * m_stretches(a) ) * ( dSa_db(b,b) - dSa_db(a,b));
    }
    else
        frac = ( Sa(b)-Sa(a) ) / (math::pow(m_stretches(b),2) - math::pow(m_stretches(a),2));


    // Compute result
    if (!m_compressible)
    {
        GISMO_ASSERT( ( (a < 2) && (b < 2) && (c < 2) && (d < 2) ) , "Index out of range. a="<<a<<", b="<<b<<", c="<<c<<", d="<<d);
        if ( ( (a==b) && (c==d)) )
            tmp = 1/m_stretches(c) * dSa_db(a,c) + 1/(math::pow(m_stretches(a),2) * math::pow(m_stretches(c),2)) * ( math::pow(m_stretches(2),2) * d2Psi_dab(2,2) + 2*dPsi_da(2)*m_stretches(2) );
        else if (( (a==d) && (b==c) && (a!=b) ) || ( ( (a==c) && (b==d) && (a!=b)) ))
            tmp = frac;

        // return 1/m_stretches(c) * dSa_db(a,c) * delta(a,b) * delta(c,d) + frac * (delta(a,c)*delta(b,d) + delta(a,d)*delta(b,c)) * (1-delta(a,b))
                // + delta(a,b)*delta(c,d)*1/(math::pow(m_stretches(a),2) * math::pow(m_stretches(c),2)) * ( math::pow(m_stretches(2),2) * d2Psi_dab(2,2) + 2*dPsi_da(2)*m_stretches(2) );
    }
    else
    {
        GISMO_ASSERT( ( (a < 3) && (b < 3) && (c < 3) && (d < 3) ) , "Index out of range. a="<<a<<", b="<<b<<", c="<<c<<", d="<<d);
        if ( ( (a==b) && (c==d)) )
            tmp = 1/m_stretches(c) * dSa_db(a,c);
        else if (( (a==d) && (b==c) && (a!=b) ) || ( ( (a==c) && (b==d) && (a!=b)) ))
            tmp = frac;
        // return 1/m_stretches(c) * dSa_db(a,c) * delta(a,b) * delta(c,d) + frac * (delta(a,c)*delta(b,d) + delta(a,d)*delta(b,c)) * (1-delta(a,b));
    }
    return tmp;
}


// template<class T>
// T gsMaterialMatrix<T>::dPsi(const index_t a) const
// {
//     T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));

//     GISMO_ASSERT( ( (a < 3) ) , "Index out of range. a="<<a);
//     GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

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

// template<class T>
// T gsMaterialMatrix<T>::d2Psi(const index_t a, const index_t b) const
// {
//     T mu = m_parvals.at(0) / (2. * (1. + m_parvals.at(1)));

//     GISMO_ASSERT( ( (a < 3) && (b < 3) ) , "Index out of range. a="<<a<<", b="<<b);
//     GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

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

template<class T>
void gsMaterialMatrix<T>::computeMetricDeformed() const
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

template<class T>
void gsMaterialMatrix<T>::computeMetricUndeformed() const
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

template<class T>
void gsMaterialMatrix<T>::getMetric(index_t k, T z) const
{
    this->getMetricDeformed(k,z);
    this->getMetricUndeformed(k,z);

    T ratio = m_Gcov_def.determinant() / m_Gcov_ori.determinant();
    GISMO_ASSERT(ratio > 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_Gcov_ori.determinant());
    m_J0_sq = ratio;
}

template<class T>
void gsMaterialMatrix<T>::getMetricDeformed(index_t k, T z) const
{
    GISMO_ASSERT(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_ncov_def_mat.cols()!=0,"Is the basis initialized?");

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

template<class T>
void gsMaterialMatrix<T>::getMetricUndeformed(index_t k, T z) const
{
    GISMO_ASSERT(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");

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


template<class T>
std::pair<gsVector<T>,gsMatrix<T>> gsMaterialMatrix<T>::evalStretch(const gsMatrix<T> & C) const
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

template<class T>
void gsMaterialMatrix<T>::computeStretch(const gsMatrix<T> & C) const
{
    std::pair<gsVector<T>,gsMatrix<T>> result = evalStretch(C);
    m_stretches = result.first;
    m_stretchvec = result.second;
}


} // end namespace
