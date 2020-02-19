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
    To Do:
    - Make computeMetric() depending on the material model.
*/



#pragma once

#include <gsThinShell2/gsMaterialMatrix.h>

namespace gismo
{

// Linear material models; no strain computation
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & par1,
                                        const gsFunction<T> & par2
                                        )
                                        :
                                        m_patches(&mp),
                                        m_thickness(&thickness),
                                        m_par1(&par1),
                                        m_par2(&par2),
                                        m_piece(nullptr)
{
    initialize();
}
// Linear material models; no strain computation
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & par1,
                                        const gsFunction<T> & par2,
                                        const gsFunction<T> & density
                                        )
                                        :
                                        m_patches(&mp),
                                        m_thickness(&thickness),
                                        m_par1(&par1),
                                        m_par2(&par2),
                                        m_density(&density),
                                        m_piece(nullptr)
{
    initialize();
}

// Linear material models
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & mp_def,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & par1,
                                        const gsFunction<T> & par2
                                        )
                                        :
                                        m_patches(&mp),
                                        m_defpatches(&mp_def),
                                        m_thickness(&thickness),
                                        m_par1(&par1),
                                        m_par2(&par2),
                                        m_piece(nullptr)
{
    initialize();
}
// Linear material models
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & mp_def,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & par1,
                                        const gsFunction<T> & par2,
                                        const gsFunction<T> & density
                                        )
                                        :
                                        m_patches(&mp),
                                        m_defpatches(&mp_def),
                                        m_thickness(&thickness),
                                        m_par1(&par1),
                                        m_par2(&par2),
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
}


template <class T>
void gsMaterialMatrix<T>::computePoints(const gsMatrix<T> & u, bool deformed) const
{
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
        m_par1->eval_into(m_map.values[0], m_par1mat);
        m_par2->eval_into(m_map.values[0], m_par2mat);
    }
}

template <class T>
void gsMaterialMatrix<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->getOptions();
    if ((m_outputType==1) || (m_outputType==2)) // for matrix and vector
    {
        if (m_integration==0) // DD method
            this->eval_into_AP(u,result);
            // this->eval_into_DD(u,result);
        else if (m_integration==1) // AP method
            this->eval_into_AP(u,result);
            // this->eval_into_AP(u,result);
        else if (m_integration==2) // NP method
            this->eval_into_AP(u,result);
    }
    // This is for density and thickness output
    else if (m_outputType==0)
        this->eval_into_dens(u,result);
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
void gsMaterialMatrix<T>::eval_into_NP(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    // Define the moment to take
    if      (m_outputType==1) // output is a vector
        if      (m_output==0) // A matrix
            m_moment=0;
        else if (m_output==1) // B matrix
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
    else
        GISMO_ERROR("Something went wrong");

    this->computePoints(u);
    if (m_material==0)
    {
        result = integrateZ(u,m_moment);
    }
    else if (m_material==1)
        result = eval_Composite(u,m_moment); // thickness integration is done inside
    else
        result = integrateZ(u,m_moment);
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
        gsDebug<<"No eval3D function for composites available. Something went wrong...";

    else if ((m_material == 2) && (m_compressible)) // NHK
        result = eval_Compressible(i, z);
    else if ((m_material == 2) && (!m_compressible)) // NHK
        result = eval_Incompressible(i, z);
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
        gsWarn<<"no function available.";
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
template <class T>
gsMatrix<T> gsMaterialMatrix<T>::multiplyZ(const gsMatrix<T>& u, const index_t moment) const
{
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

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval_Incompressible(const gsMatrix<T> &u) const
{
    // Input: 
    //      u: evaluation points
    gsMatrix<T> result(this->targetDim(), u.cols());
    result.setZero();

    for( index_t k=0; k < u.cols(); ++k ) // evaluation points
    {
        // Evaluate material properties on the quadrature point
        m_par1val = m_par1mat(0,k);
        m_par2val = m_par2mat(0,k);

        if (m_outputType==2)
        {
            // NOTE: compute metric A and metric B only!! G does not make sense to compute.

            // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
            this->getMetric(i,z.at(j)); // on point i, on height z(0,j)

            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j,3,3);
            /*
                C = C1111,  C1122,  C1112
                    symm,   C2222,  C2212
                    symm,   symm,   C1212
            */
            if (m_output==0)
            {
                C(0,0)          = Aijkl(0,0,0,0); // C1111
                C(1,1)          = Aijkl(1,1,1,1); // C2222
                C(2,2)          = Aijkl(0,1,0,1); // C1212
                C(1,0) = C(0,1) = Aijkl(0,0,1,1); // C1122
                C(2,0) = C(0,2) = Aijkl(0,0,0,1); // C1112
                C(2,1) = C(1,2) = Aijkl(1,1,0,1); // C2212
            }
            else if (m_output==1)
            {
                C(0,0)          = Bijkl(0,0,0,0); // C1111
                C(1,1)          = Bijkl(1,1,1,1); // C2222
                C(2,2)          = Bijkl(0,1,0,1); // C1212
                C(1,0) = C(0,1) = Bijkl(0,0,1,1); // C1122
                C(2,0) = C(0,2) = Bijkl(0,0,0,1); // C1112
                C(2,1) = C(1,2) = Bijkl(1,1,0,1); // C2212
            }
            else if (m_output==2)
            {
                C(0,0)          = Cijkl(0,0,0,0); // C1111
                C(1,1)          = Cijkl(1,1,1,1); // C2222
                C(2,2)          = Cijkl(0,1,0,1); // C1212
                C(1,0) = C(0,1) = Cijkl(0,0,1,1); // C1122
                C(2,0) = C(0,2) = Cijkl(0,0,0,1); // C1112
                C(2,1) = C(1,2) = Cijkl(1,1,0,1); // C2212
            }
            else if (m_output==3)
            {
                C(0,0)          = Dijkl(0,0,0,0); // C1111
                C(1,1)          = Dijkl(1,1,1,1); // C2222
                C(2,2)          = Dijkl(0,1,0,1); // C1212
                C(1,0) = C(0,1) = Dijkl(0,0,1,1); // C1122
                C(2,0) = C(0,2) = Dijkl(0,0,0,1); // C1112
                C(2,1) = C(1,2) = Dijkl(1,1,0,1); // C2212
            }
        }
        else if (m_outputType==1)
        {
            // this->computeMetric(i,z.at(j),true,true);
            this->getMetric(i,z.at(j)); // on point i, on height z(0,j)
            if (m_output==0)
            {
                result(0,j) = Sij(0,0);
                result(1,j) = Sij(1,1);
                result(2,j) = Sij(0,1);   
            }
            else if (m_output==1)
            {
                result(0,j) = Mij(0,0);
                result(1,j) = Mij(1,1);
                result(2,j) = Mij(0,1);   
            }
        }
        else
            GISMO_ERROR("no vector or matrix produced");
    }

    return result;
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
        m_par1val = m_par1mat(0,i);
        m_par2val = m_par2mat(0,i);

        // this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)
        this->getMetric(i,z.at(j)); // on point i, on height z(0,j)

        // Define objects
        gsMatrix<T,3,3> c, cinv;
        T S33, C3333, dc33;
        S33 = 0.0;
        dc33 = 0.0;
        C3333 = 1.0;

        index_t itmax = 100;
        T tol = 1e-6;

        // Initialize c
        c.setZero();
        c.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        c(2,2) = 1.0; // c33
        cinv.block(0,0,2,2) = m_Gcon_def.block(0,0,2,2);
        cinv(2,2) = 1.0;

        for (index_t it = 0; it < itmax; it++)
        {
            dc33 = -2. * S33 / C3333;
            c(2,2) += dc33;

            GISMO_ASSERT(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33);
            cinv(2,2) = 1.0/c(2,2);

            m_J = m_J0 * math::sqrt( c(2,2) );

            S33     = Sij(2,2,c,cinv);
            C3333   = Cijkl3D(2,2,2,2,c,cinv);

            if (abs(S33) < tol)
            {
                // gsInfo<<"Converged in "<<it<<" iterations, abs(S33) = "<<abs(S33)<<" and tolerance = "<<tol<<"\n";
                if (m_outputType==2)
                    {
                        gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j,3,3);
                        /*
                            C = C1111,  C1122,  C1112
                                symm,   C2222,  C2212
                                symm,   symm,   C1212
                        */
                        if (m_output==0)
                        {
                            C(0,0)          = Aijkl(0,0,0,0,c,cinv); // C1111
                            C(1,1)          = Aijkl(1,1,1,1,c,cinv); // C2222
                            C(2,2)          = Aijkl(0,1,0,1,c,cinv); // C1212
                            C(1,0) = C(0,1) = Aijkl(0,0,1,1,c,cinv); // C1122
                            C(2,0) = C(0,2) = Aijkl(0,0,0,1,c,cinv); // C1112
                            C(2,1) = C(1,2) = Aijkl(1,1,0,1,c,cinv); // C2212
                        }
                        else if (m_output==1)
                        {
                            C(0,0)          = Bijkl(0,0,0,0,c,cinv); // C1111
                            C(1,1)          = Bijkl(1,1,1,1,c,cinv); // C2222
                            C(2,2)          = Bijkl(0,1,0,1,c,cinv); // C1212
                            C(1,0) = C(0,1) = Bijkl(0,0,1,1,c,cinv); // C1122
                            C(2,0) = C(0,2) = Bijkl(0,0,0,1,c,cinv); // C1112
                            C(2,1) = C(1,2) = Bijkl(1,1,0,1,c,cinv); // C2212
                        }
                        else if (m_output==2)
                        {
                            C(0,0)          = Cijkl(0,0,0,0,c,cinv); // C1111
                            C(1,1)          = Cijkl(1,1,1,1,c,cinv); // C2222
                            C(2,2)          = Cijkl(0,1,0,1,c,cinv); // C1212
                            C(1,0) = C(0,1) = Cijkl(0,0,1,1,c,cinv); // C1122
                            C(2,0) = C(0,2) = Cijkl(0,0,0,1,c,cinv); // C1112
                            C(2,1) = C(1,2) = Cijkl(1,1,0,1,c,cinv); // C2212
                        }
                        else if (m_output==3)
                        {
                            C(0,0)          = Dijkl(0,0,0,0,c,cinv); // C1111
                            C(1,1)          = Dijkl(1,1,1,1,c,cinv); // C2222
                            C(2,2)          = Dijkl(0,1,0,1,c,cinv); // C1212
                            C(1,0) = C(0,1) = Dijkl(0,0,1,1,c,cinv); // C1122
                            C(2,0) = C(0,2) = Dijkl(0,0,0,1,c,cinv); // C1112
                            C(2,1) = C(1,2) = Dijkl(1,1,0,1,c,cinv); // C2212
                        }
                    }
                    else if (m_outputType==1)
                    {
                        if (m_output==0)
                        {
                            result(0,j) = Sij(0,0,c,cinv);
                            result(1,j) = Sij(1,1,c,cinv);
                            result(2,j) = Sij(0,1,c,cinv);   
                        }
                        else if (m_output==1)
                        {
                            result(0,j) = Mij(0,0,c,cinv);
                            result(1,j) = Mij(1,1,c,cinv);
                            result(2,j) = Mij(0,1,c,cinv);   
                        }
                    }
                else
                    GISMO_ERROR("no vector or matrix produced");

                    break;
            }


            else if (it == itmax - 1)
            {
                gsInfo<<"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";
                // std::terminate();
            }
        }
    }
    return result;
}

/*
    Available class members:
        - m_par1val
        - m_par2val
        - m_metric
        - m_metric_def
        - m_J0
        - m_J
*/
// template<class T>
// T gsMaterialMatrix<T>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const { Cijkl(i,j,k,l,NULL,NULL); }

template<class T>
T gsMaterialMatrix<T>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ASSERT( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    T tmp = 0.0;
    if (m_material==0) // svk
    {
        T lambda, mu, Cconstant;

        mu = m_par1val / (2.*(1. + m_par2val));
        lambda = m_par1val * m_par2val / ( (1. + m_par2val)*(1.-2.*m_par2val)) ;
        Cconstant = 2*lambda*mu/(lambda+2*mu);

        tmp = Cconstant*m_Acon_ori(i,j)*m_Acon_ori(k,l) + mu*(m_Acon_ori(i,k)*m_Acon_ori(j,l) + m_Acon_ori(i,l)*m_Acon_ori(j,k));
    }
    else if (m_material==1)
        gsWarn<<"Compressible material matrix  requested, but not needed. How?";
    else if (m_material==9)
    {
        T mu = m_par1val / (2.*(1. + m_par2val));
        tmp = mu*1./math::pow(m_J0,2.)*(2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k));
    }

    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else // general version
    {
        tmp = 4.0 * d2Psi(i,j,k,l) + 4.0 * d2Psi(2,2,2,2)*math::pow(m_J0,-4.0)*m_Gcon_def(i,j)*m_Gcon_def(k,l)
                - 4.0 * d2Psi(2,2,i,j)*math::pow(m_J0,-2.0)*m_Gcon_def(k,l) - 4.0 * d2Psi(2,2,k,l)*math::pow(m_J0,-2.0)*m_Gcon_def(i,j)
                + 2.0 * dPsi(2,2) * math::pow(m_J0,-2.)*(2.*m_Gcon_def(i,j)*m_Gcon_def(k,l) + m_Gcon_def(i,k)*m_Gcon_def(j,l) + m_Gcon_def(i,l)*m_Gcon_def(j,k));

    }

    return tmp;
}

// Consensation of the 3D tensor for compressible materials
template<class T>
T gsMaterialMatrix<T>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ASSERT(c.cols()==c.rows(),"Matrix c must be square");
    GISMO_ASSERT(c.cols()==3,"Matrix c must be 3x3");
    GISMO_ASSERT(cinv.cols()==cinv.rows(),"Matrix cinv must be square");
    GISMO_ASSERT(cinv.cols()==3,"Matrix cinv must be 3x3");
    GISMO_ASSERT( ( (i <2) && (j <2) && (k <2) && (l <2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(m_compressible,"Material model is not compressible?");

    return Cijkl3D(i,j,k,l,c,cinv) - ( Cijkl3D(i,j,2,2,c,cinv) * Cijkl3D(2,2,k,l,c,cinv) ) / Cijkl3D(2,2,2,2,c,cinv);
}

// 3D tensor for compressible materials
template<class T>
T gsMaterialMatrix<T>::Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ASSERT( ( (i <3) && (j <3) && (k <3) && (l <3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);

    T tmp = 0.0;
    if (m_material==0 || m_material==1) // svk
        gsWarn<<"Compressible material matrix requested, but not needed. How?";
    else if (m_material==9)
    {
        T mu = m_par1val / (2 * (1 + m_par2val));
        T K  = m_par1val / ( 3 - 6 * m_par2val);
        T traceC =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                    m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                    m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                    m_Gcov_def(1,1)*m_Gcon_ori(1,1) +
                    c(2,2);

        tmp = 1.0 / 9.0 * mu * math::pow( m_J , -2.0/3.0 ) * ( traceC * ( 2.0*cinv(i,j)*cinv(k,l) + 3.0*cinv(i,k)*cinv(j,l) + 3.0*cinv(i,l)*cinv(j,k) )
                            - 6.0 *( m_Gcon_ori(i,j)*cinv(k,l) + cinv(i,j)*m_Gcon_ori(k,l) ) ) + K * ( m_J*m_J*cinv(i,j)*cinv(k,l) - 0.5*(m_J*m_J-1)*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) ) );
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else
        tmp = 4.0 * d2Psi(i,j,k,l,c,cinv);

    return tmp;
}



/*
    Available class members:
        - m_par1val
        - m_par2val
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
    if (m_material==0)
    {
        gsMatrix<T> strain;
        // --------------------------
        // Saint Venant Kirchhoff
        // --------------------------
        stress = 0.5 * (m_Acov_def - m_Acov_ori);
        tmp =     Cijkl(i,j,0,0) * strain(0,0) + Cijkl(i,j,0,1) * strain(0,1)
                + Cijkl(i,j,1,0) * strain(1,0) + Cijkl(i,j,1,1) * strain(1,1);
    }
    else if (m_material==1)
    {
        // --------------------------
        // Composite
        // --------------------------
        gsWarn<<"Incompressible material stress tensor requested, but not needed. How?";
    }
    else if (m_material==9)
    {
        // --------------------------
        // Neo-Hoookean
        // --------------------------
        // T mu = m_par1val / (2. * (1. + m_par2val));

        // gsDebugVar(mu);
        // gsDebugVar(mu * (m_Gcon_ori(i,j) - 1./math::pow(m_J0,2.) * m_Gcon_def(i,j) ));
        // tmp = mu * (m_Gcon_ori(i,j) - 1./math::pow(m_J0,2.) * m_Gcon_def(i,j) );
    }
    else
    {
        // --------------------------
        // Generalized
        // --------------------------
        // tmp = 2.0 * dPsi(i,j) - 2.0 * dPsi(2,2) * math::pow(m_J0,-2.0)*m_Gcon_def(i,j);
    }
    return tmp;
}

template<class T>
T gsMaterialMatrix<T>::dSij(const index_t i, const index_t j) const
{
    T tmp = 0.0;
    if (m_material==0)
    {
        gsMatrix<T> strain;
        // --------------------------
        // Saint Venant Kirchhoff
        // --------------------------
        stress = 0.5 * (m_Gcov_def - m_Gcov_ori);
        tmp =     Cijkl(i,j,0,0) * strain(0,0) + Cijkl(i,j,0,1) * strain(0,1)
                + Cijkl(i,j,1,0) * strain(1,0) + Cijkl(i,j,1,1) * strain(1,1);
    }
    else if (m_material==1)
    {
        // --------------------------
        // Composite
        // --------------------------
        gsWarn<<"Incompressible material stress tensor requested, but not needed. How?";
    }
    else if (m_material==9)
    {
        // --------------------------
        // Neo-Hoookean
        // --------------------------
        // T mu = m_par1val / (2. * (1. + m_par2val));

        // gsDebugVar(mu);
        // gsDebugVar(mu * (m_Gcon_ori(i,j) - 1./math::pow(m_J0,2.) * m_Gcon_def(i,j) ));
        // tmp = mu * (m_Gcon_ori(i,j) - 1./math::pow(m_J0,2.) * m_Gcon_def(i,j) );
    }
    else
    {
        // --------------------------
        // Generalized
        // --------------------------
        // tmp = 2.0 * dPsi(i,j) - 2.0 * dPsi(2,2) * math::pow(m_J0,-2.0)*m_Gcon_def(i,j);
    }
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
    if (m_material==0 || m_material==1)
        gsWarn<<"Incompressible material stress tensor requested, but not needed. How?";
    else if (m_material==9)
    {
        // T mu = m_par1val / (2 * (1 + m_par2val));
        // T K  = m_par1val / ( 3 - 6 * m_par2val);
        // T traceC =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
        //             m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
        //             m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
        //             m_Gcov_def(1,1)*m_Gcon_ori(1,1) +
        //             c(2,2);

        // tmp =  mu * math::pow( m_J , -2.0/3.0 ) * ( m_Gcon_ori(i,j) - 1.0/3.0 * traceC * cinv(i,j) ) + 0.5 * K * ( m_J*m_J - 1 ) * cinv(i,j);
    }
    else
        // tmp = 2.0 * dPsi(i,j,c,cinv);

    return tmp;
}

// ---------------------------------------------------------------------------------------------------------------------------------
//                                          INCOMPRESSIBLE
// ---------------------------------------------------------------------------------------------------------------------------------
template<class T>
T gsMaterialMatrix<T>::dPsi(const index_t i, const index_t j) const
{
    T mu = m_par1val / (2. * (1. + m_par2val));

    GISMO_ASSERT( ( (i < 3) && (j < 3) ) , "Index out of range. i="<<i<<", j="<<j);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material==2)
    {
        return 0.5 * mu * m_Gcon_ori(i,j);
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else if (m_material==4)
    {
        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = math::pow(m_J0,-2.0);
        computeStretch(C);

        T result = 0.0;
        for (index_t a = 0; a != m_stretches.rows(); a++)
            result += dPsi(a)*(
                        ( m_gcov_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcov_ori.col(j).dot(m_stretchvec.col(a)) )
                        );
        return  result;
    }
    else
        GISMO_ERROR("Material model not implemented.");
}

template<class T>
T gsMaterialMatrix<T>::d2Psi(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    T mu = m_par1val / (2. * (1. + m_par2val));

    GISMO_ASSERT( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material==2)
    {
        return 0.0;
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else if (m_material==4)
    {
        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = math::pow(m_J0,-2.0);
        computeStretch(C);

        T result = 0.0;
        for (index_t a = 0; a != m_stretches.rows(); a++)
            for (index_t b = 0; b != m_stretches.rows(); b++)
                result += d2Psi(a,b)*(
                            ( m_gcov_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcov_ori.col(j).dot(m_stretchvec.col(a)) )*
                            ( m_gcov_ori.col(k).dot(m_stretchvec.col(b)) )*( m_gcov_ori.col(l).dot(m_stretchvec.col(b)) )
                            );
        return  result;
    }
    else
        GISMO_ERROR("Material model not implemented.");
}
// ---------------------------------------------------------------------------------------------------------------------------------
//                                          COMPRESSIBLE
// ---------------------------------------------------------------------------------------------------------------------------------
template<class T>
T gsMaterialMatrix<T>::dPsi(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T mu = m_par1val / (2. * (1. + m_par2val));
    T K  = m_par1val / ( 3 - 6 * m_par2val);
    T traceC =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                m_Gcov_def(1,1)*m_Gcon_ori(1,1) +
                c(2,2);

    GISMO_ASSERT( ( (i < 3) && (j < 3) ) , "Index out of range. i="<<i<<", j="<<j);
    GISMO_ASSERT(m_compressible,"Material model is not compressible?");

    if (m_material==2)
    {
        T tmp = 0.5 * mu * math::pow(m_J,-2./3.) * ( m_Gcon_ori(i,j) - 1.0/3.0 * traceC * cinv(i,j) ) + 0.25 * K * (m_J*m_J - 1.0) * cinv(i,j);
        return tmp;
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else
        GISMO_ERROR("Material model not implemented (Cijkl.");
}

template<class T>
T gsMaterialMatrix<T>::d2Psi(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
   T mu = m_par1val / (2 * (1 + m_par2val));
   T K  = m_par1val / ( 3 - 6 * m_par2val);
   T traceC =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
               m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
               m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
               m_Gcov_def(1,1)*m_Gcon_ori(1,1) +
               c(2,2);

    GISMO_ASSERT( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(m_compressible,"Material model is not compressible?");

    if (m_material==2)
    {
        T tmp = 1.0 / 36.0 * mu * math::pow( m_J , -2.0/3.0 ) * ( traceC * ( 2.0*cinv(i,j)*cinv(k,l) + 3.0*cinv(i,k)*cinv(j,l) + 3.0*cinv(i,l)*cinv(j,k) )
                    - 6.0 *( m_Gcon_ori(i,j)*cinv(k,l) + cinv(i,j)*m_Gcon_ori(k,l) ) )
                    + 0.25 * K * ( m_J*m_J*cinv(i,j)*cinv(k,l) - 0.5*(m_J*m_J-1.0)*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) ) );
        return tmp;
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else
        GISMO_ERROR("Material model not implemented.");
}
// ---------------------------------------------------------------------------------------------------------------------------------
//                                          STRETCHES
// ---------------------------------------------------------------------------------------------------------------------------------
template<class T>
T gsMaterialMatrix<T>::dPsi(const index_t a) const
{
    T mu = m_par1val / (2. * (1. + m_par2val));

    GISMO_ASSERT( ( (a < 3) ) , "Index out of range. a="<<a);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material==2)
    {
        return mu * m_stretches(a,0);
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else if (m_material==4)
    {

    }
    else
        GISMO_ERROR("Material model not implemented.");
}

template<class T>
T gsMaterialMatrix<T>::d2Psi(const index_t a, const index_t b) const
{
    T mu = m_par1val / (2. * (1. + m_par2val));

    GISMO_ASSERT( ( (a < 3) && (b < 3) ) , "Index out of range. a="<<a<<", b="<<b);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material==2)
    {
        return ( a==b ? mu : 0.0 );
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
    }
    else
        GISMO_ERROR("Material model not implemented.");
}
// ---------------------------------------------------------------------------------------------------------------------------------
//                                          Metric Computations
// ---------------------------------------------------------------------------------------------------------------------------------

template<class T>
void gsMaterialMatrix<T>::computeMetricDeformed() const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;

    m_Acov_def_mat.resize(4,m_map_def.points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.points.cols());    m_Acon_def_mat.setZero();
    m_Bcov_def_mat.resize(4,m_map_def.points.cols());    m_Bcov_def_mat.setZero();

    m_acov_def_mat.resize(2*3,m_map_def.points.cols());    m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*3,m_map_def.points.cols());    m_acon_def_mat.setZero();
    m_ncov_def_mat.resize(2*3,m_map_def.points.cols());    m_ncov_def_mat.setZero();

    gsMatrix<T> tmp;

    for (index_t k=0; k!= m_map_def.points.cols(); k++)
    {
        // covariant basis vectors
        m_acov_def_mat.reshapeCol(k,3,2)   = m_map_def.jacobian(k);

        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,3,2);
        // Construct metric tensor a = [dcd1*dcd1, dcd1*dcd2; dcd2*dcd1, dcd2*dcd2]
        tmp                         = acov.transpose() * acov;

        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricAcov = m_Acov_def_mat.reshapeCol(k,2,2);
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();
        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        // contravariant basis vectors
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

        // Derivative of the normal
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

    m_Acov_ori_mat.resize(4,m_map.points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map.points.cols());    m_Acon_ori_mat.setZero();
    m_Bcov_ori_mat.resize(4,m_map.points.cols());    m_Bcov_ori_mat.setZero();

    m_acov_ori_mat.resize(2*3,m_map.points.cols());    m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*3,m_map.points.cols());    m_acon_ori_mat.setZero();
    m_ncov_ori_mat.resize(2*3,m_map.points.cols());    m_ncov_ori_mat.setZero();

    gsMatrix<T> tmp;

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        // covariant basis vectors
        m_acov_ori_mat.reshapeCol(k,3,2)   = m_map.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,3,2);
        // Construct metric tensor a = [dcd1*dcd1, dcd1*dcd2; dcd2*dcd1, dcd2*dcd2]
        tmp                         = acov.transpose() * acov;

        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricAcov = m_Acov_ori_mat.reshapeCol(k,2,2);
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();
        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        // contravariant basis vectors
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

        // Derivative of the normal
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

    // this->getBasisDeformed(k,z);
    // this->getBasisUndeformed(k,z);
    this->getBasis(k,z);

    m_J0 = math::sqrt( m_Gcov_def.determinant() / m_Gcov_ori.determinant() );
}

template<class T>
void gsMaterialMatrix<T>::getMetricDeformed(index_t k, T z) const
{
    GISMO_ASSERT(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");

    m_Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    m_Acon_def = m_Acon_def_mat.reshapeCol(k,2,2);
    m_Bcov_def = m_Bcov_def_mat.reshapeCol(k,2,2);

    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def - 2.0 * z * m_Bcov_def;
    m_Gcov_def(2,2) = 1.0;

    m_Gcon_def = m_Gcov_def.inverse();
}

template<class T>
void gsMaterialMatrix<T>::getMetricUndeformed(index_t k, T z) const
{

    GISMO_ASSERT(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ASSERT(m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");

    m_Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    m_Acon_ori = m_Acon_ori_mat.reshapeCol(k,2,2);
    m_Bcov_ori = m_Bcov_ori_mat.reshapeCol(k,2,2);

    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori - 2.0 * z * m_Bcov_ori;
    m_Gcov_ori(2,2) = 1.0;

    m_Gcon_ori = m_Gcov_ori.inverse();
}

template<class T>
void gsMaterialMatrix<T>::getBasis(index_t k, T z) const
{
    this->getBasisDeformed(k,z);
    this->getBasisUndeformed(k,z);

    // Approx previously computed
    // gsDebugVar(m_Gcov_def);
    // Approx based on basis
    // gsDebugVar(m_gcov_def.transpose()*m_gcov_def);
    // Approx minus the truncated (quadratic) part
    // gsDebugVar(m_gcov_def.transpose()*m_gcov_def - (z * m_ncov_def).transpose()*(z * m_ncov_def));
}

template<class T>
void gsMaterialMatrix<T>::getBasisDeformed(index_t k, T z) const
{
    GISMO_ASSERT(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_ncov_def_mat.cols()!=0,"Is the basis initialized?");

    m_acov_def = m_acov_def_mat.reshapeCol(k,3,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,3,2);
    m_ncov_def = m_ncov_def_mat.reshapeCol(k,3,2);

    m_gcov_def = m_acov_def + z * m_ncov_def;
}

template<class T>
void gsMaterialMatrix<T>::getBasisUndeformed(index_t k, T z) const
{
    GISMO_ASSERT(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ASSERT(m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");

    m_acov_ori = m_acov_ori_mat.reshapeCol(k,3,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,3,2);
    m_ncov_ori = m_ncov_ori_mat.reshapeCol(k,3,2);

    m_gcov_ori = m_acov_ori + z * m_ncov_ori;
}

template<class T>
void gsMaterialMatrix<T>::computeBasisUndeformed() const
{
    m_acov_ori_mat.resize(2*3,m_map.points.cols());    m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*3,m_map.points.cols());    m_acon_ori_mat.setZero();
    m_ncov_ori_mat.resize(2*3,m_map.points.cols());    m_ncov_ori_mat.setZero();

    gsMatrix<T,2,2> mixedB;

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,3,2);
        acov = m_map.jacobian(k);

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_ori_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_ori_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

        gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_ori_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);

        // for testing
        gsMatrix<T> metricB(2,2);
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                metricB(i,j) = -ncov.col(j).transpose()*acov.col(i);
        // for testing
        gsMatrix<T,2,2> metricAcon2;
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                metricAcon2(i,j) = acon.col(i).transpose()*acon.col(j);

    }
}

template<class T>
void gsMaterialMatrix<T>::computeBasisDeformed() const
{
    m_acov_def_mat.resize(2*3,m_map.points.cols());    m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*3,m_map.points.cols());    m_acon_def_mat.setZero();
    m_ncov_def_mat.resize(2*3,m_map.points.cols());    m_ncov_def_mat.setZero();

    gsMatrix<T,2,2> mixedB;

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,3,2);
        acov = m_map.jacobian(k);

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_def_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

        gsMatrix<T,2,2> metricAcon2;
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                metricAcon2(i,j) = acon.col(i).transpose()*acon.col(j);

        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

        gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_def_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);

        gsMatrix<T> metricB(2,2);
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                metricB(i,j) = -ncov.col(j).transpose()*acov.col(i);

    }
}


template<class T>
void gsMaterialMatrix<T>::computeStretch(const gsMatrix<T> & C) const
{
    m_stretches.resize(3,1);
    m_stretchvec.resize(3,3);
    gsMatrix<> evs;
    Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

    eigSolver.compute(C);
    evs = eigSolver.eigenvalues();
    for (index_t r=0; r!=3; r++)
        m_stretches(r,0) = math::sqrt(evs(r,0));

    gsDebugVar(eigSolver.eigenvalues());
    gsDebugVar(eigSolver.eigenvectors());
}

// template<class T>
// gsMatrix<T> gsMaterialMatrix<T>::basisTransform(const gsMatrix<T> & basis1, const gsMatrix<T> & basis2) const
// {
//     GISMO_ASSERT(basis1.cols()==basis2.cols(),"Number of basis vectors should be the same! basis1.cols() = "<<basis1.cols()<<", basis2.cols() = "<<basis2.cols());
//     GISMO_ASSERT(basis1.rows()==basis2.rows(),"Basis dimensions should be the same! basis1.rows() = "<<basis1.rows()<<", basis2.rows() = "<<basis2.rows());
//     index_t nc = basis1.cols();

//     gsMatrix<T> result(nc,nc);

//     for (index_t i = 1; i != nc; i++)
//         for (index_t j = 1; j != nc; j++)
//         {

//         }


// }

} // end namespace