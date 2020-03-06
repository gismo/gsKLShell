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

int delta(const int a, const int b)
{
    return (a==b) ? 1 : 0;
}

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
        gsInfo<<"No eval3D function for composites available. Something went wrong...";
    else if ((m_material == 2) && (m_compressible)) // NHK
        result = eval_Compressible(i, z);
    else if ((m_material == 2) && (!m_compressible)) // NHK
        result = eval_Incompressible(i, z);

    else if ((m_material == 4) && (m_compressible)) // NHK
        result = eval_Compressible(i, z);
    else if ((m_material == 4) && (!m_compressible)) // NHK
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

// integrates the material matrix over its last coordinate (thickness)
template <class T>
gsMatrix<T> gsMaterialMatrix<T>::integrateZ(const gsMatrix<T>& u, const index_t moment) const
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
                res += quWeights.at(k) * math::pow(quNodes(0,k),moment) * evalPoints(i,k);
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
        m_par1val = m_par1mat(0,i);
        m_par2val = m_par2mat(0,i);

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

            gsDebugVar(C);
        }
        else if (m_outputType==1)
        {
            // this->computeMetric(i,z.at(j),true,true);
            this->getMetric(i,z.at(j)); // on point i, on height z(0,j)

            result(0,j) = Sij(0,0);
            result(1,j) = Sij(1,1);
            result(2,j) = Sij(0,1);
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
        - m_par1val
        - m_par2val
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
    else if (m_material==4)
    {
        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = math::pow(m_J0,-2.0);
        computeStretch(C);

        tmp = 0.0;
        index_t maxIdx = 2;
        gsDebug<<"g ori = \n"<<m_gcon_ori<<"\n";
        gsDebug<<"eivec = \n"<<m_stretchvec<<"\n";
        for (index_t a = 0; a != maxIdx; a++)
            for (index_t b = 0; b != maxIdx; b++)
                for (index_t c = 0; c != maxIdx; c++)
                    for (index_t d = 0; d != maxIdx; d++)
                    {
                        if (((a==b)&&(c==d)) || (((a==c)&&(b==d)) && (a!=b)) || (((a==d)&&(b==c)) && (a!=b)))
                        {
                            // gsDebug<<"i = "<<i<<"\tj = "<<j<<"\tk = "<<k<<"\tl = "<<l<<"\ta = "<<a<<"\tb = "<<b<<"\tc = "<<c<<"\td = "<<d<<"\n";

                            // gsDebugVar(Cabcd(a,b,c,d));

                            // gsDebugVar(( m_gcov_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcov_ori.col(j).dot(m_stretchvec.col(b)) )*
                            //            ( m_gcov_ori.col(k).dot(m_stretchvec.col(c)) )*( m_gcov_ori.col(l).dot(m_stretchvec.col(d)) ));
                            // gsDebugVar(d2Psi(a,b));
                            // gsDebug<<"i = "<<i<<"\tj = "<<j<<"\tk = "<<k<<"\tl = "<<l<<"\ta = "<<a<<"\tb = "<<b<<"\n";
                            // gsDebugVar(( ( m_gcov_ori.col(i).dot(m_stretchvec.col(a)) ) ));
                            // gsDebugVar(( ( m_gcov_ori.col(j).dot(m_stretchvec.col(a)) ) ));
                            // gsDebugVar(( ( m_gcov_ori.col(k).dot(m_stretchvec.col(b)) ) ));
                            // gsDebugVar(( ( m_gcov_ori.col(l).dot(m_stretchvec.col(b)) ) ));

                            // gsDebugVar(( ( m_gcov_ori.col(i) ) ));
                            // gsDebugVar(( ( m_gcov_ori.col(j) ) ));
                            // gsDebugVar(( ( m_gcov_ori.col(k) ) ));
                            // gsDebugVar(( ( m_gcov_ori.col(l) ) ));

                            // gsDebugVar(( ( (m_stretchvec.col(b)) ) ));
                            // gsDebugVar(( ( (m_stretchvec.col(b)) ) ));

                            T fac = ( m_gcov_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcov_ori.col(j).dot(m_stretchvec.col(b)) )*
                                    ( m_gcov_ori.col(k).dot(m_stretchvec.col(c)) )*( m_gcov_ori.col(l).dot(m_stretchvec.col(d)) );

                            gsDebug<<"a = "<<a<<"; b = "<<b<<"; c = "<<c<<"; d = "<<d<<"; Cabcd = "<<Cabcd(a,b,c,d)<<"; fac = "<<fac<<"; prod = "<<fac*Cabcd(a,b,c,d)<<"\n";
                            gsDebug<<"i = "<<i<<"; j = "<<j<<"; k = "<<k<<"; l = "<<l<<"\n";
                            gsDebugVar(m_gcon_ori.col(i).dot(m_stretchvec.col(a)));
                            gsDebugVar(m_gcon_ori.col(j).dot(m_stretchvec.col(b)));
                            gsDebugVar(m_gcon_ori.col(k).dot(m_stretchvec.col(c)));
                            gsDebugVar(m_gcon_ori.col(l).dot(m_stretchvec.col(d)));

                            tmp +=  Cabcd(a,b,c,d)*(
                                    ( m_gcon_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcon_ori.col(j).dot(m_stretchvec.col(b)) )*
                                    ( m_gcon_ori.col(k).dot(m_stretchvec.col(c)) )*( m_gcon_ori.col(l).dot(m_stretchvec.col(d)) )
                                    );

                            // gsDebug<<"----------------------------------------------------------------------\n";
                        }
                        else
                            continue;
                    }
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
        gsMatrix<T> stress;
        // --------------------------
        // Saint Venant Kirchhoff
        // --------------------------
        // if (m_moment==0)
        // {
        //     GISMO_ASSERT( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
        //     stress = 0.5*(m_Acov_def - m_Acov_ori);
        // }
        // else if (m_moment==2)
        // {
        //     GISMO_ASSERT( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
        //     // tmp = (m_Bcov_ori - m_Bcov_def);
        //     stress = (m_Bcov_def - m_Bcov_ori);
        // }
        // else
        // {
        //     GISMO_ERROR("Warning: no material model known in simplification");
        // }
        stress = 0.5 * (m_Gcov_def - m_Gcov_ori);

        tmp =     Cijkl(i,j,0,0) * stress(0,0) + Cijkl(i,j,0,1) * stress(0,1)
                + Cijkl(i,j,1,0) * stress(1,0) + Cijkl(i,j,1,1) * stress(1,1);
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
        T mu = m_par1val / (2. * (1. + m_par2val));

        // gsDebugVar(mu);
        // gsDebugVar(mu * (m_Gcon_ori(i,j) - 1./math::pow(m_J0,2.) * m_Gcon_def(i,j) ));
        tmp = mu * (m_Gcon_ori(i,j) - 1./math::pow(m_J0,2.) * m_Gcon_def(i,j) );
    }
    else if (m_material==4)
    {
        gsMatrix<T> C(3,3);
        C.setZero();
        C.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        C(2,2) = math::pow(m_J0,-2.0);
        computeStretch(C);

        tmp = 0.0;
        for (index_t a = 0; a != 2; a++)
        {

            T fac = ( m_gcov_ori.col(i).dot(m_stretchvec.col(a)) )*
                    ( m_gcov_ori.col(j).dot(m_stretchvec.col(a)) );

            // gsDebug<<"a = "<<a<<"; Sa = "<<Sa(a)<<"; fac = "<<fac<<"; prod = "<<fac*Sa(a)<<"\n";
            // gsDebug<<"i = "<<i<<"; j = "<<j<<"\n";
            // gsDebug<<"g ori = \n"<<m_gcov_ori<<"\n";
            // gsDebug<<"eivec = \n"<<m_stretchvec<<"\n";
            tmp += Sa(a)*(
                        ( m_gcov_ori.col(i).dot(m_stretchvec.col(a)) )*( m_gcov_ori.col(j).dot(m_stretchvec.col(a)) )
                        );
        }
    }
    else
    {
        // --------------------------
        // Generalized
        // --------------------------
        tmp = 2.0 * dPsi(i,j) - 2.0 * dPsi(2,2) * math::pow(m_J0,-2.0)*m_Gcon_def(i,j);
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
        T mu = m_par1val / (2 * (1 + m_par2val));
        T K  = m_par1val / ( 3 - 6 * m_par2val);
        T traceC =  m_Gcov_def(0,0)*m_Gcon_ori(0,0) +
                    m_Gcov_def(0,1)*m_Gcon_ori(0,1) +
                    m_Gcov_def(1,0)*m_Gcon_ori(1,0) +
                    m_Gcov_def(1,1)*m_Gcon_ori(1,1) +
                    c(2,2);

        tmp =  mu * math::pow( m_J , -2.0/3.0 ) * ( m_Gcon_ori(i,j) - 1.0/3.0 * traceC * cinv(i,j) ) + 0.5 * K * ( m_J*m_J - 1 ) * cinv(i,j);
    }
    else
        tmp = 2.0 * dPsi(i,j,c,cinv);

    return tmp;
}

// Stretch-based formulation

// template<class T>
// T gsMaterialMatrix<T>::Sa_S(const index_t a) const
// {
//     T tmp = 0.0;
//     if (m_material==4)
//     {
//         // --------------------------
//         // Neo-Hoookean
//         // --------------------------
//         T mu = m_par1val / (2. * (1. + m_par2val));

//         // dPsidLambda_a
//         T dpsi = delta(i,j) * mu * m_stretches(i);
//         // S_a
//         tmp = (1. / m_stretches(i)) * dpsi;
//     }
//     else
//     {
//         GISMO_ERROR("Not available");
//     }
//     return tmp;
// }

// template<class T>
// T gsMaterialMatrix<T>::dSij_S(const index_t i, const index_t j) const
// {
//     T tmp = 0.0;
//     if (m_material==0)
//     {
//         tmp = 0.0;
//     }
//     else if (m_material==1)
//     {
//         // --------------------------
//         // Composite
//         // --------------------------
//         gsWarn<<"Incompressible material stress tensor requested, but not needed. How?";
//     }
//     else if (m_material==4)
//     {
//         // --------------------------
//         // Neo-Hoookean
//         // --------------------------
//         T mu = m_par1val / (2. * (1. + m_par2val));

//         // dPsidLambda_a
//         T DpsiDa = mu * (math::pow(m_stretches(0),2)*math::pow(m_stretches(1),2)*math::pow(m_stretches(2),2)) / m_stretches(i);
//         T D2psiDab = DpsiDa/m_stretches(j);
//         if (i!=j)
//             D2psiDab *= 2;
//         // S_a
//         tmp = (1. / m_stretches(i)) * D2psiDab + (1. / math::pow(m_stretches(i),2)) * DpsiDa;
//     }
//     else
//     {
//         tmp = 0.0;
//     }
//     return tmp;
// }

// template<class T>
// T gsMaterialMatrix<T>::Cijkl_S(const index_t i, const index_t j, const index_t k, const index_t l) const
// {
//     T tmp = 0.0;
//     if (m_material==0)
//     {
//         tmp = 0.0;
//     }
//     else if (m_material==1)
//     {
//         // --------------------------
//         // Composite
//         // --------------------------
//         gsWarn<<"Incompressible material stress tensor requested, but not needed. How?";
//     }
//     else if (m_material==4)
//     {
//         // --------------------------
//         // Neo-Hoookean
//         // --------------------------
//         if ((i==j) && (k==l))
//             tmp = 1. / m_stretches(j) * dSij_S(i,j);
//         else if ( ( ((i==k) && (j==l)) && (i!=j)) || ( ((i==l) && (j==k)) && (i!=j)) )
//         {
//             if (m_stretches(i)==m_stretches(j))
//                 tmp = 1./2. * ( 1. / m_stretches(j) * dSij_S(j,j) - 1. / m_stretches(i) * dSij_S(j,j) );
//             else
//                 tmp = (Sij_S(j,j)-Sij_S(i,i)) / (m_stretches(j) * m_stretches(j) - m_stretches(i) * m_stretches(i));
//         }
//         else
//             tmp = 0.0;
//     }
//     else
//     {
//         tmp = 0.0;
//     }
//     return tmp;
// }

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
        T S33_old;
        S33 = 0.0;
        dc33 = 0.0;
        C3333 = 1.0;

        index_t itmax = 10;
        T tol = 1e-10;

        // Initialize c
        c.setZero();
        c.block(0,0,2,2) = m_Gcov_def.block(0,0,2,2);
        c(2,2) = math::pow(m_J0,-2.0);; // c33
        cinv.block(0,0,2,2) = m_Gcon_def.block(0,0,2,2);
        cinv(2,2) = 1.0/c(2,2);

        m_J = m_J0 * math::sqrt( c(2,2) );
        S33 = Sij(2,2,c,cinv);
        S33_old = (S33 == 0.0) ? 1.0 : S33;
        C3333   = Cijkl3D(2,2,2,2,c,cinv);

        for (index_t it = 0; it < itmax; it++)
        {
            dc33 = -2. * S33 / C3333;
            c(2,2) += dc33;

            GISMO_ASSERT(c(2,2)>= 0,"ERROR! c(2,2) = " << c(2,2) << " C3333=" << C3333 <<" S33=" << S33<<" dc33 = "<<dc33);
            cinv(2,2) = 1.0/c(2,2);

            m_J = m_J0 * math::sqrt( c(2,2) );

            S33     = Sij(2,2,c,cinv);
            C3333   = Cijkl3D(2,2,2,2,c,cinv);

            // gsDebugVar(abs(S33/S33_old));

            if (abs(S33/S33_old) < tol)
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
                else
                    GISMO_ERROR("no vector or matrix produced");

                    break;
            }
            GISMO_ASSERT(it != itmax-1,"Error: Method did not converge, S33 = "<<S33/S33_old<<" and tolerance = "<<tol<<"\n");
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
    T mu = m_par1val / (2. * (1. + m_par2val));

    GISMO_ASSERT( ( (i < 3) && (j < 3) ) , "Index out of range. i="<<i<<", j="<<j);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material == 2)
    {
        return 0.5 * mu * m_Gcon_ori(i,j);
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
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

    if (m_material == 2)
    {
        return 0.0;
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
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

    if (m_material == 2)
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

    if (m_material == 2)
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
T gsMaterialMatrix<T>::dPsi_da(const index_t a) const
{
    T tmp = 0.0;
    if (m_material==4)
    {
        T mu = m_par1val / (2 * (1 + m_par2val));
        tmp  = mu/2 * 2*m_stretches(a);
    }
    else
        GISMO_ERROR("not available...");

    return tmp;
}

template<class T>
T gsMaterialMatrix<T>::d2Psi_dab(const index_t a, const index_t b) const
{
    T tmp = 0.0;
    if (m_material==4)
    {
        T mu = m_par1val / (2 * (1 + m_par2val));
        tmp  = mu/2 * 2 * delta(a,b);
    }
    else
        GISMO_ERROR("not available...");

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
    return ( 1.0 - delta(a,b) ) / ( m_stretches(a) * m_stretches(b) );
}

template<class T>
T gsMaterialMatrix<T>::p() const
{
    return m_stretches(2) * dPsi_da(2);
}

template<class T>
T gsMaterialMatrix<T>::dp_da(const index_t a) const
{
    return m_stretches(2) * d2Psi_dab(3,a) + delta(a,2) * dPsi_da(2);
}

template<class T>
T gsMaterialMatrix<T>::Sa(const index_t a) const
{
    return 1.0/m_stretches(a) * (dPsi_da(a) - p() * dJ_da(a) );
}

template<class T>
T gsMaterialMatrix<T>::dSa_db(const index_t a, const index_t b) const
{
    return 1.0/m_stretches(a) * ( d2Psi_dab(a,b) - dp_da(a)*dJ_da(b) - dp_da(b)*dJ_da(a) - p() * d2J_dab(a,b) ) - delta(a,b) / math::pow(m_stretches(a),2) * (dPsi_da(a) - p() * dJ_da(a));
}

template<class T>
T gsMaterialMatrix<T>::Cabcd(const index_t a, const index_t b, const index_t c, const index_t d) const
{
    GISMO_ASSERT( ( (a < 2) && (b < 2) && (c < 2) && (d < 2) ) , "Index out of range. a="<<a<<", b="<<b<<", c="<<c<<", d="<<d);

    T tmp = 0.0;

    // gsDebug<<"a = "<<a<<"; b = "<<b<<"; m_stretches(a) "<<m_stretches(a)<<"; m_stretches(b) "<<m_stretches(b)<<"\n";

    if ( (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-12) && (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) > 1e-14) )
        gsInfo<<"Warning: difference in stretches is close to machine precision?\n";

    if (abs((m_stretches(a) - m_stretches(b)) / m_stretches(a)) < 1e-14)
        tmp = 1.0 / (2.0 * m_stretches(a) ) * ( dSa_db(b,b) - dSa_db(a,b));
    else
        tmp = ( Sa(b)-Sa(a) ) / (m_stretches(b) - m_stretches(a));

    return 1/m_stretches(c) * dSa_db(a,c) * delta(a,b) * delta(c,d) + tmp * (delta(a,c)*delta(b,d) + delta(a,d)*delta(b,c)) * (1-delta(a,b))
            + delta(a,b)*delta(c,d)*1/(math::pow(m_stretches(a),2) * math::pow(m_stretches(c),2)) * ( math::pow(m_stretches(2),2) * d2Psi_dab(2,2) + 2*dPsi_da(2)*m_stretches(2) );
}


// template<class T>
// T gsMaterialMatrix<T>::dPsi(const index_t a) const
// {
//     T mu = m_par1val / (2. * (1. + m_par2val));

//     GISMO_ASSERT( ( (a < 3) ) , "Index out of range. a="<<a);
//     GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

//     if (m_material==9)
//     {
//         return mu * m_stretches(a,0);
//     }
//     else if (m_material==3)
//     {
//         // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
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
//     T mu = m_par1val / (2. * (1. + m_par2val));

//     GISMO_ASSERT( ( (a < 3) && (b < 3) ) , "Index out of range. a="<<a<<", b="<<b);
//     GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

//     if (m_material==9)
//     {
//         return ( a==b ? mu : 0.0 );
//     }
//     else if (m_material==3)
//     {
//         // return 2.0*m_par1val*m_J0*m_G(i,j)*m_G(k,l) + m_G(i,k)*m_G(j,l) + m_G(i,l)*m_G(j,k);
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

    gsMatrix<T,2,2> tmp;
    tmp.setZero();
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

    gsMatrix<T,3,1> normal = m_map_def.normal(k).normalized();

    m_gcov_def.block(0,0,3,2) = m_acov_def + z * m_ncov_def;
    m_gcov_def.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_def.col(c) = m_Gcon_def(c,0) * m_gcov_def.col(0)
                            + m_Gcon_def(c,1) * m_gcov_def.col(1)
                            + m_Gcon_def(c,2) * m_gcov_def.col(2);
    }

    // // Debugging
    // gsMatrix<T> result(3,3);
    // for (index_t c=0; c!=2; c++)
    //     for (index_t r=0; r!=2; r++)
    //     {
    //         result(r,c) = m_gcon_def.col(r).transpose() * m_gcon_def.col(c);
    //     }
    // gsDebugVar(result);
    // gsDebugVar(m_Gcon_ori);

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

    gsMatrix<T,3,1> normal = m_map.normal(k).normalized();

    m_gcov_ori.block(0,0,3,2) = m_acov_ori + z * m_ncov_ori;
    m_gcov_ori.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_ori.col(c) = m_Gcon_ori(c,0) * m_gcov_ori.col(0)
                            + m_Gcon_ori(c,1) * m_gcov_ori.col(1)
                            + m_Gcon_ori(c,2) * m_gcov_ori.col(2);
    }

    // // Debugging
    // gsMatrix<T> result1(3,3);
    // gsMatrix<T> result2(3,3);
    // for (index_t c=0; c!=2; c++)
    //     for (index_t r=0; r!=2; r++)
    //     {
    //         result1(r,c) = m_gcon_ori.col(r).transpose() * m_gcon_ori.col(c);
    //         result2(r,c) = m_gcov_ori.col(r).transpose() * m_gcov_ori.col(c);
    //         gsDebugVar( z*z*m_ncov_ori.col(r).transpose()*m_ncov_ori.col(c));
    //     }



    // gsDebugVar(result1);
    // gsDebugVar(result2);
    // gsDebugVar(m_Gcon_ori);
    // gsDebugVar(m_Gcov_ori);
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
    // m_stretches.resize(3,1);
    // m_stretchvec.resize(3,3);
    // gsMatrix<T,3,3> evs, evecs;
    // Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

    // eigSolver.compute(C);
    // evs = eigSolver.eigenvalues();
    // evecs = eigSolver.eigenvectors();
    // for (index_t r=0; r!=3; r++)
    //     m_stretches(r,0) = math::sqrt(evs(r,0));

    // m_stretchvec = eigSolver.eigenvectors();

    // gsDebugVar(m_stretchvec);
    // gsDebugVar(m_stretches);


    m_stretches.resize(3,1);
    m_stretchvec.resize(3,3);
    gsMatrix<T,2,1> evs;
    Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

    eigSolver.compute(C.block(0,0,2,2));
    evs = eigSolver.eigenvalues();
    for (index_t r=0; r!=2; r++)
        m_stretches(r,0) = math::sqrt(evs(r,0));
    m_stretches(2,0) = math::sqrt(C(2,2));

    m_stretchvec.setZero();
    m_stretchvec.block(0,0,2,2) = eigSolver.eigenvectors();
    m_stretchvec(2,2) = 1.0;
    // gsDebugVar(m_stretchvec);
    // gsDebugVar(m_stretches);
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