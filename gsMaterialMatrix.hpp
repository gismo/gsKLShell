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
                    "2 : Neo-Hookean Compressible,\n"
                    "3 : Neo-Hookean Incompressible,\n"
                    "4 : Neo-Hookean Compressible,\n"
                    "5 : Neo-Hookean Incompressible,\n"
                    "6 : Neo-Hookean Compressible,\n"
                    "7 : Neo-Hookean Incompressible"
                    ,material_law::SvK_Isotropic);

    m_options.addInt("Compressibility","Specifies whether the material is modelled compressibile or incompressible",compressibility::incompressible);

    m_options.addInt("NumGauss","Number of Gaussian points through thickness",10);
}

template <class T>
short_t gsMaterialMatrix<T>::domainDim() const { return 2; }

template <class T>
short_t gsMaterialMatrix<T>::targetDim() const
{
    if (m_output==2)
        return 9;
    else if (m_output==1)
        return 3;
    else if (m_output==0)
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
    m_output = 2;
}


template <class T>
void gsMaterialMatrix<T>::computePoints(const gsMatrix<T> & u, bool deformed) const
{
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    if (m_defpatches->nPieces()!=0)
    {
        m_map_def.flags = m_map.flags;
        m_map_def.points = u;
        static_cast<const gsFunction<T>&>(m_defpatches->piece(0)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
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

    if ((m_output==1) || (m_output==2)) // for matrix and vector
    {
        this->computePoints(u);

        if (m_material==0)
            result = integrateZ(u);
        else if (m_material==1)
            result = eval_Composite(u); // thickness integration is done inside
        else
            result = integrateZ(u);
    }
    else if (m_output==0)
    {
        m_map.flags = NEED_VALUE;
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
    else 
        GISMO_ERROR("Output type unknown");
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
    else
        gsWarn<<"no function available.";
    // else if ((m_material == 3) && (m_compressible)) // NHK
    //     result = eval_Compressible(i, z);
    // else if ((m_material == 3) && (!m_compressible)) // NHK
    //     result = eval_Incompressible(i, z);


    // switch (m_material)
    // {
    //     case 0 :
    //         break;
    //     case 1 :
    //         gsDebug<<"No eval3D function for composites available. Something went wrong...";
    //         break;
        // case 2 :
        //     result = evalThickness(i, z);
        //     break;
        // case 3 :
        //     result = evalThickness(i, z);
        //     break;
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
gsMatrix<T> gsMaterialMatrix<T>::multiplyZ(const gsMatrix<T>& u) const
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
gsMatrix<T> gsMaterialMatrix<T>::integrateZ(const gsMatrix<T>& u) const
{
    // Input: points in R2
    // Ouput: results in targetDim
    m_numGauss = m_options.getInt("NumGauss");

    gsMatrix<T> result(9,1);
    result.resize(this->targetDim(),u.cols());
    result.setZero();

    m_gauss = gsGaussRule<T>(m_numGauss);
    m_points3D.resize(1,m_numGauss);

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
        m_points3D.bottomRows(1) = quNodes;

        m_evalPoints = this->eval3D(j, m_points3D);
        // gsDebugVar(m_evalPoints);
        for (index_t i=0; i!=this->targetDim(); ++i) // components
        {
            res = 0.0;
            for (index_t k = 0; k != m_numGauss; ++k) // compute integral
                res += quWeights.at(k) * math::pow(quNodes(0,k),m_moment) * m_evalPoints(i,k);
            result(i,j) = res;
        }

    }
    return result;
}

template<class T>
void gsMaterialMatrix<T>::computeMetric(index_t k, T z, bool computedeformed, bool computefull) const
{
    // k: index of quadrature point
    // computedeformed: compute quantities also on deformed shape
    // computeall: compute metric of full coordinate system rather than only in-plane basis

    m_metricA.resize(2,2);
    // Construct metric tensor a = [dcd1*dcd1, dcd1*dcd2; dcd2*dcd1, dcd2*dcd2]
    m_metricA        = m_map.jacobian(k);
    m_metricA        = m_metricA.transpose() * m_metricA;
    if (computefull)
    {
        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        m_deriv2    = m_map.deriv2(k);
        m_deriv2.resize(3,3);
        m_normal    = m_map.normal(k).normalized();

        m_metricB.resize(2,2);
        m_metricB(0,0) = m_deriv2.row(0).dot(m_normal);
        m_metricB(1,1) = m_deriv2.row(1).dot(m_normal);
        m_metricB(0,1) = m_metricB(1,0) = m_deriv2.row(2).dot(m_normal);
        // Construct metric of coordinate system g = [a_ij - 2*theta3*b_ij]

        m_metricG.resize(3,3);
        m_metricG.setZero();
        m_metricG.block(0,0,2,2)        = m_metricA - 2 * z * m_metricB;
        m_metricG(2,2) = 1.0;
    }
    if (computedeformed)
    {
        GISMO_ASSERT(m_defpatches->nPieces()!=0,"Deformed multipatch is empty; cannot compute strains!");
        m_metricA_def.resize(2,2);
        // Construct metric tensor a = [dcd1*dcd1, dcd1*dcd2; dcd2*dcd1, dcd2*dcd2]
        m_metricA_def    = m_map_def.jacobian(k);
        m_metricA_def    = m_metricA_def.transpose() * m_metricA_def;
        if (computefull)
        {
            // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
            m_deriv2_def = m_map_def.deriv2(k);
            m_deriv2_def.resize(3,3); // gives [d11 c1, d11c2, d11c3; d22c1, d22c1, d22c3; d12c1, d12c2, d12c3]
            m_normal_def = m_map_def.normal(k).normalized();

            m_metricB_def.resize(2,2);
            m_metricB_def(0,0) = m_deriv2_def.row(0).dot(m_normal_def);
            m_metricB_def(1,1) = m_deriv2_def.row(1).dot(m_normal_def);
            m_metricB_def(0,1) = m_metricB_def(1,0) = m_deriv2_def.row(2).dot(m_normal_def);
            // Construct metric of coordinate system g = [a_ij - 2*theta3*b_ij]
            m_metricG_def = m_metricA_def - 2 * z * m_metricB_def;

            m_metricG_def.resize(3,3);
            m_metricG_def.setZero();
            m_metricG_def.block(0,0,2,2)= m_metricA_def - 2 * z * m_metricB_def;
            m_metricG_def(2,2) = 1.0;

            m_J0 = math::sqrt( m_metricG_def.determinant() / m_metricG.determinant() );
            m_J0 = math::pow( m_J0, -2 );
        }
    }
}

// template <class T>
// gsMatrix<T> gsMaterialMatrix<T>::eval3D_Linear(const gsMatrix<T>& u, const gsMatrix<T>& z) const
// {
//     // gsInfo<<"TO DO: evaluate moments using thickness";
//     // Input: in-plane points in R2 (u)
//     //        out-of-plane points (through thickness) in R1 (z)
//     // Output: material matrix in R9 (cols stacked in rows; every row corresponds to a point)
//     gsMatrix<T> result(this->targetDim(), u.cols() * z.cols());
//     result.setZero();

//     index_t sz = z.cols();
//     for( index_t i=0; i < u.cols(); ++i )
//         {
//             gsDebugVar(u.col(i));
//             for( index_t j=0; j < z.cols(); ++j )
//             {
//                 // Evaluate material properties on the quadrature point
//                 m_par1val = m_par1mat(0,i);
//                 m_par2val = m_par2mat(0,i);

//                 if (m_matrix)
//                 {
//                     this->computeMetric(i,z(0,j)); // on point i, on height z(0,j)

//                     gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j+i*sz,3,3);
//                     /*
//                         C = C1111,  C1122,  C1112
//                             symm,   C2222,  C2212
//                             symm,   symm,   C1212
//                     */
//                     C(0,0)          = Cijkl(0,0,0,0); // C1111
//                     C(1,1)          = Cijkl(1,1,1,1); // C2222
//                     C(2,2)          = Cijkl(0,1,0,1); // C1212
//                     C(1,0) = C(0,1) = Cijkl(0,0,1,1); // C1122
//                     C(2,0) = C(0,2) = Cijkl(0,0,0,1); // C1112
//                     C(2,1) = C(1,2) = Cijkl(1,1,0,1); // C2212
//                 }
//                 else
//                 {
//                     this->computeMetric(i,z(0,j),true,true);

//                     gsDebugVar(m_metricA_def);
//                     // gsDebugVar(m_metricB_def);

//                     result(0,j+i*sz) = Sij(0,0);
//                     result(1,j+i*sz) = Sij(1,1);
//                     result(2,j+i*sz) = Sij(0,1);
//                 }
//             }
//             gsDebugVar(result);
//         }

//     return result;
// }

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

        if (m_output==2)
        {
            this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)

            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j,3,3);
            /*
                C = C1111,  C1122,  C1112
                    symm,   C2222,  C2212
                    symm,   symm,   C1212
            */
            C(0,0)          = Cijkl_i(0,0,0,0); // C1111
            C(1,1)          = Cijkl_i(1,1,1,1); // C2222
            C(2,2)          = Cijkl_i(0,1,0,1); // C1212
            C(1,0) = C(0,1) = Cijkl_i(0,0,1,1); // C1122
            C(2,0) = C(0,2) = Cijkl_i(0,0,0,1); // C1112
            C(2,1) = C(1,2) = Cijkl_i(1,1,0,1); // C2212
        }
        else if (m_output==1)
        {
            this->computeMetric(i,z.at(j),true,true);
            result(0,j) = Sij_i(0,0);
            result(1,j) = Sij_i(1,1);
            result(2,j) = Sij_i(0,1);
        }
        else 
            GISMO_ERROR("no vector or matrix produced");
    }
    return result;
}

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval_Composite(const gsMatrix<T>& u) const
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

        switch (m_moment)
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

        if (m_output==2)
            result.reshapeCol(k,3,3) = Cmat;
        else if (m_output==1)
        {
            gsMatrix<T> Eij(2,2);
            computeMetric(k,0.0,true,true); // height is set to 0, but we do not need m_metricG

            if (m_moment==0)
                Eij = 0.5*(m_metricA_def - m_metricA);
            else if (m_moment==2)
                Eij = (m_metricB - m_metricB_def);
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
// template<class T>
// T gsMaterialMatrix<T>::Cijkl(const index_t i, const index_t j, const index_t k, const index_t l) const { Cijkl(i,j,k,l,NULL,NULL); }

template<class T>
T gsMaterialMatrix<T>::Cijkl_i(const index_t i, const index_t j, const index_t k, const index_t l) const
{
    GISMO_ASSERT( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(!m_compressible,"Material model is not incompressible?");

    if (m_material==0) // svk
    {
        T lambda, mu, Cconstant;

        lambda = m_par1val * m_par2val / ( (1. + m_par2val)*(1.-2.*m_par2val)) ;
        mu     = m_par1val / (2.*(1. + m_par2val)) ;

        Cconstant = 2*lambda*mu/(lambda+2*mu);

        return Cconstant*m_metricA(i,j)*m_metricA(k,l) + mu*(m_metricA(i,k)*m_metricA(j,l) + m_metricA(i,l)*m_metricA(j,k));
    }
    else if (m_material==1)
        gsWarn<<"Compressible material matrix  requested, but not needed. How?";
    else if (m_material==2)
    {
        return 2.0*m_par1val*m_J0*m_metricG(i,j)*m_metricG(k,l) + m_metricG(i,k)*m_metricG(j,l) + m_metricG(i,l)*m_metricG(j,k);
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_metricG(i,j)*m_metricG(k,l) + m_metricG(i,k)*m_metricG(j,l) + m_metricG(i,l)*m_metricG(j,k);
    }
    else
        GISMO_ERROR("Material model not implemented (Cijkl_i.");
}

template<class T>
T gsMaterialMatrix<T>::Cijkl_c(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    GISMO_ASSERT(c.cols()==c.rows(),"Matrix c must be square");
    GISMO_ASSERT(c.cols()==3,"Matrix c must be 3x3");
    GISMO_ASSERT(cinv.cols()==cinv.rows(),"Matrix cinv must be square");
    GISMO_ASSERT(cinv.cols()==3,"Matrix cinv must be 3x3");
    GISMO_ASSERT( ( (i <=2) && (j <=2) && (k <=2) && (l <=2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
    GISMO_ASSERT(m_compressible,"Material model is not compressible?");

    if (m_material==0 || m_material==1) // svk
        gsWarn<<"Compressible material matrix requested, but not needed. How?";
    else if (m_material==2)
    {
        return 1.0 / 9.0 * m_par1val * math::pow( m_J , -2.0/3.0 ) * ( c.trace() * ( 2*cinv(i,j)*cinv(k,l) + 3*cinv(i,k)*cinv(j,l) + 3*cinv(i,l)*cinv(j,k) )
                            - 6*m_metricG_def(i,j)*cinv(k,l) + cinv(i,j)*m_metricG(k,l) ) + m_par2val * ( m_J*m_J*cinv(i,j)*cinv(k,l) - 0.5*(m_J*m_J-1)*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) ) );
    }
    else if (m_material==3)
    {
        // return 2.0*m_par1val*m_J0*m_metricG(i,j)*m_metricG(k,l) + m_metricG(i,k)*m_metricG(j,l) + m_metricG(i,l)*m_metricG(j,k);
    }
    else
        GISMO_ERROR("Material model not implemented (Cijkl_c.");

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
T gsMaterialMatrix<T>::Sij_i(const index_t i, const index_t j) const
{
    gsMatrix<T> tmp;
    if (m_material==0)
    {
        if (m_moment==0)
        {
            GISMO_ASSERT( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
            tmp = 0.5*(m_metricA_def - m_metricA);
        }
        else if (m_moment==2)
        {
            GISMO_ASSERT( ( (i < 2) && (j < 2) ) , "Index out of range. i="<<i<<", j="<<j);
            tmp = (m_metricB - m_metricB_def);
        }
        else
        {
            gsDebug<<"Warning: no material model known in simplification";
        }
        return    Cijkl_i(i,j,0,0) * tmp(0,0) + Cijkl_i(i,j,0,1) * tmp(0,1)
                + Cijkl_i(i,j,1,0) * tmp(1,0) + Cijkl_i(i,j,1,1) * tmp(1,1);
    }
    else if (m_material==1)
        gsWarn<<"Incompressible material stress tensor requested, but not needed. How?";
    else if (m_material==2)
    {
            return m_par1val * (m_metricG(i,j) - math::pow(m_J0,-2.) * m_metricG_def(i,j) );
    }
    else
        GISMO_ERROR("Material model not implemented (Sij_i).");
}

template<class T>
T gsMaterialMatrix<T>::Sij_c(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
{
    T tmp;
    if (m_material==0 || m_material==1)
        gsWarn<<"Incompressible material stress tensor requested, but not needed. How?";
    else if (m_material==2)
    {
            gsDebugVar(m_par1val);
            gsDebugVar(m_J);
            gsDebugVar(m_metricG(i,j));
            gsDebugVar(c);
            gsDebugVar(cinv);
            gsDebugVar(m_par2val);


            tmp =  m_par1val * math::pow( m_J , -2.0/3.0 ) * ( m_metricG(i,j) - 1.0/3.0 * c.trace() * cinv(i,j) ) + 0.5 * m_par2val * ( m_J*m_J - 1 ) * cinv(i,j);
            return tmp;
    }
    else
        GISMO_ERROR("Material model not implemented (Sij_i).");
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

        this->computeMetric(i,z.at(j),true,true); // on point i, on height z(0,j)

        // Define objects
        gsMatrix<T,3,3> c, cinv;
        T S33, C3333, dc33, traceC;

        // Initialize c
        c.setZero();
        c.block(0,0,2,2) = m_metricG_def.block(0,0,2,2);
        c(2,2) = 1.0; // c33
        cinv = c.inverse();
        // note: can also just do c = jacGdef because the normal has length one and hence c(2,2) is 1. CHECK!

        index_t itmax = 2;
        T tol = 1e-6;
        S33 = 0.0;
        C3333 = 1.0;

        for (index_t it = 0; it < itmax; it++)
        {
            dc33 = -2. * S33 / C3333;
            c(2,2) += dc33;
            cinv(2,2) = 1.0/c(2,2);

            traceC = c.trace();
            m_J = m_J0 * math::sqrt( c(2,2) );

            S33     = Sij_c(2,2,c,cinv);
            C3333   = Cijkl_c(2,2,2,2,c,cinv);
            if (S33 < tol)
            {
                // gsInfo<<"Converged in "<<it<<" iterations, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";
                if (m_output==2)
                    {
                        /*
                            C =     C1111,  C1122,  C1112
                                    symm,   C2222,  C2212
                                    symm,   symm,   C1212
                            Here, Cabcd = Cijkl - Cab33*C33cd / C3333;
                            a,b,c,d = 1,2; i,j,k,l = 1...3;
                        */
                        gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(j,3,3);
                        C(0,0) = Cijkl_c(0,0,0,0,c,cinv) - ( Cijkl_c(0,0,2,2,c,cinv) * Cijkl_c(2,2,0,0,c,cinv) ) / (Cijkl_c(2,2,2,2,c,cinv)); // C1111
                        C(0,1) =
                        C(1,0) = Cijkl_c(0,0,1,1,c,cinv) - ( Cijkl_c(0,0,2,2,c,cinv) * Cijkl_c(2,2,1,1,c,cinv) ) / (Cijkl_c(2,2,2,2,c,cinv)); // C1122
                        C(0,2) =
                        C(2,0) = Cijkl_c(0,0,0,1,c,cinv) - ( Cijkl_c(0,0,2,2,c,cinv) * Cijkl_c(2,2,0,1,c,cinv) ) / (Cijkl_c(2,2,2,2,c,cinv)); // C1112
                        C(1,1) = Cijkl_c(1,1,1,1,c,cinv) - ( Cijkl_c(1,1,2,2,c,cinv) * Cijkl_c(2,2,1,1,c,cinv) ) / (Cijkl_c(2,2,2,2,c,cinv)); // C2222
                        C(1,2) =
                        C(2,1) = Cijkl_c(1,1,0,1,c,cinv) - ( Cijkl_c(1,1,2,2,c,cinv) * Cijkl_c(2,2,0,1,c,cinv) ) / (Cijkl_c(2,2,2,2,c,cinv)); // C2212
                        C(2,2) = Cijkl_c(0,1,0,1,c,cinv) - ( Cijkl_c(0,1,2,2,c,cinv) * Cijkl_c(2,2,0,1,c,cinv) ) / (Cijkl_c(2,2,2,2,c,cinv)); // C1212
                    }
                    else if (m_output==1)
                    {
                        result(0,j) = Sij_c(0,0,c,cinv); // S11
                        result(1,j) = Sij_c(1,1,c,cinv); // S22
                        result(2,j) = Sij_c(0,1,c,cinv); // S12
                    }
                    else
                        GISMO_ERROR("no vector or matrix produced");
                }
            else if (it == itmax - 1)
            {
                gsInfo<<"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";
                // std::terminate();
            }
        }
        return result;
    }
}

} // end namespace