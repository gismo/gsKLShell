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
    - Finish Compressible and Incompressible
    - Add option to assemble stress tensor. Note: this requires computation of the strain as well!!!!!!!
    - Add option to the class for different material models, to avoid making lambda functions. Note; use references to m_defG etc
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
    m_model = 0;
    m_moment = 0;
    m_matrix = true;
    m_numGauss = 10;

    m_map.flags = m_map_def.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
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
    m_model = 0;
    m_moment = 0;
    m_matrix = true;
    m_numGauss = 10;

    m_map.flags = m_map_def.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
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
    m_model = 1;
    m_moment = 0;
    m_matrix = true;
    m_numGauss = 10;

    m_map.flags = m_map_def.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
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
//     // m_model = 0;
//     m_moment = 0;
//     m_map.flags = m_map_def.flags =  NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;

//     // normal and deriv2 only needed when ComputeFull
// }


template <class T>
gsOptionList gsMaterialMatrix<T>::defaultOptions()
{
    gsOptionList opt;
    // to do
    // opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",material_law::saint_venant_kirchhoff);
    return opt;
}

template <class T>
short_t gsMaterialMatrix<T>::domainDim() const { return 2; }

template <class T>
short_t gsMaterialMatrix<T>::targetDim() const
{
    if (m_matrix)
        return 9;
    else
        return 3;
}

// template <class T>
// void gsMaterialMatrix<T>::computePoints(const gsMatrix<T> & u) const
// {
//     m_map.points = m_map_def.points = u;

//     // NOTE 1: if the input \a u is considered to be in physical coordinates
//     // then we first need to invert the points to parameter space
//     // m_patches.patch(0).invertPoints(u, m_map.points, 1e-8) which is not exact (!),
//     // otherwise we just use the input paramteric points

//     static_cast<const gsFunction<T>&>(m_patches->piece(0)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
//     static_cast<const gsFunction<T>&>(m_patches->piece(0)).computeMap(m_map_def);

//     // NOTE 2: in the case that parametric value is needed it suffices
//     // to evaluate Youngs modulus and Poisson's ratio at
//     // \a u instead of _tmp.values[0].

//     m_thickness->eval_into(m_map.values[0], m_Tmat);



// }


template <class T>
void gsMaterialMatrix<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_map.points = u;
    m_map_def.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    static_cast<const gsFunction<T>&>(m_defpatches->piece(0)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    // Compute the thickness
    m_thickness->eval_into(m_map.values[0], m_Tmat);
    m_par1->eval_into(m_map.values[0], m_par1mat);
    m_par2->eval_into(m_map.values[0], m_par2mat);

    switch (m_model)
    {
        case 0 : // linear model
            // result = multiplyZ(u); // the matrix does not depend on thickness so we can multiply
            result = integrateZ(u); // the matrix does not depend on thickness so we can multiply
            break;
        case 1 : // composite model
            result = eval_Composite(u); // thickness integration is done inside
            break;

        // case 2 : eval3D_Incompressible(u,result);
        //     break;
        // case 3 : eval3D_Compressible(u,result);
        //     break;
    }
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
    switch (m_model)
    {
        case 0 :
            result = evalThickness(i, z);
            break;
        case 1 :
            gsDebug<<"No eval3D function for composites available. Something went wrong...";
            break;

        // case 1 : eval3D_Incompressible(u,result);
        //     break;
        // case 2 : eval3D_Incompressible(u,result);
        //     break;
    }
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
        // auxilary variable since eval
        gsMatrix<T> z(1,1);
        z.setZero();

        T fac;
        for (index_t i = 0; i != u.cols(); ++i) // points
        {
            m_evalPoints = this->evalThickness(i,z);
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
        m_metricG     = m_metricA - 2 * z * m_metricB;
    }
    if (computedeformed)
    {
        GISMO_ASSERT(NULL!=m_defpatches->size(),"Deformed multipatch is empty; cannot compute strains!");
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
gsMatrix<T> gsMaterialMatrix<T>::evalThickness(const index_t i, const gsMatrix<T>& z) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: j index in-plane point
    //        z out-of-plane coordinate (through thickness) in R1 (z)
    gsMatrix<T> result(this->targetDim(), z.cols());
    result.setZero();

    for( index_t j=0; j < z.cols(); ++j )
    {
        // Evaluate material properties on the quadrature point
        m_par1val = m_par1mat(0,i);
        m_par2val = m_par2mat(0,i);

        if (m_matrix)
        {
            this->computeMetric(i,z.at(j)); // on point i, on height z(0,j)

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
        else
        {
            this->computeMetric(i,z.at(j),true,true);
            result(0,j) = Sij(0,0);
            result(1,j) = Sij(1,1);
            result(2,j) = Sij(0,1);
        }
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

    m_map.points = m_map_def.points = u.topRows(2);
    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    static_cast<const gsFunction<T>&>(m_defpatches->piece(m_pIndex)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

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

        if (m_matrix)
            result.reshapeCol(k,3,3) = Cmat;
        else
        {
            gsMatrix<T> Eij(2,2);
            computeMetric(k,0.0,true,true); // height is set to 0, but we do not need m_metricG

            if (m_moment==0)
            {
                Eij = 0.5*(m_metricA - m_metricA_def);
                gsDebugVar(Eij);
            }
            else if (m_moment==2)
            {
                Eij = (m_metricB - m_metricB_def);
                gsDebugVar(Eij);
            }
            else
            {
                gsDebug<<"Warning: no material model known in simplification";
            }
            Eij(0,1) += Eij(1,0);
            std::swap(Eij(1,0), Eij(1,1));
            Eij.resize(4,1);
            Eij.conservativeResize(3,1);
            result = Cmat * Eij;
        }







        // else
        // {
        //     // result.col(k) = Cmat *
        // }
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
    switch (m_model)
    {
        case 0: // linear St. Venant-Kirchhoff
            GISMO_ASSERT( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);

            T lambda, mu, Cconstant;

            lambda = m_par1val * m_par2val / ( (1. + m_par2val)*(1.-2.*m_par2val)) ;
            mu     = m_par1val / (2.*(1. + m_par2val)) ;

            Cconstant = 2*lambda*mu/(lambda+2*mu);

            return Cconstant*m_metricA(i,j)*m_metricA(k,l) + mu*(m_metricA(i,k)*m_metricA(j,l) + m_metricA(i,l)*m_metricA(j,k));

            break;

        case 1: // composite material model
            // NO IMPLEMENTATION
            break;

        case 2: // nonlinear Incompressible Neo-Hookean
            GISMO_ASSERT( ( (i < 2) && (j < 2) && (k < 2) && (l < 2) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);
            return 2.0*m_par1val*m_J0*m_metricG(i,j)*m_metricG(k,l) + m_metricG(i,k)*m_metricG(j,l) + m_metricG(i,l)*m_metricG(j,k);
            break;

        // case 2: // nonlinear Compressible Neo-Hookean
        //     GISMO_ASSERT( ( (i < 3) && (j < 3) && (k < 3) && (l < 3) ) , "Index out of range. i="<<i<<", j="<<j<<", k="<<k<<", l="<<l);

        //     return 2*m_par1*m_J0*m_metric(i,j)*m_metric(k,l) + m_metric(i,k)*m_metric(j,l) + m_metric(i,l)*m_metric(j,k);
        //     break;
    }

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
template<class T>
T gsMaterialMatrix<T>::Sij(const index_t i, const index_t j) const
{

    // CHECK IF QUANTITIES HAVE BEEN COMPUTED!

    gsMatrix<T> tmp;
    switch (m_model)
    {
        case 0: // linear St. Venant-Kirchhoff
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
            return    Cijkl(i,j,0,0) * tmp(0,0) + Cijkl(i,j,0,1) * tmp(0,1)
                    + Cijkl(i,j,1,0) * tmp(1,0) + Cijkl(i,j,1,1) * tmp(1,1);
            break;

        // composite material model
        // Note: Provides the strain tensor.
        case 1:
            // NO IMPLEMENTATION
            break;

        case 2: // nonlinear Incompressible Neo-Hookean
            return m_par1val * (m_metricG(i,j) - math::pow(m_J0,-2.) * m_metricG_def(i,j) );
            break;
    }
}

template<class T>
gsMatrix<T> gsMaterialMatrix<T>::S() const
{

    // CHECK IF QUANTITIES HAVE BEEN COMPUTED!

    gsMatrix<T> tmp;
    switch (m_model)
    {
        case 0: // linear St. Venant-Kirchhoff

            break;

        // composite material model
        // Note: Provides the strain tensor.
        case 1:
            if (m_moment==0)
            {
                tmp = 0.5*(m_metricA - m_metricA_def);
            }
            else if (m_moment==2)
            {
                tmp = (m_metricB - m_metricB_def);
            }
            else
            {
                gsDebug<<"Warning: no material model known in simplification";
            }
            return tmp;
            break;

        case 2: // nonlinear Incompressible Neo-Hookean

            break;
    }
}

template<class T>
gsMatrix<T> gsMaterialMatrix<T>::eval3D_Incompressible(const gsMatrix<T>& u) const
{
    m_result.resize(this->targetDim(),u.cols());

    // NOTE 1: if the input \a u is considered to be in physical coordinates
    // then we first need to invert the points to parameter space
    // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
    // otherwise we just use the input paramteric points
    m_map.points     = u.topRows(2);
    m_map_def.points = u.topRows(2);

    static_cast<const gsFunction<T>&>( m_patches->piece(m_pIndex)     ).computeMap(m_map);
    static_cast<const gsFunction<T>&>( m_defpatches->piece(m_pIndex)  ).computeMap(m_map_def);

    // NOTE 2: in the case that parametric value is needed it suffices
    // to evaluate Youngs modulus and Poisson's ratio at
    // \a u instead of _tmp.values[0].
    m_par1->eval_into(m_map.values[0], m_par1mat);
    // m_par2->eval_into(m_map.values[0], m_par2mat);

    m_result.resize( targetDim() , u.cols() );
    for( index_t k=0; k< u.cols(); ++k )
    {
        this->computeMetric(k, u(2,k), true, true);

        // Evaluate material properties on the quadrature point
        m_par1val = m_par1mat(0,k);
        m_J0 = math::sqrt( m_metricG_def.determinant() / m_metricG.determinant() );
        m_J0 = math::pow( m_J0, -2 );

        if (m_matrix)
        {
            /*
                C =     C1111,  C1122,  C1112
                        symm,   C2222,  C2212
                        symm,   symm,   C1212
            */
            gsAsMatrix<T, Dynamic, Dynamic> C = m_result.reshapeCol(k,3,3);
            C(0,0) = Cijkl(0,0,0,0); // C1111
            C(1,0) = C(0,1) = Cijkl(0,0,1,1); // C1122
            C(2,0) = C(0,2) = Cijkl(0,0,0,1); // C1112
            C(1,1) = Cijkl(1,1,1,1); // C2222
            C(2,1) = C(1,2) = Cijkl(1,1,0,1); // C2212
            C(2,2) = Cijkl(0,1,0,1); // C1212
            break;
        }
        else
        {
            m_result(0,k) = Sij(0,0);
            m_result(1,k) = Sij(1,1);
            m_result(2,k) = Sij(0,1);
            break;
        }
    }
    return m_result;
};

// void eval3D_Compressible(const gsMatrix<T>& u, gsMatrix<T>& result) const
// {
//     // NOTE 1: if the input \a u is considered to be in physical coordinates
//     // then we first need to invert the points to parameter space
//     // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
//     // otherwise we just use the input paramteric points
//     _tmp.points     = u.topRows(2);
//     _tmp_def.points = u.topRows(2);

//     static_cast<const gsFunction<T>&>( _mp->piece(0)     ).computeMap(_tmp);
//     static_cast<const gsFunction<T>&>( _mp_def->piece(0) ).computeMap(_tmp_def);

//     // NOTE 2: in the case that parametric value is needed it suffices
//     // to evaluate Youngs modulus and Poisson's ratio at
//     // \a u instead of _tmp.values[0].
//     _par1->eval_into(_tmp.values[0], par1mat);
//     _par2->eval_into(_tmp.values[0], par2mat);

//     result.resize( targetDim() , u.cols() );
//     for( index_t k=0; k< u.cols(); ++k )
//     {
//         // Material parameters
//         mu = par1mat(0,k);
//         K = par1mat(0,k);

//         // Define objects
//         gsMatrix<T,3,3> c, cinv;
//         T S33, C3333, dc33, traceC;

//         // Construct metric tensor a = [dcd1*dcd1, dcd1*dcd2; dcd2*dcd1, dcd2*dcd2]
//         jacGdef = _tmp_def.jacobian(k);
//         a_def   = jacGdef.transpose() * jacGdef;
//         jacGori = _tmp.jacobian(k);
//         a_ori   = jacGori.transpose() * jacGori;

//         // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
//         gsAsConstMatrix<T,3,3> deriv2def( _tmp_def.deriv2(0).data(),3,3 ); // gives [d11 c1, d11c2, d11c3; d22c1, d22c1, d22c3; d12c1, d12c2, d12c3]
//         gsAsConstMatrix<T,3,3> deriv2ori(     _tmp.deriv2(0).data(),3,3 );
//         n_def = _tmp_def.normal(k).normalized();
//         n_ori = _tmp.normal(k).normalized();
//         b_def.resize(2,2);
//         b_ori.resize(2,2);

//         b_def(0,0) = deriv2def.row(0).dot(n_def);
//         b_def(1,1) = deriv2def.row(1).dot(n_def);
//         b_def(0,1) = b_def(1,0) = deriv2def.row(2).dot(n_def);

//         b_ori(0,0) = deriv2ori.row(0).dot(n_ori);
//         b_ori(1,1) = deriv2ori.row(1).dot(n_ori);
//         b_ori(0,1) = b_ori(1,0) = deriv2ori.row(2).dot(n_ori);

//         // Construct basis of coordinate system g = [a_ij - 2*theta3*b_ij]
//         g_def = g_ori = gsMatrix<T>::Zero(3,3);
//         g_def.block(0,0,2,2) = a_def - 2 * u(2,k) * b_def;
//         g_ori.block(0,0,2,2) = a_ori - 2 * u(2,k) * b_ori;
//         g_def(2,2) = g_ori(2,2) = 1.0;

//         // Initialize c
//         c.setZero();
//         c.block(0,0,2,2) = g_def.block(0,0,2,2);
//         c(2,2) = 1.0; // c33
//         cinv = c.inverse();
//         // note: can also just do c = jacGdef because the normal has length one and hence c(2,2) is 1. CHECK!

//         J0 = math::sqrt( g_def.determinant() / g_ori.determinant() );
//         J = J0 * math::sqrt( c(2,2) );

//         index_t itmax = 20;
//         T tol = 1e-6;
//         S33 = 0.0;
//         C3333 = 1.0;

//         // Define lambda function for C
//         std::function<T (index_t i, index_t j, index_t k, index_t l)> Cijkl;
//         Cijkl = [=](index_t i, index_t j, index_t k, index_t l)
//         {
//             T res = 1.0 / 9.0 * mu * math::pow( J , -2.0/3.0 ) * ( traceC * ( 2*cinv(i,j)*cinv(k,l) + 3*cinv(i,k)*cinv(j,l) + 3*cinv(i,l)*cinv(j,k) )
//                             - 6*g_ori(i,j)*cinv(k,l) + cinv(i,j)*g_ori(k,l) ) + K * ( J*J*cinv(i,j)*cinv(k,l) - 0.5*(J*J-1)*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) ) );
//             return res;
//         };

//         // Define lambda function for S
//         std::function<T (index_t i, index_t j)> Sij;
//         Sij = [=](index_t i, index_t j)
//         {
//             T res = mu * math::pow( J , -2.0/3.0 ) * ( g_ori(i,j) - 1.0/3.0 * traceC * cinv(i,j) ) + 0.5 * K * ( J*J - 1 ) * cinv(i,j)
//             return res;
//         };

//         for (index_t it = 0; it < itmax; it++)
//         {
//             dc33 = -2. * S33 / C3333;
//             c(2,2) += dc33;
//             cinv(2,2) = 1.0/c(2,2);

//             traceC = c.trace();
//             J = J0 * math::sqrt( c(2,2) );

//             S33     = Sij(2,2);
//             C3333   = Cijkl(2,2,2,2);

//             if (S33 < tol)
//             {
//                 // gsInfo<<"Converged in "<<it<<" iterations, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";
//                 switch (mat)
//                 {
//                     case 0:
//                         result(0,k) = Sij(0,0); // S11
//                         result(1,k) = Sij(1,1); // S22
//                         result(2,k) = Sij(0,1); // S12
//                         break;
//                     case 1:
//                         /*
//                             C =     C1111,  C1122,  C1112
//                                     symm,   C2222,  C2212
//                                     symm,   symm,   C1212
//                             Here, Cabcd = Cijkl - Cab33*C33cd / C3333;
//                             a,b,c,d = 1,2; i,j,k,l = 1...3;
//                         */
//                         gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(k,3,3);
//                         C(0,0) = Cijkl(0,0,0,0) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,0) ) / (Cijkl(2,2,2,2)); // C1111
//                         C(0,1) =
//                         C(1,0) = Cijkl(0,0,1,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C1122
//                         C(0,2) =
//                         C(2,0) = Cijkl(0,0,0,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1112
//                         C(1,1) = Cijkl(1,1,1,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C2222
//                         C(1,2) =
//                         C(2,1) = Cijkl(1,1,0,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C2212
//                         C(2,2) = Cijkl(0,1,0,1) - ( Cijkl(0,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1212
//                         break;
//                 }
//                 break;
//             }
//             else if (it == itmax - 1)
//             {
//                 gsInfo<<"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";
//                 // std::terminate();
//             }
//         }
//     }
// }

} // end namespace