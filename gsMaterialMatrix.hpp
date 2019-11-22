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

#pragma once

#include <gsThinShell2/gsMaterialMatrix.h>

namespace gismo
{

// Linear material models
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & YoungsModulus,
                                        const gsFunction<T> & PoissonRatio
                                        )
                                        :
                                        m_patches(&mp),
                                        m_thickness(&thickness),
                                        m_YoungsModulus(&YoungsModulus),
                                        m_PoissonRatio(&PoissonRatio),
                                        m_piece(nullptr)
{
    m_model = 0;
    m_moment = 0;
    m_map.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
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
    m_map.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
}

// Nonlinear material models
template<class T>
gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
                                        const gsFunctionSet<T> & mp_def,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & YoungsModulus,
                                        const gsFunction<T> & PoissonRatio
                                        )
                                        :
                                        m_patches(&mp),
                                        m_defpatches(&mp_def),
                                        m_thickness(&thickness),
                                        m_YoungsModulus(&YoungsModulus),
                                        m_PoissonRatio(&PoissonRatio),
                                        m_piece(nullptr)
{
    // m_model = 0;
    m_moment = 0;
    m_map.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
}


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
short_t gsMaterialMatrix<T>::targetDim() const { return 9; }

template <class T>
void gsMaterialMatrix<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    switch (m_model)
    {
        case 0 : // linear model
            result = multiplyZ(u, m_moment); // the matrix does not depend on thickness so we can multiply
            break;
        case 1 : // composite model
            result = eval_Composite(u, m_moment); // thickness integration is done inside
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
void gsMaterialMatrix<T>::eval3D_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    switch (m_model)
    {
        case 0 :
            result = eval3D_Linear(u);
            break;
        case 1 :
            gsDebug<<"No eval3D function for composites available. Something went wrong...";
            break;

        // case 1 : eval3D_Incompressible(u,result);
        //     break;
        // case 2 : eval3D_Incompressible(u,result);
        //     break;
    }
}

// multiplies the material matrix by the thickness on all points
template <class T>
gsMatrix<T> gsMaterialMatrix<T>::multiplyZ(const gsMatrix<T>& u, int moment) const
{
    // Input: points in R2
    // Ouput: results in targetDim
    gsMatrix<T> result(this->targetDim(),u.cols());

    if (m_moment % 2 != 0)  //then the moment is odd
        result.setZero();
    else                    // then the moment is even
    {
        m_map.points = u;
        static_cast<const gsFunction<T>&>(m_patches->piece(0)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

        // Compute the thickness
        m_thickness->eval_into(m_map.values[0], m_Tmat);

        // Define evaluation points
        m_points.resize(u.rows() + 1, u.cols() );
        m_points.topRows(u.rows()) = u;
        m_points.bottomRows(1).setZero(); // we dont use this, since the function we evaluate is by choice not dependent on thickness
        this->eval3D_into(m_points,result);

        T fac;
        for (index_t j = 0; j != u.cols(); ++j) // points
        {
            fac = 2.0/(m_moment+1) * math::pow( m_Tmat(0,j) / 2.0 , m_moment + 1);
            result.col(j) *= fac;
        }

    }
    return result;
}

// integrates the material matrix over its last coordinate (thickness)
template <class T>
gsMatrix<T> gsMaterialMatrix<T>::integrateZ(const gsMatrix<T>& u, int moment) const
{
    // Input: points in R2
    // Ouput: results in targetDim
    gsMatrix<T> result(9,1);

    // NOTE 1: if the input \a u is considered to be in physical coordinates
    // then we first need to invert the points to parameter space
    // m_patches.patch(0).invertPoints(u, m_map.points, 1e-8) which is not exact (!),
    // otherwise we just use the input paramteric points
    m_map.points = u;

    static_cast<const gsFunction<T>&>(m_patches->piece(0)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    // NOTE 2: in the case that parametric value is needed it suffices
    // to evaluate Youngs modulus and Poisson's ratio at
    // \a u instead of _tmp.values[0].
    m_thickness->eval_into(m_map.values[0], m_Tmat);

    result.resize(this->targetDim(),u.cols());

    m_numGauss = 10;
    m_gauss = gsGaussRule<T>(m_numGauss);
    m_points = gsMatrix<T>(u.rows() + 1,m_numGauss);

    T res;
    for (index_t j = 0; j != u.cols(); ++j) // points
    {
        // Compute values of in-plane points
        m_points.topRows(u.rows()) = u.col(j).replicate(1,m_numGauss);

        // set new integration point
        m_tHalf = m_Tmat(0,j)/2.0;
        m_gauss.mapTo(-m_tHalf,m_tHalf,m_quNodes,m_quWeights);
        m_points.bottomRows(1) = m_quNodes;

        this->eval3D_into(m_points,m_evalPoints);
        for (index_t i=0; i!=this->targetDim(); ++i) // components
        {
            res = 0.0;
            for (index_t k = 0; k != m_numGauss; ++k) // compute integral
                res += m_quWeights.at(k) * math::pow(m_quNodes(0,k),moment) * m_evalPoints(i,k);
            result(i,j) = res;
        }
    }
    return result;
}

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval3D_Linear(const gsMatrix<T>& u) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: points in R3
    // Output: material matrix in R9 (cols stacked in rows; every row corresponds to a point)
    gsMatrix<T> result(9,u.cols());

    // NOTE 1: if the input \a u is considered to be in physical coordinates
    // then we first need to invert the points to parameter space
    // m_patches.patch(0).invertPoints(u, m_map.points, 1e-8) which is not exact (!),
    // otherwise we just use the input paramteric points
    m_map.points = u.topRows(2);

    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    // NOTE 2: in the case that parametric value is needed it suffices
    // to evaluate Youngs modulus and Poisson's ratio at
    // \a u instead of _tmp.values[0].
    m_YoungsModulus->eval_into(m_map.values[0], m_Emat);
    m_PoissonRatio->eval_into(m_map.values[0], m_Nmat);

    for( index_t i=0; i< u.cols(); ++i )
    {
        gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

        F0.leftCols(2) = m_map.jacobian(i);
        F0.col(2)      = m_map.normal(i).normalized();
        F0 = F0.inverse();
        F0 = F0 * F0.transpose(); //3x3

        // Evaluate material properties on the quadrature point
        m_E = m_Emat(0,i);
        m_nu = m_Nmat(0,i);
        m_lambda = m_E * m_nu / ( (1. + m_nu)*(1.-2.*m_nu)) ;
        m_mu     = m_E / (2.*(1. + m_nu)) ;

        m_Cconstant = 2*m_lambda*m_mu/(m_lambda+2*m_mu);

        C(0,0) = m_Cconstant*F0(0,0)*F0(0,0) + 1*m_mu*(2*F0(0,0)*F0(0,0));
        C(1,1) = m_Cconstant*F0(1,1)*F0(1,1) + 1*m_mu*(2*F0(1,1)*F0(1,1));
        C(2,2) = m_Cconstant*F0(0,1)*F0(0,1) + 1*m_mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
        C(1,0) =
        C(0,1) = m_Cconstant*F0(0,0)*F0(1,1) + 1*m_mu*(2*F0(0,1)*F0(0,1));
        C(2,0) =
        C(0,2) = m_Cconstant*F0(0,0)*F0(0,1) + 1*m_mu*(2*F0(0,0)*F0(0,1));
        C(2,1) = C(1,2) = m_Cconstant*F0(0,1)*F0(1,1) + 1*m_mu*(2*F0(0,1)*F0(1,1));

        //gsDebugVar(C);
    }

    return result;
}

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval_Composite(const gsMatrix<T>& u, int moment) const
{
    gsMatrix<T> result(this->targetDim(), u.cols());
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

    // Initialize and result
    // gsMatrix<T> tmp(3,3);
    // tmp.setZero();
    result.resize( 9, 1 );


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
        m_Dmat(0,0) = m_E1 / (1-m_nu12*m_nu21);
        m_Dmat(1,1) = m_E2 / (1-m_nu12*m_nu21);;
        m_Dmat(2,2) = m_G12;
        m_Dmat(0,1) = m_nu21*m_E1 / (1-m_nu12*m_nu21);
        m_Dmat(1,0) = m_nu12*m_E2 / (1-m_nu12*m_nu21);
        m_Dmat(2,0) = m_Dmat(0,2) = m_Dmat(2,1) = m_Dmat(1,2) = 0.0;

        // Make transformation matrix
        m_Transform(0,0) = m_Transform(1,1) = math::pow(math::cos(m_phi),2);
        m_Transform(0,1) = m_Transform(1,0) = math::pow(math::sin(m_phi),2);
        m_Transform(2,0) = m_Transform(0,2) = m_Transform(2,1) = m_Transform(1,2) = math::sin(m_phi) * math::cos(m_phi);
        m_Transform(2,0) *= -2.0;
        m_Transform(2,1) *= 2.0;
        m_Transform(1,2) *= -1.0;
        m_Transform(2,2) = math::pow(math::cos(m_phi),2) - math::pow(math::sin(m_phi),2);

        // Compute laminate stiffness matrix
        m_Dmat = m_Transform.transpose() * m_Dmat * m_Transform;

        m_z = math::abs(m_z_mid - (m_t/2.0 + m_t_temp) ); // distance from mid-plane of plate

        switch (moment)
        {
            case 0:
                result.reshape(3,3) += m_Dmat * m_t; // A
                break;
            case 1:
                result.reshape(3,3) += m_Dmat * m_t*m_z; // B
                break;
            case 2:
                result.reshape(3,3) += m_Dmat * ( m_t*math::pow(m_z,2) + math::pow(m_t,3)/12.0 ); // D
                break;
        }

        m_t_temp += m_t;
    }

    GISMO_ASSERT(m_t_tot==m_t_temp,"Total thickness after loop is wrong. t_temp = "<<m_t_temp<<" and sum(thickness) = "<<m_t_tot);
    // gsDebugVar(tmp);

    m_map.points = u.topRows(2);

    result = result.replicate(1,u.cols());
    // TRANSFORMATION OF THE MATRIX FROM LOCAL CARTESIAN TO CURVILINEAR
    // Compute covariant bases in deformed and undeformed configuration
    static_cast<const gsFunction<T>&>(m_patches->piece(m_pIndex)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    for (index_t k = 0; k != u.cols(); ++k)
    {
        m_covBasis.resize(3,3);
        m_covBasis.leftCols(2) = m_map.jacobian(k);
        m_covBasis.col(2)      = m_map.normal(k).normalized();
        m_covMetric = m_covBasis.transpose() * m_covBasis;

        m_conMetric = m_covMetric.inverse();

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

        // Covariant to local cartesian
        m_Transform(0,0) = (m_e1.dot(m_a1))*(m_a1.dot(m_e1));
        m_Transform(0,1) = (m_e1.dot(m_a2))*(m_a2.dot(m_e2));
        m_Transform(0,2) = 2*(m_e1.dot(m_a1))*(m_a2.dot(m_e1));
        // Row 1
        m_Transform(1,0) = (m_e2.dot(m_a1))*(m_a1.dot(m_e2));
        m_Transform(1,1) = (m_e2.dot(m_a2))*(m_a2.dot(m_e2));
        m_Transform(1,2) = 2*(m_e2.dot(m_a1))*(m_a2.dot(m_e2));
        // Row 2
        m_Transform(2,0) = (m_e1.dot(m_a1))*(m_a1.dot(m_e2));
        m_Transform(2,1) = (m_e1.dot(m_a2))*(m_a2.dot(m_e2));
        m_Transform(2,2) = (m_e1.dot(m_a1))*(m_a2.dot(m_e2)) + (m_e1.dot(m_a2))*(m_a1.dot(m_e2));

        m_Transform = m_Transform.inverse(); // !!!!
        result.reshapeCol(k,3,3) = m_Transform * result.reshapeCol(k,3,3);

        // Contravariant to local cartesian
        m_Transform(0,0) = (m_e1.dot(m_ac1))*(m_ac1.dot(m_e1));
        m_Transform(0,1) = (m_e1.dot(m_ac2))*(m_ac2.dot(m_e2));
        m_Transform(0,2) = 2*(m_e1.dot(m_ac1))*(m_ac2.dot(m_e1));
        // Row 1
        m_Transform(1,0) = (m_e2.dot(m_ac1))*(m_ac1.dot(m_e2));
        m_Transform(1,1) = (m_e2.dot(m_ac2))*(m_ac2.dot(m_e2));
        m_Transform(1,2) = 2*(m_e2.dot(m_ac1))*(m_ac2.dot(m_e2));
        // Row 2
        m_Transform(2,0) = (m_e1.dot(m_ac1))*(m_ac1.dot(m_e2));
        m_Transform(2,1) = (m_e1.dot(m_ac2))*(m_ac2.dot(m_e2));
        m_Transform(2,2) = (m_e1.dot(m_ac1))*(m_ac2.dot(m_e2)) + (m_e1.dot(m_ac2))*(m_ac1.dot(m_e2));

        result.reshapeCol(k,3,3) = result.reshapeCol(k,3,3) * m_Transform;
    }


    return result;
}

} // end namespace