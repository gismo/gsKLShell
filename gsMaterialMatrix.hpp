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
    m_map.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
}

// // Nonlinear material models
// template<class T>
// gsMaterialMatrix<T>::gsMaterialMatrix(  const gsFunctionSet<T> & mp,
//                                         const gsFunctionSet<T> & mp_def,
//                                         const gsFunction<T> & thickness,
//                                         const gsFunction<T> & YoungsModulus,
//                                         const gsFunction<T> & PoissonRatio
//                                         )
//                                         :
//                                         m_patches(&mp),
//                                         m_defpatches(&mp_def),
//                                         m_thickness(&thickness),
//                                         m_YoungsModulus(&YoungsModulus),
//                                         m_PoissonRatio(&PoissonRatio),
//                                         m_piece(nullptr)
// {
//     m_map.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
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
short_t gsMaterialMatrix<T>::targetDim() const { return 9; }

template <class T>
void gsMaterialMatrix<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = integrateZ(u);
}

template <class T>
void gsMaterialMatrix<T>::eval3D_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    switch (m_model)
    {
        case 0 :
            result = eval3D_Linear(u);
            break;
        // case 1 : eval3D_Incompressible(u,result);
        //     break;
        // case 2 : eval3D_Incompressible(u,result);
        //     break;
    }
}

template <class T>
gsMatrix<T> gsMaterialMatrix<T>::eval3D_Linear(const gsMatrix<T>& u) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: points in R3
    // Output: material matrix in R9 (cols stacked in rows)
    gsMatrix<T> result(9,1);

    // NOTE 1: if the input \a u is considered to be in physical coordinates
    // then we first need to invert the points to parameter space
    // m_patches.patch(0).invertPoints(u, m_map.points, 1e-8) which is not exact (!),
    // otherwise we just use the input paramteric points
    m_map.points = u.topRows(2);

    static_cast<const gsFunction<T>&>(m_patches->piece(0)).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    // NOTE 2: in the case that parametric value is needed it suffices
    // to evaluate Youngs modulus and Poisson's ratio at
    // \a u instead of _tmp.values[0].
    m_YoungsModulus->eval_into(m_map.values[0], m_Emat);
    m_PoissonRatio->eval_into(m_map.values[0], m_Nmat);

    result.resize( targetDim() , u.cols() );
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

// evaluates the material matrix over its last coordinate (thickness)
template <class T>
gsMatrix<T> gsMaterialMatrix<T>::integrateZ(const gsMatrix<T>& u, int moment) const
{
    // Input: points in R3
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

} // end namespace