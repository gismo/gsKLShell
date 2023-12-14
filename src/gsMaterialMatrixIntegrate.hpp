/** @file gsMaterialMatrixIntegrate.hpp

    @brief Provides an evaluator for material matrices for thin shells

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

#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsAssembler/gsGaussRule.h>

namespace gismo
{

// Linear material models
template <class T, enum MaterialOutput out>
gsMaterialMatrixIntegrateSingle<T,out>::gsMaterialMatrixIntegrateSingle(index_t patch,
                                                                        gsMaterialMatrixBase<T> * materialMatrix, //??
                                                                        const gsFunctionSet<T> * deformed
                                                                        )
:
m_pIndex(patch),
m_materialMat(materialMatrix)
{
    m_materialMat->setDeformed(deformed);
    // m_materialMat = new gsMaterialMatrix(materialMatrix);
}

// Linear material models
template <class T, enum MaterialOutput out>
gsMaterialMatrixIntegrateSingle<T,out>::gsMaterialMatrixIntegrateSingle(index_t patch,
                                                                        gsMaterialMatrixBase<T> * materialMatrix, //??
                                                                        const gsFunctionSet<T> * undeformed,
                                                                        const gsFunctionSet<T> * deformed
                                                                        )
:
m_pIndex(patch),
m_materialMat(materialMatrix)
{
    m_materialMat->setUndeformed(undeformed);
    m_materialMat->setDeformed(deformed);
    // m_materialMat = new gsMaterialMatrix(materialMatrix);
}

template <class T, enum MaterialOutput out>
void gsMaterialMatrixIntegrateSingle<T,out>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
/// Non-parallel evaluation
// #pragma omp critical (gsMaterialMatrixIntegrateSingle_eval_into)
    this->eval_into_impl<out>(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Density, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->density_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorN ||
                        _out==MaterialOutput::CauchyVectorN, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    if (m_materialMat->isVecIntegrated() == MatIntegration::NotIntegrated)
        this->integrateZ_into(u,getMoment(),result);
    else if (m_materialMat->isVecIntegrated() == MatIntegration::Integrated)
        result = this->_eval(u);
    else if (m_materialMat->isVecIntegrated() == MatIntegration::Constant)
        this->multiplyZ_into(u,0,result);
    else if (m_materialMat->isVecIntegrated() == MatIntegration::Linear)
        this->multiplyLinZ_into(u,getMoment(),result);
    else
        GISMO_ERROR("Integration status unknown");
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorM ||
                        _out==MaterialOutput::CauchyVectorM, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    if (m_materialMat->isVecIntegrated() == MatIntegration::NotIntegrated)
        this->integrateZ_into(u,getMoment(),result);
    else if (m_materialMat->isVecIntegrated() == MatIntegration::Integrated)
        result = this->_eval(u);
    else if (m_materialMat->isVecIntegrated() == MatIntegration::Constant)
        this->multiplyZ_into(u,2,result);
    else if (m_materialMat->isVecIntegrated() == MatIntegration::Linear)
        this->multiplyLinZ_into(u,getMoment(),result);
    else
        GISMO_ERROR("Integration status unknown");
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<   _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB
                        || _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    if (m_materialMat->isMatIntegrated() == MatIntegration::NotIntegrated)
        this->integrateZ_into(u,getMoment(),result);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Integrated)
        result = this->_eval(u);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Constant)
        this->multiplyZ_into(u,getMoment(),result);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Linear)
        this->multiplyLinZ_into(u,getMoment(),result);
    else
        GISMO_ERROR("Integration status unknown");
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStressN || _out==MaterialOutput::PStressM, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    if (m_materialMat->isMatIntegrated() == MatIntegration::NotIntegrated)
        this->integrateZ_into(u,getMoment(),result);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Integrated)
        result = this->_eval(u);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Constant)
        this->multiplyZ_into(u,getMoment(),result);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Linear)
        this->multiplyLinZ_into(u,getMoment(),result);
    else
        GISMO_ERROR("Integration status unknown");
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStrainN || _out==MaterialOutput::PStrainM, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    if (m_materialMat->isMatIntegrated() == MatIntegration::NotIntegrated)
        this->integrateZ_into(u,getMoment(),result);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Integrated)
        result = this->_eval(u);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Constant)
        this->multiplyZ_into(u,getMoment(),result);
    else if (m_materialMat->isMatIntegrated() == MatIntegration::Linear)
        this->multiplyLinZ_into(u,getMoment(),result);
    else
        GISMO_ERROR("Integration status unknown");
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Stretch, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretch_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::StretchDir, void>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretchDir_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
void gsMaterialMatrixIntegrateSingle<T,out>::integrateZ_into(const gsMatrix<T>& u, const index_t moment, gsMatrix<T> & result) const
{
    // Input: points in R2
    // Ouput: results in targetDim

    // Perform integration
    index_t numGauss = m_materialMat->options().askInt("NumGauss",4);

    result.resize(this->targetDim(),u.cols());
    result.setZero();

    gsGaussRule<T> gauss(numGauss);
    gsMatrix<T> z(numGauss,u.cols());
    gsMatrix<T> w(numGauss,u.cols());
    gsMatrix<T> vals;
    // m_points3D.resize(1,numGauss);

    gsMatrix<T> Tmat;
    m_materialMat->thickness_into(m_pIndex,u,Tmat);
    T Thalf;
    // pre-compute Z
    for (index_t k = 0; k != u.cols(); ++k) // for all points
    {
        gsMatrix<T> quNodes(1,numGauss);
        gsVector<T> quWeights(numGauss);
        // set new integration point
        // Thalf = Tmat(0,k)/2.0;
        Thalf = 0.5;
        gauss.mapTo(-Thalf,Thalf,quNodes,quWeights);
        w.col(k)=quWeights;
        z.col(k)=quNodes.transpose();
    }
    vals = this->_eval3D(u,z);
    T res;
    for (index_t k = 0; k != u.cols(); ++k) // for all points
    {
        for (index_t i=0; i!=this->targetDim(); ++i) // components
        {
            res = 0.0;
            for (index_t j = 0; j != numGauss; ++j) // compute integral
                res += w(j, k) * math::pow(z(j, k) * Tmat(0, k), moment) * vals(i, j * u.cols() + k) * Tmat(0, k); // int_z=[-1/2,1/2] (xi*t)^moment * <quantity(xi)> * t
            result(i,k) = res;
        }

    }
}


/*
 * WARNING: This function assumes the function to be integrated: f(z,...) = g(,...) + z h(...)
    Note: z is assumed to be [-1,1]
 */
template <class T, enum MaterialOutput out>
void gsMaterialMatrixIntegrateSingle<T,out>::multiplyLinZ_into(const gsMatrix<T>& u, const index_t moment, gsMatrix<T>& result) const
{
    // Input: points in R2
    // Ouput: results in targetDim
    result.resize(this->targetDim(),u.cols());

    gsMatrix<T> z(2,u.cols());
    gsMatrix<T> Tmat;
    m_materialMat->thickness_into(m_pIndex,u,Tmat);

    T Thalf;
    // pre-compute Z
    for (index_t k = 0; k != u.cols(); ++k) // for all points
    {
        Thalf = 1.0 / 2.0;
        // Thalf = Tmat(0, k) / 2.0;
        z.col(k)<<Thalf,-Thalf;
    }

    T fac = (moment % 2 == 0) ? 0. : 1.;

    gsMatrix<T> vals = this->_eval3D(u,z);
    for (index_t k = 0; k != u.cols(); ++k) // for all points
    {
            // 1/(alpha+1) * [ (t/2)^(alpha+1) * g(...)  - (-t/2)^(alpha+1) * g(...) ]
            // 1/(alpha+2) * [ (t/2)^(alpha+1) * h(...)  - (-t/2)^(alpha+1) * h(...) ]
            result.col(k) = 1.0/(moment+fac+1) *
                            ( math::pow( z(0,k) * Tmat(0, k), moment + 1) * vals.col(0*u.cols() + k)
                                - math::pow( z(1,k) * Tmat(0, k), moment + 1) * vals.col(1*u.cols() + k) );
    }

}

template <class T, enum MaterialOutput out>
void gsMaterialMatrixIntegrateSingle<T,out>::multiplyZ_into(const gsMatrix<T>& u, index_t moment, gsMatrix<T>& result) const
{
    // Input: points in R2
    // Ouput: results in targetDim
    result.resize(this->targetDim(),u.cols());

    if (moment % 2 != 0)  //then the moment is odd
        result.setZero();
    else                    // then the moment is even
    {
        T fac;
        gsMatrix<T> Tmat;
        m_materialMat->thickness_into(m_pIndex,u,Tmat);
        T Thalf;
        gsMatrix<T> vals = this->_eval(u);
        for (index_t k = 0; k != u.cols(); ++k) // for all points
        {
            Thalf = Tmat(0, k) / 2.0;
            fac = 2.0/(moment+1) * math::pow( Thalf , moment + 1);
            result.col(k) = vals.col(k) * fac;
        }

    }
}

template <class T, enum MaterialOutput out>
gsMatrix<T> gsMaterialMatrixIntegrateSingle<T,out>::_eval(const gsMatrix<T>& u) const
{
    gsMatrix<T> Z(1,u.cols());
    Z.setZero();
    return this->_eval3D(u,Z);
}

template <class T, enum MaterialOutput out>
gsMatrix<T> gsMaterialMatrixIntegrateSingle<T,out>::_eval3D(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return this->eval3D_impl<out>(u,Z);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<   _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB
                        || _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD, gsMatrix<T>>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return m_materialMat->eval3D_matrix(m_pIndex,u,Z,_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorN ||
                        _out==MaterialOutput::VectorM   , gsMatrix<T>>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return m_materialMat->eval3D_vector(m_pIndex,u,Z,_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::CauchyVectorN ||
                        _out==MaterialOutput::CauchyVectorM, gsMatrix<T>>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return m_materialMat->eval3D_CauchyVector(m_pIndex,u,Z,_out);
}


template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStress  ||
                        _out==MaterialOutput::PStressN ||
                        _out==MaterialOutput::PStressM, gsMatrix<T>>::type
gsMaterialMatrixIntegrateSingle<T,out>::eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return m_materialMat->eval3D_pstress(m_pIndex,u,Z,_out);
}

} // end namespace
