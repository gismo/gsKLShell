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

#include <gsKLShell/gsMaterialMatrixEval.h>

namespace gismo
{

// Linear material models
template <class T, enum MaterialOutput out>
gsMaterialMatrixEval<T,out>::gsMaterialMatrixEval( gsMaterialMatrixBase<T> * materialMatrix //??
                                        )
:
m_materialMat(materialMatrix),
m_piece(nullptr)
{
    // m_materialMat = new gsMaterialMatrix(materialMatrix);
}

template <class T, enum MaterialOutput out>
short_t gsMaterialMatrixEval<T,out>::domainDim() const { return 2; }

template <class T, enum MaterialOutput out>
short_t gsMaterialMatrixEval<T,out>::targetDim() const
{
    return targetDim_impl<out>();
}


template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Density, short_t>::type
gsMaterialMatrixEval<T,out>::targetDim_impl() const
{
    return 1;
}


template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorN || _out==MaterialOutput::VectorM, short_t>::type
gsMaterialMatrixEval<T,out>::targetDim_impl() const
{
    return 3;
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<   _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB
                        || _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD, short_t>::type
gsMaterialMatrixEval<T,out>::targetDim_impl() const
{
    return 9;
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStressN || _out==MaterialOutput::PStressM, short_t>::type
gsMaterialMatrixEval<T,out>::targetDim_impl() const
{
    return 2;
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Stretch, short_t>::type
gsMaterialMatrixEval<T,out>::targetDim_impl() const
{
    return 3;
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::StretchDir, short_t>::type
gsMaterialMatrixEval<T,out>::targetDim_impl() const
{
    return 9;
}



    // if (m_outputType==2)
    //     return 9;
    // else if (m_outputType==1)
    //     return 3;
    // else if (m_outputType==0)
    //     return 1;
    // else if (m_outputType==9)
    //     return 3;
    // else if (m_outputType==10)
    //     return 2;
    // else if (m_outputType==11)
    //     return 9;
    // else
    // {
    //     GISMO_ERROR("This option is unknown");
    //     return 1;
    // }

template <class T, enum MaterialOutput out>
void gsMaterialMatrixEval<T,out>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->eval_into_impl<out>(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Density, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->setOutputType(0);
    m_materialMat->density_into(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorN || _out==MaterialOutput::VectorM, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->integrateZ_into(u,getMoment(),result);
    // m_materialMat->setOutputType(1);
    // m_materialMat->setOutput(getMoment());
    // m_materialMat->eval_into_NP(u,result);
    // m_materialMat->eval_into_vector(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<   _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB
                        || _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->integrateZ_into(u,getMoment(),result);
    // m_materialMat->setOutputType(2);
    // m_materialMat->setOutput(getMoment());
    // m_materialMat->eval_into_NP(u,result);
    // m_materialMat->eval_into_matrix(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStressN || _out==MaterialOutput::PStressM, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->integrateZ_into(u,getMoment(),result);
    // m_materialMat->setOutputType(10);
    // m_materialMat->setOutput(getMoment());
    // m_materialMat->eval_into_NP(u,result);
    // m_materialMat->eval_into_pstress(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Stretch, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretch_into(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::StretchDir, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretchDir_into(u,result);
}

template <class T, enum MaterialOutput out>
void gsMaterialMatrixEval<T,out>::integrateZ_into(const gsMatrix<T>& u, const index_t moment, gsMatrix<T> & result) const
{
    // Input: points in R2
    // Ouput: results in targetDim

/*
    !!!!!!!!!!!!!!!!!!!
    if ((m_output==1) && (m_outputType==1))
        m_moment = 2; // NEEDED SINCE m_moment=2 IS FOR THE OUTPUT OF THE M TENSOR, WHICH IN FACT HAS MOMENT 2. THIS IS BY CHOICE OF THE COMPUTATION OF THE STRAINS IN THE Sij() FUNCTION
*/
    // Perform integration
    index_t numGauss = m_materialMat->options().askInt("NumGauss",4);

    result.resize(this->targetDim(),u.cols());
    result.setZero();

    gsGaussRule<T> gauss = gsGaussRule<T>(numGauss);
    gsMatrix<T> z(numGauss,u.cols());
    gsMatrix<T> w(numGauss,u.cols());
    gsMatrix<T> vals;
    // m_points3D.resize(1,numGauss);

    gsMatrix<T> Tmat;
    m_materialMat->thickness_into(u,Tmat);
    T Thalf;
    // pre-compute Z
    for (index_t k = 0; k != u.cols(); ++k) // for all points
    {
        gsMatrix<T> quNodes(1,numGauss);
        gsVector<T> quWeights(numGauss);
        // set new integration point
        Thalf = Tmat(0,k)/2.0;
        gauss.mapTo(-Thalf,Thalf,quNodes,quWeights);
        w.col(k)=quWeights;
        z.col(k)=quNodes.transpose();
    }

    vals = this->eval3D(u,z);

    T res;
    for (index_t k = 0; k != u.cols(); ++k) // for all points
    {
        for (index_t i=0; i!=this->targetDim(); ++i) // components
        {
            res = 0.0;
            for (index_t j = 0; j != numGauss; ++j) // compute integral
                res += w(j,k) * math::pow(z(j,k),moment) * vals(i,j*u.cols() + k);
            result(i,k) = res;
        }

    }
}

template <class T, enum MaterialOutput out>
gsMatrix<T> gsMaterialMatrixEval<T,out>::eval3D(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return this->eval3D_impl<out>(u,Z);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<   _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB
                        || _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD, gsMatrix<T>>::type
gsMaterialMatrixEval<T,out>::eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return m_materialMat->eval3D_matrix(u,Z);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorN || _out==MaterialOutput::VectorM, gsMatrix<T>>::type
gsMaterialMatrixEval<T,out>::eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return m_materialMat->eval3D_vector(u,Z);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStressN || _out==MaterialOutput::PStressM, gsMatrix<T>>::type
gsMaterialMatrixEval<T,out>::eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
{
    return m_materialMat->eval3D_pstress(u,Z);
}


} // end namespace
