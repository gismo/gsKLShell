/** @file gsMaterialMatrixEval.hpp

    @brief Provides an evaluator for material matrices for thin shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsAssembler/gsGaussRule.h>

namespace gismo
{
template <class T, enum MaterialOutput out>
gsMaterialMatrixEvalSingle<T,out>::gsMaterialMatrixEvalSingle(    index_t patch,
                                                        gsMaterialMatrixBase<T> * materialMatrix,
                                                        const gsFunctionSet<T> * deformed,
                                                        const gsMatrix<T> z
                                                   )
:
m_pIndex(patch),
m_materialMat(materialMatrix),
m_piece(nullptr),
m_z(z)
{
    m_materialMat->setDeformed(deformed);

    if (m_z.cols()==0 || m_z.rows()==0)
    {
        m_z.resize(1,1);
        m_z.setZero();
    }
    GISMO_ASSERT(z.cols()==1,"Z coordinates should be provided row-wise in one column");

    // m_materialMat = new gsMaterialMatrix(materialMatrix);
}

template <class T, enum MaterialOutput out>
void gsMaterialMatrixEvalSingle<T,out>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
/// Non-parallel evaluation
#pragma omp critical (gsMaterialMatrixEvalSingle_eval_into)
    this->eval_into_impl<out>(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Density, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->density_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorN || _out==MaterialOutput::VectorM || _out==MaterialOutput::Generic, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = m_materialMat->eval3D_vector(m_pIndex,u,m_z.replicate(1,u.cols()),_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<   _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB
                        || _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = m_materialMat->eval3D_matrix(m_pIndex,u,m_z.replicate(1,u.cols()),_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStressN || _out==MaterialOutput::PStressM, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = m_materialMat->eval3D_pstress(m_pIndex,u,m_z.replicate(1,u.cols()),_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Stretch, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretch_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::StretchDir, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretchDir_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Transformation, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->transform_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Thickness, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->thickness_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Parameters, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->parameters_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Deformation, void>::type
gsMaterialMatrixEvalSingle<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->deformation_into(m_pIndex,u,result);
}


} // end namespace
