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

/*
    To Do [updated 16-06-2020]:
    - Make beta (compressible materials) and material parameters universal for all integration points over the thickness. So get them out of the dPsi functions etc and move them into the integration loops as global variables.

*/



#pragma once

#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsAssembler/gsGaussRule.h>

namespace gismo
{

// Linear material models
template <class T, enum MaterialOutput out>
gsMaterialMatrixEval<T,out>::gsMaterialMatrixEval(  gsMaterialMatrixBase<T> * materialMatrix, //??
                                                    const gsFunctionSet<T> & deformed,
                                                    const gsMatrix<T> z
                                                   )
:
m_materialMat(materialMatrix),
m_z(z),
m_piece(nullptr)
{
    m_pIndex = 0;
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
void gsMaterialMatrixEval<T,out>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    this->eval_into_impl<out>(u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Density, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->density_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::VectorN || _out==MaterialOutput::VectorM || _out==MaterialOutput::Generic, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = m_materialMat->eval3D_vector(m_pIndex,u,m_z.replicate(1,u.cols()),_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<   _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB
                        || _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = m_materialMat->eval3D_matrix(m_pIndex,u,m_z.replicate(1,u.cols()),_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::PStressN || _out==MaterialOutput::PStressM, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = m_materialMat->eval3D_pstress(m_pIndex,u,m_z.replicate(1,u.cols()),_out);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Stretch, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretch_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::StretchDir, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->stretchDir_into(m_pIndex,u,result);
}

template <class T, enum MaterialOutput out>
template <enum MaterialOutput _out>
typename std::enable_if<_out==MaterialOutput::Transformation, void>::type
gsMaterialMatrixEval<T,out>::eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_materialMat->transform_into(m_pIndex,u,result);
}


} // end namespace
