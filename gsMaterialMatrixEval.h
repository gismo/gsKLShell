/** @file gsMaterialMatrixEval.h

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

#include <gsCore/gsFunction.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{


/**
 * @brief      This class serves as the evaluator of material matrices, based on \ref gsMaterialMatrixBase
 *
 * @tparam     T     Real tyoe
 * @tparam     out   Output type (see \ref MaterialOutput)
 *
 * @ingroup    MaterialMatrix
 */
template <class T, enum MaterialOutput out>
class gsMaterialMatrixEval : public gsFunction<T>
{
public:

    /// Constructor
    gsMaterialMatrixEval(  gsMaterialMatrixBase<T> * materialMatrix,
                                const gsFunctionSet<T> & deformed,
                                const gsMatrix<T> z);

    /// Domain dimension, always 2 for shells
    short_t domainDim() const {return 2;}

    /**
     * @brief      Target dimension
     *
     * For a scalar (e.g. density) the target dimension is 1, for a vector (e.g. stress tensor in Voight notation) the target dimension is 3 and for a matrix (e.g. the material matrix) the target dimension is 9, which can be reshaped to a 3x3 matrix.
     *
     * @return     Returns the target dimension depending on the specified type (scalar, vector, matrix etc.)
     */
    short_t targetDim() const { return targetDim_impl<out>(); }


private:
    /// Implementation of \ref targetDim for densities
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Density   , short_t>::type targetDim_impl() const { return 1; };

    /// Implementation of \ref targetDim for stress tensors
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN ||
                            _out==MaterialOutput::VectorM ||
                            _out==MaterialOutput::Generic  , short_t>::type targetDim_impl() const { return 3; };

    /// Implementation of \ref targetDim for material tensors
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA ||
                            _out==MaterialOutput::MatrixB ||
                            _out==MaterialOutput::MatrixC ||
                            _out==MaterialOutput::MatrixD   , short_t>::type targetDim_impl() const { return 9; };

    /// Implementation of \ref targetDim for principal stress fields
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressN ||
                            _out==MaterialOutput::PStressM  , short_t>::type targetDim_impl() const { return 3; };

    /// Implementation of \ref targetDim for principal stretch fields
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Stretch   , short_t>::type targetDim_impl() const { return 3; };

    /// Implementation of \ref targetDim for principal stress directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::StretchDir, short_t>::type targetDim_impl() const { return 9; };

    /// Implementation of \ref targetDim for principal stress directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Transformation, short_t>::type targetDim_impl() const { return 9; };


public:
    /// Implementation of piece, see \ref gsFunction
    const gsFunction<T> & piece(const index_t p) const
    {
        m_piece = new gsMaterialMatrixEval(*this);
        m_piece->setPatch(p);
        return *m_piece;
    }

protected:
    /// Sets the patch index
    void setPatch(index_t p) {m_pIndex = p; }

public:
    /// Destructor
    ~gsMaterialMatrixEval() { delete m_piece; }

    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

private:
    /// Specialisation of \ref eval_into for densities
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Density   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the membrane stress tensor N, M and generic stress
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN ||
                            _out==MaterialOutput::VectorM ||
                            _out==MaterialOutput::Generic , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the moments of the material matrices
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA ||
                            _out==MaterialOutput::MatrixB ||
                            _out==MaterialOutput::MatrixC ||
                            _out==MaterialOutput::MatrixD   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the membrane and flexural principle stresses
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressN ||
                            _out==MaterialOutput::PStressM  , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the stretches
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Stretch   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the stretch directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::StretchDir, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the basis transformation
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Transformation, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

protected:
    gsMaterialMatrixBase<T> * m_materialMat;
    gsMatrix<T> m_z;
    mutable gsMaterialMatrixEval<T,out> * m_piece;
    index_t m_pIndex;



};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixEval.hpp)
#endif

