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
#include <gsKLShell/gsMaterialMatrixContainer.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

template <class T, enum MaterialOutput out> class gsMaterialMatrixEvalSingle;

template <class T, enum MaterialOutput out>
class gsMaterialMatrixEval : public gsFunction<T>
{
public:
    /// Constructor
    gsMaterialMatrixEval(   const gsMaterialMatrixContainer<T> & materialMatrices,
                            const gsFunctionSet<T> * deformed,
                            const gsMatrix<T> z)
    :
    m_materialMatrices(materialMatrices),
    m_deformed(deformed),
    m_z(z)
    {
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
            m_pieces.push_back(new gsMaterialMatrixEvalSingle<T,out>(p,m_materialMatrices.piece(p),m_deformed,m_z));
    }

    /// Constructor
    gsMaterialMatrixEval(   gsMaterialMatrixBase<T> * materialMatrix,
                            const gsFunctionSet<T> * deformed,
                            const gsMatrix<T> z)
    :
    m_materialMatrices(deformed->nPieces()),
    m_deformed(deformed),
    m_z(z)
    {
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
        {
            m_materialMatrices.add(materialMatrix);
            m_pieces.push_back(new gsMaterialMatrixEvalSingle<T,out>(p,m_materialMatrices.piece(p),m_deformed,m_z));
        }
    }

    /// Destructor
    ~gsMaterialMatrixEval()
    {
        freeAll(m_pieces);
    }

    /// Domain dimension, always 2 for shells
    short_t domainDim() const {return 2;}

    /**
     * @brief      Target dimension
     *
     * For a scalar (e.g. density) the target dimension is 1, for a vector (e.g. stress tensor in Voight notation) the target dimension is 3 and for a matrix (e.g. the material matrix) the target dimension is 9, which can be reshaped to a 3x3 matrix.
     *
     * @return     Returns the target dimension depending on the specified type (scalar, vector, matrix etc.)
     */
    short_t targetDim() const { return this->piece(0).targetDim(); }

    /// Implementation of piece, see \ref gsFunction
    const gsFunction<T> & piece(const index_t p) const
    {
        return *m_pieces[p];
    }

    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

protected:
    gsMaterialMatrixContainer<T> m_materialMatrices;
    const gsFunctionSet<T> * m_deformed;
    gsMatrix<T> m_z;
    mutable std::vector<gsMaterialMatrixEvalSingle<T,out> *> m_pieces;
};

/**
 * @brief      This class serves as the evaluator of material matrices, based on \ref gsMaterialMatrixBase
 *
 * @tparam     T     Real type
 * @tparam     out   Output type (see \ref MaterialOutput)
 *
 * @ingroup    KLShell
 */
template <class T, enum MaterialOutput out>
class gsMaterialMatrixEvalSingle : public gsFunction<T>
{
public:

    /// Constructor
    gsMaterialMatrixEvalSingle( index_t patch,
                                gsMaterialMatrixBase<T> * materialMatrix,
                                const gsFunctionSet<T> * deformed,
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
                            _out==MaterialOutput::PStressM  , short_t>::type targetDim_impl() const { return 2; };

    /// Implementation of \ref targetDim for principal stretch fields
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Stretch   , short_t>::type targetDim_impl() const { return 3; };

    /// Implementation of \ref targetDim for principal stress directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::StretchDir, short_t>::type targetDim_impl() const { return 9; };

    /// Implementation of \ref targetDim for principal stress directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Transformation, short_t>::type targetDim_impl() const { return 9; };

    /// Implementation of \ref targetDim for the thickness
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Thickness, short_t>::type targetDim_impl() const { return 1; };

    /// Implementation of \ref targetDim for the parameters
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Parameters, short_t>::type targetDim_impl() const { return m_materialMat->numParameters(); };

    /// Implementation of \ref targetDim for principal stress directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Deformation, short_t>::type targetDim_impl() const { return 9; };

protected:
    /// Sets the patch index
    void setPatch(index_t p) {m_pIndex = p; }

public:
    /// Destructor
    ~gsMaterialMatrixEvalSingle() { }

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

    /// Specialisation of \ref eval_into for the thickness
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Thickness, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the parameters
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Parameters, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the deformation gradient
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Deformation, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;


protected:
    index_t m_pIndex;
    gsMaterialMatrixBase<T> * m_materialMat;
    gsMatrix<T> m_z;
};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixEval.hpp)
#endif

