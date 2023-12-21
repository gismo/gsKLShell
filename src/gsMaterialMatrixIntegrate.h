/** @file gsMaterialMatrixIntegrate.h

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
#include <gsKLShell/src/gsMaterialMatrixContainer.h>

namespace gismo
{

template <class T, enum MaterialOutput out> class gsMaterialMatrixIntegrateSingle;

template <class T, enum MaterialOutput out>
class gsMaterialMatrixIntegrate : public gsFunction<T>
{
public:
    /// Constructor
    gsMaterialMatrixIntegrate(  const gsMaterialMatrixContainer<T> & materialMatrices,
                            const gsFunctionSet<T> * deformed)
    :
    m_materialMatrices(materialMatrices),
    m_deformed(deformed)
    {
        for (index_t p = 0; p!=m_materialMatrices.size(); p++)
        {
            GISMO_ASSERT(materialMatrices.piece(p)!=nullptr,"Material matrix "<<p<<" is incomplete!");
            GISMO_ASSERT(materialMatrices.piece(p)!=NULL,"Material matrix "<<p<<" is incomplete!");
            GISMO_ASSERT(materialMatrices.piece(p)->initialized(),"Material matrix "<<p<<" is incomplete!");
        }

        this->_makePieces(deformed);
    }

    /// Constructor
    gsMaterialMatrixIntegrate(  gsMaterialMatrixBase<T> * materialMatrix,
                            const gsFunctionSet<T> * deformed)
    :
    m_materialMatrices(deformed->nPieces()),
    m_deformed(deformed)
    {
        GISMO_ASSERT(materialMatrix->initialized(),"Material matrix is incomplete!");
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
            m_materialMatrices.set(p,materialMatrix);
        this->_makePieces(deformed);
    }

    /// Constructor
    gsMaterialMatrixIntegrate(  const gsMaterialMatrixContainer<T> & materialMatrices,
                            const gsFunctionSet<T> * undeformed,
                            const gsFunctionSet<T> * deformed)
    :
    m_materialMatrices(materialMatrices),
    m_deformed(deformed)
    {
        this->_makePieces(undeformed,deformed);
    }

    /// Constructor
    gsMaterialMatrixIntegrate(  gsMaterialMatrixBase<T> * materialMatrix,
                            const gsFunctionSet<T> * undeformed,
                            const gsFunctionSet<T> * deformed)
    :
    m_materialMatrices(deformed->nPieces()),
    m_deformed(deformed)
    {
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
        {
            m_materialMatrices.set(p,materialMatrix);
            this->_makePieces(undeformed,deformed);
        }
    }

    /// Destructor
    ~gsMaterialMatrixIntegrate()
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
    void _makePieces(const gsFunctionSet<T> * deformed)
    {
        m_pieces.resize(deformed->nPieces());
        for (size_t p = 0; p!=m_pieces.size(); ++p)
            m_pieces.at(p) = new gsMaterialMatrixIntegrateSingle<T,out>(p,m_materialMatrices.piece(p),deformed);
    }

    void _makePieces(const gsFunctionSet<T> * undeformed, const gsFunctionSet<T> * deformed)
    {
        m_pieces.resize(m_deformed->nPieces());
        for (size_t p = 0; p!=m_pieces.size(); ++p)
            m_pieces.at(p) = new gsMaterialMatrixIntegrateSingle<T,out>(p,m_materialMatrices.piece(p),undeformed,deformed);
    }

protected:
    gsMaterialMatrixContainer<T> m_materialMatrices;
    const gsFunctionSet<T> * m_deformed;
    mutable std::vector<gsMaterialMatrixIntegrateSingle<T,out> *> m_pieces;
};

/**
 * @brief      This class serves as the integrator of material matrices, based on \ref gsMaterialMatrixBase
 *
 * @tparam     T     Real type
 * @tparam     out   Output type (see \ref MaterialOutput)
 *
 * @ingroup    KLShell
 */
template <class T, enum MaterialOutput out>
class gsMaterialMatrixIntegrateSingle : public gsFunction<T>
{
public:

    /// Constructor
    gsMaterialMatrixIntegrateSingle(  index_t patch,
                                gsMaterialMatrixBase<T> * materialMatrix,
                                const gsFunctionSet<T> * deformed);

    /// Constructor
    gsMaterialMatrixIntegrateSingle(  index_t patch,
                                gsMaterialMatrixBase<T> * materialMatrix,
                                const gsFunctionSet<T> * undeformed,
                                const gsFunctionSet<T> * deformed);

    gsMaterialMatrixIntegrateSingle(){};

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
                            _out==MaterialOutput::CauchyVectorN ||
                            _out==MaterialOutput::VectorM ||
                            _out==MaterialOutput::CauchyVectorM   , short_t>::type targetDim_impl() const { return 3; };

    /// Implementation of \ref targetDim for material tensors
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA ||
                            _out==MaterialOutput::MatrixB ||
                            _out==MaterialOutput::MatrixC ||
                            _out==MaterialOutput::MatrixD   , short_t>::type targetDim_impl() const { return 9; };

    /// Implementation of \ref targetDim for principal stress fields
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressN ||
                            _out==MaterialOutput::PStressM ||
                            _out==MaterialOutput::PStrainN ||
                            _out==MaterialOutput::PStrainM   , short_t>::type targetDim_impl() const { return 2; };

    /// Implementation of \ref targetDim for principal stretch fields
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Stretch   , short_t>::type targetDim_impl() const { return 3; };

    /// Implementation of \ref targetDim for principal stress directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::StretchDir, short_t>::type targetDim_impl() const { return 9; };

protected:
    /// Sets the patch index
    void setPatch(index_t p) {m_pIndex = p; }

public:
    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

private:
    /// Specialisation of \ref eval_into for densities
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Density   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the membrane stress tensor N
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN ||
                            _out==MaterialOutput::CauchyVectorN   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the flexural stress tensor M
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorM ||
                            _out==MaterialOutput::CauchyVectorM   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

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

    /// Specialisation of \ref eval_into for the membrane and flexural principle stresses
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStrainN ||
                            _out==MaterialOutput::PStrainM  , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the stretches
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Stretch   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Specialisation of \ref eval_into for the stretch directions
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::StretchDir, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

protected:

    /**
     * @brief      Gets the moment based on the output type
     *
     * for
     *      VectorN, the moment is 0
     *      VectorM, the moment is 1 (unless MatIntegrated==Constant or MatIntegrated==Integrated)
     *      MatrixA, the moment is 0
     *      MatrixB, the moment is 1
     *      MatrixC, the moment is 1
     *      MatrixD, the moment is 2
     *      PStressN,the moment is 0
     *      PStressM,the moment is 1
     *
     * @return     The moment.
     */
    T getMoment() const { return getMoment_impl<out>(); }

private:
    /// Implementation of \ref getMoment for MaterialOutput::VectorN; the moment is 0
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN ||
                            _out==MaterialOutput::PStressN||
                            _out==MaterialOutput::CauchyVectorN, T>::type getMoment_impl() const { return 0; };

    /// Implementation of \ref getMoment for MaterialOutput::VectorM; the moment is 1
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorM ||
                            _out==MaterialOutput::PStressM||
                            _out==MaterialOutput::CauchyVectorM, T>::type getMoment_impl() const { return 1; };

    /// Implementation of \ref getMoment for MaterialOutput::MatrixA; the moment is 0
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA, T>::type getMoment_impl() const { return 0; };

    /// Implementation of \ref getMoment for MaterialOutput::MatrixB; the moment is 1
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixB, T>::type getMoment_impl() const { return 1; };

    /// Implementation of \ref getMoment for MaterialOutput::MatrixC; the moment is 1
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixC, T>::type getMoment_impl() const { return 1; }; // must be 1

    /// Implementation of \ref getMoment for MaterialOutput::MatrixD; the moment is 2
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixD, T>::type getMoment_impl() const { return 2; }; // must be 2

    /// Implementation of \ref getMoment for MaterialOutput other than VectorN, VectorM, MatrixA, MatrixB, MatrixC, MatrixD, PStressN, PStressM
    template<enum MaterialOutput _out>
    typename std::enable_if<!(_out==MaterialOutput::VectorN         || _out==MaterialOutput::VectorM        ||
                              _out==MaterialOutput::CauchyVectorN   || _out==MaterialOutput::CauchyVectorM  ||
                              _out==MaterialOutput::MatrixA         || _out==MaterialOutput::MatrixB        ||
                              _out==MaterialOutput::MatrixC         || _out==MaterialOutput::MatrixD        ||
                              _out==MaterialOutput::PStressN        || _out==MaterialOutput::PStressM  )
                                                 , index_t>::type getMoment_impl() const { GISMO_NO_IMPLEMENTATION};

protected:

    /**
     * @brief      Integrates through-thickness using Gauss integration
     *
     * @param[in]  u       evaluation points (in plane)
     * @param[in]  moment  moment to be taken
     * @param      result  the result
     */
    void integrateZ_into(const gsMatrix<T> & u, const index_t moment, gsMatrix<T> & result) const;

    /**
     * @brief      Uses the top and bottom parts of the thickness to compute the integral exactly
     *
     * This function assumes that the function \f$ h(z,...) \f$ to be integrated is of the form \f$ h(z,...) = f(...) + z g(...)\f$ (no higher-order terms!). This implies that \f$f(...)\f$ is the even part and \f$z g(...)\f$ is the odd part of h. Writing out the thickness integral for a moment \f$\alpha\f$ gives:
      \f{eqnarray*}{
            \int_{-t/2}^{t/2} z^\alpha h(z,...)\:\text{d}\z &= \int_{-t/2}^{t/2} z^\alpha f(...) + z^{\alpha+1} g(...)\:\text{d}\z \\
             & = \frac{z^{\alpha+1}}{\alpha+1} f(...) + \frac{z^{\alpha+2}}{\alpha+2} g(...) \bigg\vert_{t/2}^{t/2} \\
             & = \frac{z^{\alpha+1}}{\alpha+1} f(...) + \frac{z^{\alpha+1}}{\alpha+2} zg(...) \bigg\vert_{t/2}^{t/2} \\
             & = \frac{1}{\alpha+1} \left[ \left(\frac{t}{2}\right)^{\alpha+1} - \left(-\frac{t}{2}\right)^{\alpha+1} \right]f(...) + \frac{1}{\alpha+2}\left[ \left(\frac{t}{2}\right)^{\alpha+2} - \left(-\frac{t}{2}\right)^{\alpha+2} \right] g(...) \bigg\vert_{t/2}^{t/2} \\
      \f}
      From this we observe that for \f$\alpha\f$ odd, the part of \f$f\f$ contributes, whereas for \f$\alpha\f$ even, the \f$g\f$ part contributes. Remember our assumption on the form of \f$h(z,...)\f$, then we can integrate this function by evaluating the following:
      \f{eqnarray*}{
            \int_{-t/2}^{t/2} z^\alpha h(z,...)\:\text{d}\z =
            \begin{dcases}
                \frac{1}{\alpha+1} \left[ \left(\frac{t}{2}\right)^{\alpha+1} - \left(-\frac{t}{2}\right)^{\alpha+1} \right]h(...) & \text{} if $\alpha$ is odd} \\
                \frac{1}{\alpha+2}\left[ \left(\frac{t}{2}\right)^{\alpha+1} - \left(-\frac{t}{2}\right)^{\alpha+1} \right] h(...) & \text{} if $\alpha$ is even}
            \end{dcases}
      \f}
     *
     * @param[in]  u       evaluation points (in plane)
     * @param[in]  moment  moment to be taken
     * @param      result  the result
     */
    void multiplyLinZ_into(const gsMatrix<T> & u, const index_t moment, gsMatrix<T> & result) const;

    /**
     * @brief      Multiplies the evaluation at z=0 with  2.0/(moment+1) * Thalf^(moment + 1)
     *
     * WARNING: Recommended use only for constant functions over thickness.
     *
     * The function assumes that moments are handled inside the evaluation. For example, when we want to integrate \f$f(z,...) = g(...) + zh(...)\f$ with different moments, then for moment 0, the function assumes that \f$g(...)\f$ is evaluated, for moment 1, the function assumes that nothing is returned and for moment 2, it assumes that \f$h\f$ is returned.
     *
     * The function is equivalent to \ref multiplyLinZ_into when moment = 2 here and moment = 1 in \ref multiplyLinZ_into (given that the correct function is returned).
     *
     *
     * @param[in]  u       evaluation points (in plane)
     * @param[in]  moment  moment to be taken
     * @param      result  the result
     */
    void multiplyZ_into(const gsMatrix<T> & u, index_t moment, gsMatrix<T> & result) const;

protected:

    /**
     * @brief      Integrateuates the base class in 3D
     *
     * @param[in]  u        The evaluation points (in plane)
     * @param[in]  Z        The through-thickness coordinate
     *
     * @return     Matrix ordered over Z and over u within
     */
    gsMatrix<T> _eval3D(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;


    /**
     * @brief      Integrateuates the base class in 3D, but on Z=0
     *
     * This function is primarily used in cases where the base class already integrates the material matrix or vectors
     *
     * @param[in]  u     The evaluation points (in plane)
     *
     * @return     Matrix with results
     */
    gsMatrix<T> _eval  (const gsMatrix<T>& u) const;

private:
    /// Specialisation of \ref eval3D for vectors
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN ||
                            _out==MaterialOutput::VectorM         , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;

    /// Specialisation of \ref eval3D for vectors
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::CauchyVectorN ||
                            _out==MaterialOutput::CauchyVectorM   , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;


    /// Specialisation of \ref eval3D for matrix
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA ||
                            _out==MaterialOutput::MatrixB ||
                            _out==MaterialOutput::MatrixC ||
                            _out==MaterialOutput::MatrixD   , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;

    /// Specialisation of \ref eval3D for vectors
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStress  ||
                            _out==MaterialOutput::PStressN ||
                            _out==MaterialOutput::PStressM   , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;

    /// Specialisation of \ref eval3D for other types
    template<enum MaterialOutput _out>
    typename std::enable_if<!(_out==MaterialOutput::VectorN         || _out==MaterialOutput::VectorM        ||
                              _out==MaterialOutput::CauchyVectorN   || _out==MaterialOutput::CauchyVectorM  ||
                              _out==MaterialOutput::MatrixA         || _out==MaterialOutput::MatrixB        ||
                              _out==MaterialOutput::MatrixC         || _out==MaterialOutput::MatrixD        ||
                              _out==MaterialOutput::PStressN        || _out==MaterialOutput::PStressM  )
                                                            , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
    { GISMO_NO_IMPLEMENTATION};

protected:
    index_t m_pIndex;
    gsMaterialMatrixBase<T> * m_materialMat;
    gsMatrix<T> m_z;
};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixIntegrate.hpp)
#endif

