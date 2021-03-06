/** @file gsMaterialMatrixEval.h

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

#include <gsCore/gsFunction.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

template <class T, enum MaterialOutput out>
class gsMaterialMatrixEval : public gsFunction<T>
{
public:

    // GISMO_CLONE_FUNCTION(gsMaterialMatrixBase)

    gsMaterialMatrixEval( gsMaterialMatrixBase<T> * materialMatrix);

    short_t domainDim() const;// { return 2; }

    // template OUT
    short_t targetDim() const;// { return getMoment_impl<out>(); }

private:
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Density   , short_t>::type targetDim_impl() const;// { return 1; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN ||
                            _out==MaterialOutput::VectorM   , short_t>::type targetDim_impl() const;// { return 3; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA ||
                            _out==MaterialOutput::MatrixB ||
                            _out==MaterialOutput::MatrixC ||
                            _out==MaterialOutput::MatrixD   , short_t>::type targetDim_impl() const;// { return 9; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressN ||
                            _out==MaterialOutput::PStressM  , short_t>::type targetDim_impl() const;// { return 2; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Stretch   , short_t>::type targetDim_impl() const;// { return 3; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::StretchDir, short_t>::type targetDim_impl() const;// { return 9; };

public:
    const gsFunction<T> & piece(const index_t p) const
    {
        m_piece = new gsMaterialMatrixEval(*this);
        m_piece->setPatch(p);
        return *m_piece;
    }

protected:
    void setPatch(index_t p) {m_pIndex = p; }

public:

    ~gsMaterialMatrixEval() { delete m_piece; }


    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

private:
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Density   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorM   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA ||
                            _out==MaterialOutput::MatrixB ||
                            _out==MaterialOutput::MatrixC ||
                            _out==MaterialOutput::MatrixD   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressN ||
                            _out==MaterialOutput::PStressM  , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::Stretch   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::StretchDir, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

public:

    T getMoment() const { return getMoment_impl<out>(); }

private:
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN, T>::type getMoment_impl() const { return 0; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorM, T>::type getMoment_impl() const
    {
        // if (m_materialMat->material()==Material::SvK)
        //     return 2;
        // else
            return 1;

    };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA, T>::type getMoment_impl() const { return 0; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixB, T>::type getMoment_impl() const { return 1; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixC, T>::type getMoment_impl() const { return 1; }; // must be 1
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixD, T>::type getMoment_impl() const { return 2; }; // must be 2
    template<enum MaterialOutput _out>
    typename std::enable_if<!(_out==MaterialOutput::VectorN || _out==MaterialOutput::VectorM ||
                              _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB ||
                              _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD ||
                              _out==MaterialOutput::PStressN|| _out==MaterialOutput::PStressM  )
                                                 , index_t>::type getMoment_impl() const { GISMO_NO_IMPLEMENTATION};
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressN, T>::type getMoment_impl() const { return 0; };
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressM, T>::type getMoment_impl() const { return 1; };

protected:
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

// private:


public:

    gsMatrix<T> eval3D(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;
    gsMatrix<T> eval  (const gsMatrix<T>& u) const;

private:
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::VectorN ||
                            _out==MaterialOutput::VectorM   , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::MatrixA ||
                            _out==MaterialOutput::MatrixB ||
                            _out==MaterialOutput::MatrixC ||
                            _out==MaterialOutput::MatrixD   , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;
    template<enum MaterialOutput _out>
    typename std::enable_if<_out==MaterialOutput::PStressN ||
                            _out==MaterialOutput::PStressM  , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const;

    template<enum MaterialOutput _out>
    typename std::enable_if<!(_out==MaterialOutput::VectorN || _out==MaterialOutput::VectorM ||
                              _out==MaterialOutput::MatrixA || _out==MaterialOutput::MatrixB ||
                              _out==MaterialOutput::MatrixC || _out==MaterialOutput::MatrixD ||
                              _out==MaterialOutput::PStressN|| _out==MaterialOutput::PStressM  )
                                                            , gsMatrix<T>>::type eval3D_impl(const gsMatrix<T>& u, const gsMatrix<T>& Z) const
    { GISMO_NO_IMPLEMENTATION};

    // // template<short_t num=0>
    //  virtual gsMaterialMatrixBase * makeMatrix(short_t num) = 0;
    // // template<short_t num=0>
    //  virtual gsMaterialMatrixBase * makeVector(short_t num) = 0;
    // // template<short_t num=0>
    //  virtual gsMaterialMatrixBase * makePrincipleStress(short_t num) = 0;

    //  virtual gsMaterialMatrixBase * makeDensity() = 0;
    //  virtual gsMaterialMatrixBase * makeStretch() = 0;
    //  virtual gsMaterialMatrixBase * makeDirections() = 0;
    //
protected:
    gsMaterialMatrixBase<T> * m_materialMat;
    mutable gsMaterialMatrixEval<T,out> * m_piece;
    index_t m_pIndex;


};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixEval.hpp)
#endif

