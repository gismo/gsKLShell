/** @file gsMaterialMatrixTFT.h

    @brief Provides linear material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{


/**
 * @brief      This class defines a linear material
 *
 * @tparam     dim   The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T     Real type
 *
 * @ingroup    MaterialMatrix
 *
 */
template <  short_t dim,
            class T,
            bool linear=false
         >
class gsMaterialMatrixTFT : public gsMaterialMatrixBaseDim<dim,T>
{
public:
    using Base = gsMaterialMatrixBaseDim<dim,T>;

    /**
     * @brief      Constructor without deformed multipatch and density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrixTFT(gsMaterialMatrixBaseDim<dim,T> * materialmatrix)
    :
    m_materialMat(materialmatrix)
    {
        // if (!dynamic_cast<const gsMaterialMatrixLinear<dim,T> *>(m_materialMat))
            // GISMO_ERROR("Material matrix must be linear");
    }

    /// Destructor
    gsMaterialMatrixTFT() { }

    /// See \ref gsMaterialMatrixBase for details
    /// Here, we use the integrated matrix and vector to make our computations!
    inline enum MatIntegration isMatIntegrated() const override
    // {return m_materialMat->isMatIntegrated(); gsWarn<<"Change this to a constant one, like for linear. So that no things out of the mid-plane are needed"; }
    {return MatIntegration::Linear;} // Multiplies with z t or its moment

    /// See \ref gsMaterialMatrixBase for details
    /// Here, we use the integrated matrix and vector to make our computations!
    inline enum MatIntegration isVecIntegrated() const override
    // {return m_materialMat->isVecIntegrated(); gsWarn<<"Change this to a constant one, like for linear. So that no things out of the mid-plane are needed"; }
    {return MatIntegration::Linear;} // Multiplies with z t or its moment

    // bool isLinear() { return isLinear_impl(); }
    // typename util::enable_if<U::Linear,bool>::type
    // eval_impl(const U & u, const index_t k)
    // { return m_materialMat; }


    /// See \ref gsMaterialMatrixBase for details
    gsOptionList & options()
    {return m_materialMat->options(); }

    /// See \ref gsMaterialMatrixBase for details
    void setOptions(gsOptionList opt)
    {m_materialMat->setOptions(opt); }

    /// See \ref gsMaterialMatrixBase for details
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { m_materialMat->density_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { m_materialMat->stretch_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { m_materialMat->stretchDir_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    { m_materialMat->thickness_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void covtransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    { m_materialMat->covtransform_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void contransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    { m_materialMat->contransform_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const
    { return m_materialMat->eval3D_pstress( patch,u,z,out); }

    /// See \ref gsMaterialMatrixBase for details
    void setParameters(const std::vector<gsFunction<T>*> &pars)
    { m_materialMat->setParameters(pars); }

    /// See \ref gsMaterialMatrixBase for details
    void info() const
    { m_materialMat->info(); }

    /// See \ref gsMaterialMatrixBase for details
    void setDeformed(const gsFunctionSet<T> * deformed)
    {
        gsDebugVar(deformed);
        m_materialMat->setDeformed(deformed);
    }

protected:
    template <bool _linear>
    typename std::enable_if< _linear, gsMatrix<T> >::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    template <bool _linear>
    typename std::enable_if<!_linear, gsMatrix<T> >::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    template <bool _linear>
    typename std::enable_if< _linear, gsMatrix<T> >::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    template <bool _linear>
    typename std::enable_if<!_linear, gsMatrix<T> >::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    static gsMatrix<T> _compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & dC);

    static gsMatrix<T> _compute_S(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S);


    // /// Computes theta
    // gsMatrix<T> eval_theta(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// Computes theta
    static gsMatrix<T> eval_theta(const gsMatrix<T> & Cs, const gsMatrix<T> & Ns);


public:
    /// Shared pointer for gsMaterialMatrixTFT
    typedef memory::shared_ptr< gsMaterialMatrixTFT > Ptr;

    /// Unique pointer for gsMaterialMatrixTFT
    typedef memory::unique_ptr< gsMaterialMatrixTFT > uPtr;

protected:
    mutable gsMaterialMatrixBaseDim<dim,T> * m_materialMat;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixTFT.hpp)
#endif
