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
            class T
         >
class gsMaterialMatrixTFT : public gsMaterialMatrixBaseDim<dim,T,true>
{
public:
    using Base = gsMaterialMatrixBaseDim<dim,T,true>;

    /**
     * @brief      Constructor without deformed multipatch and density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrixTFT(gsMaterialMatrixBaseDim<dim,T,true> * materialmatrix)
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
    inline enum MatIntegration isMatIntegrated() const {return m_materialMat->isMatIntegrated(); }

    /// See \ref gsMaterialMatrixBase for details
    /// Here, we use the integrated matrix and vector to make our computations!
    inline enum MatIntegration isVecIntegrated() const {return m_materialMat->isVecIntegrated(); }

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


public:
    /// Shared pointer for gsMaterialMatrixTFT
    typedef memory::shared_ptr< gsMaterialMatrixTFT > Ptr;

    /// Unique pointer for gsMaterialMatrixTFT
    typedef memory::unique_ptr< gsMaterialMatrixTFT > uPtr;

protected:
    mutable gsMaterialMatrixBaseDim<dim,T,true> * m_materialMat;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixTFT.hpp)
#endif
