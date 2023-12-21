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

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>

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
    /// Shared pointer for gsMaterialMatrixTFT
    typedef memory::shared_ptr< gsMaterialMatrixTFT > Ptr;

    /// Unique pointer for gsMaterialMatrixTFT
    typedef memory::unique_ptr< gsMaterialMatrixTFT > uPtr;

    /// Shared pointer to a \ref gsMaterialMatrixBase
    typedef typename gsMaterialMatrixBase<T>::Ptr material_ptr;

    /// Base class
    using Base = gsMaterialMatrixBaseDim<dim,T>;

    /// Shared pointer to a \ref gsFunctionSet
    typedef typename gsFunctionSet<T>::Ptr function_ptr;

    // Define clone functions
    GISMO_CLONE_FUNCTION(gsMaterialMatrixTFT)

public:

    /**
     * @brief      Constructs a TFT material matrix using a reference to another material matrix
     *
     * @param      materialMatrix  The original material matrix
     */
    gsMaterialMatrixTFT(gsMaterialMatrixBase<T> * materialMatrix)
    :
    gsMaterialMatrixTFT(memory::make_shared_not_owned(materialMatrix))
    {}

    /**
     * @brief      Constructs a TFT material matrix using a reference to another material matrix
     *
     * @param      materialMatrix  The original material matrix
     */
    gsMaterialMatrixTFT(gsMaterialMatrixBase<T> & materialMatrix)
    :
    gsMaterialMatrixTFT(memory::make_shared(materialMatrix.clone().release()))
    {}

    /**
     * @brief      Constructs a TFT material matrix using a reference to another material matrix
     *
     * @param      materialMatrix  The original material matrix
     */
    gsMaterialMatrixTFT(const material_ptr & materialMatrix)
    :
    m_materialMat(give(materialMatrix))
    {
        if (materialMatrix.get()->hasUndeformed())
            this->setUndeformed(materialMatrix.get()->getUndeformed());
        if (materialMatrix.get()->hasDeformed())
            this->setDeformed(materialMatrix.get()->getDeformed());

        this->setThickness(materialMatrix.get()->getThickness());

        m_options.addReal("SlackMultiplier","Multiplies the original value of the matrix for the slack state",0);
        m_options.addSwitch("Explicit","Explicit iterations; use tension field from a fixed deformed geometry that does not change when calling setDeformed",false);
    }

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
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { m_materialMat->density_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    { m_materialMat->thickness_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// Returns the principal stresses stored in the underlying material model where the tension-field is non-slack
    ///
    /// @param[in]  patch  The patch
    /// @param[in]  u      Evaluation points
    /// @param[in]  z      Evaluation thickness
    /// @param[in]  out    (not used)
    ///
    /// @return     { description_of_the_return_value }
    ///
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override
    {return m_materialMat->eval3D_CauchyPStress(patch,u,z,out);}

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_theta(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_gamma(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstretch(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override
    { return m_materialMat->eval3D_pstretch( patch,u,z ); }

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstretchDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override
    { return m_materialMat->eval3D_pstretchDir( patch,u,z ); }

    /// See \ref gsMaterialMatrixBase for details
    void setParameters(const std::vector<gsFunctionSet<T>*> &pars)
    { m_materialMat->setParameters(pars); }

    /// See \ref gsMaterialMatrixBase for details
    void info() const
    { m_materialMat->info(); }

    /// See \ref gsMaterialMatrixBase for details
    void setUndeformed(function_ptr undeformed) override
    {
        Base::setUndeformed(undeformed);
        m_materialMat->setUndeformed(m_patches);
    }

    /// See \ref gsMaterialMatrixBase for details
    void setDeformed(function_ptr deformed) override
    {
        Base::setDeformed(deformed);
        m_materialMat->setDeformed(m_defpatches);
    }

    /// See \ref gsMaterialMatrixBase for details
    void setUndeformed(const gsFunctionSet<T> * undeformed) override
    {
        Base::setUndeformed(undeformed);
        m_materialMat->setUndeformed(m_patches);
    }

    /// See \ref gsMaterialMatrixBase for details
    void setDeformed(const gsFunctionSet<T> * deformed) override
    {
        Base::setDeformed(deformed);
        m_materialMat->setDeformed(m_defpatches);
    }

    /// Updates the reference to the deformed patches used for TFT (explicit only)
    void updateDeformed(const gsFunctionSet<T> * deformed)
    {
        function_ptr fun = memory::make_shared_not_owned(deformed);
        m_defpatches0 = fun;
    }

    /// Updates the reference to the deformed patches used for TFT (explicit only)
    void updateDeformed(const function_ptr & deformed)
    {
        m_defpatches0 = deformed;
    }

    /// Computes theta
    gsMatrix<T> eval_theta(const gsMatrix<T> & Cs, const gsMatrix<T> & Ns, const gsMatrix<T> & Es) const;

    /// See \ref gsMaterialMatrixBase for details
    std::ostream &print(std::ostream &os) const override;

    /// See \ref gsMaterialMatrixBase for details
    const gsMaterialMatrixBase<T> * material() const override { return m_materialMat.get(); }

    /// See \ref gsMaterialMatrixBase for details
    gsMaterialMatrixBase<T> * material() override { return m_materialMat.get(); }

protected:
    template <bool _linear>
    typename std::enable_if< _linear, gsMatrix<T> >::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    template <bool _linear>
    typename std::enable_if<!_linear, gsMatrix<T> >::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    gsMatrix<T> _compute_TF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;
    gsMatrix<T> _compute_TF(const index_t patch, const gsVector<T> & u, const T & z)           const;

    gsMatrix<T> _compute_E(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & E) const;

    gsMatrix<T> _compute_S(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S) const;

    gsMatrix<T> _compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S) const;

    gsMatrix<T> _compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & dC) const;


    T _compute_gamma(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;
    gsVector<T> _theta_interval(const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;
    bool _check_theta_full(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;
    bool _check_theta_gamma(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;

protected:

    mutable material_ptr m_materialMat;
    using Base::m_patches;
    using Base::m_defpatches;
    function_ptr m_defpatches0;
    using Base::m_options;

    using Base::m_data;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixTFT.hpp)
#endif
