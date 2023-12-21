/** @file gsMaterialMatrixNonlinear.h

    @brief Provides hyperelastic material matrices

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
#include <gsIO/gsOptionList.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{


/**
 * @brief      This class defines hyperelastic material matrices
 *
 * @tparam     dim    The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T      Real type
 * @tparam     matId  Encoded ID for material (see \ref encodeMat_id)
 * @tparam     comp   Compressibility flag
 * @tparam     mat    Material flag (see \ref Material)
 * @tparam     imp    Implementation flag (see \ref Implementation)
 *
 * @ingroup    KLShell
 */
template <  short_t dim,
            class T,
            short_t matId,
            bool comp,
            enum Material mat = decodeMat_id<matId>::material,
            enum Implementation imp  = decodeMat_id<matId>::implementation
         >
class gsMaterialMatrix :    public gsMaterialMatrixBaseDim<dim,T>
{
public:

    GISMO_CLONE_FUNCTION(gsMaterialMatrix)

    using Base = gsMaterialMatrixBaseDim<dim,T>;

    typedef typename Base::function_ptr function_ptr;

    /**
     * @brief      Constructor without material parameters
     *
     * @param[in]  mp             Original geometry
     * @param[in]  thickness      Thickness function
     */
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness);

    /**
     * @brief      General constructor without density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars);

    /**
     * @brief      General constructor without density and multipatch
     *
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrix(   const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars);

    /**
     * @brief      General constructor without multipatch
     *
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrix(   const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars,
                        const gsFunctionSet<T> & Density);

    /**
     * @brief      Full general constructor
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars,
                        const gsFunctionSet<T> & Density);

protected:
    /**
     * @brief      Full general constructor
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrix(   const gsFunctionSet<T> * mp,
                        const gsFunctionSet<T> * thickness,
                        const std::vector<gsFunctionSet<T> *> &pars,
                        const gsFunctionSet<T> * Density);

public:
    /// Destructor
    gsMaterialMatrix() { }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isMatIntegrated() const {return MatIntegration::NotIntegrated; }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isVecIntegrated() const {return MatIntegration::NotIntegrated; }

    /// See \ref gsMaterialMatrixBase for details
    void defaultOptions() override;

    /// See \ref gsMaterialMatrixBase for details
    void pstretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    /// See \ref gsMaterialMatrixBase for details
    void pstretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_dmatrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyVector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstressDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_stress_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_detF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// Sets the YoungsModulus
    void setYoungsModulus(const gsFunctionSet<T> & YoungsModulus);

    /// Gets the YoungsModulus
    const function_ptr getYoungsModulus() const ;

    /// Sets the Poisson's Ratio
    void setPoissonsRatio(const gsFunctionSet<T> & PoissonsRatio);
    /// Gets the Poisson's Ratio
    const function_ptr getPoissonsRatio() const ;

    /// Sets the Ratio for the MR material
    void setRatio(const gsFunctionSet<T> & Ratio) { _setRatio_impl<mat>(Ratio);}
    /// Gets the Ratio for the MR material
    const function_ptr getRatio() const { return _getRatio_impl<mat>(); }

    /// Sets Mu_i
    void setMu(const index_t & i, const gsFunctionSet<T> & Mu_i) { _setMu_impl<mat>(i,Mu_i); }

    /// Gets Mu_i
    const function_ptr getMu(const index_t & i) const {return _getMu_impl<mat>(i);}

    /// Sets Alpha_i
    void setAlpha(const index_t & i, const gsFunctionSet<T> & Alpha_i) { _setAlpha_impl<mat>(i,Alpha_i); }

    /// Gets Alpha_i
    const function_ptr getAlpha(const index_t & i) const {return _getAlpha_impl<mat>(i);}

    /// See \ref gsMaterialMatrixBase for details
    std::ostream &print(std::ostream &os) const override;

public:
    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

protected:
    /**
     * @brief      Initializes the object.
     *
     * Initializes options, flags and defines the number of parameters
     *
     */
    void _initialize();

private:
    /**
     * @brief      Implementation of pstretch_into, specialization for compressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<_com, void>::type _pstretch_into_impl(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /**
     * @brief      Implementation of pstretch_into, specialization for incompressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<!_com, void>::type _pstretch_into_impl(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /**
     * @brief      Implementation of stretchDir_into, specialization for compressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<_com, void>::type _pstretchDir_into_impl(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /**
     * @brief      Implementation of stretchDir_into, specialization for incompressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<!_com, void>::type _pstretchDir_into_impl(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

protected:

    /**
     * @brief      Evalluates the incompressible material matrix
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The matrices (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval3D_Incompressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval3D_Incompressible_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;

    /**
     * @brief      Evalluates the incompressible stress vector
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The vectors (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval3D_Incompressible_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval3D_Incompressible_stress_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;

    gsMatrix<T> _eval3D_Incompressible_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /// Equivalent to  \ref _eval_Incompressible_vector but for the Cauchy Vector
    gsMatrix<T> _eval_Incompressible_CauchyVector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evaluates the C33 (through thickness) component of the deformation tensor C for compressible materials
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     C33 values (every column for a point u*z)
     */
    gsMatrix<T> _eval3D_Compressible_C33(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval3D_Compressible_C33(const gsMatrix<T> & Cmat, const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evaluates the jacobian determinant for compressible materials
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     { description_of_the_return_value }
     */
    gsMatrix<T> _eval3D_Compressible_detF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evaluates the jacobian determinant for compressible materials
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     { description_of_the_return_value }
     */
    gsMatrix<T> _eval3D_Incompressible_detF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /// Equivalent to  \ref _eval_Compressible_vector but for the Cauchy Vector
    gsMatrix<T> _eval_Compressible_CauchyVector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evalluates the compressible material matrix
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The matrices (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval3D_Compressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval3D_Compressible_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;

    /**
     * @brief      Evalluates the compressible stress vector
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The vectors (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval3D_Compressible_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval3D_Compressible_stress_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;

    gsMatrix<T> _eval3D_Compressible_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    // TODO: Add docs and add implementations to private member functions
    constexpr gsMatrix<T> dCijkl(const index_t patch, const gsVector<T> & u, const T z) const;

    constexpr T dCijkl_dCmn(const index_t i, const index_t j, const index_t k, const index_t l, const index_t m, const index_t n) const;

    constexpr T _dgconij_dCkl(const gsMatrix<T> & gcon, const index_t i, const index_t j, const index_t k, const index_t l) const;

private:

    template <enum Material _mat, bool _comp>
    constexpr typename std::enable_if< (!_comp && _mat==Material::NH), gsMatrix<T>>::type dCijkl_impl(const index_t patch, const gsVector<T> & u, const T z) const;

    template <enum Material _mat, bool _comp>
    constexpr typename std::enable_if<!(!_comp && _mat==Material::NH), gsMatrix<T>>::type dCijkl_impl(const index_t patch, const gsVector<T> & u, const T z) const;

    template <enum Material _mat, bool _comp>
    constexpr typename std::enable_if< (!_comp && _mat==Material::NH), T>::type dCijkl_dCmn_impl(const index_t i, const index_t j, const index_t k, const index_t l, const index_t m, const index_t n) const;

    template <enum Material _mat, bool _comp>
    constexpr typename std::enable_if<!(!_comp && _mat==Material::NH), T>::type dCijkl_dCmn_impl(const index_t i, const index_t j, const index_t k, const index_t l, const index_t m, const index_t n) const;

private:

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the material matrix, specialization for SvK material (incompressible)
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The material matrix (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_matrix_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the material matrix, specialization compressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The material matrix (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the material matrix, specialization for incompressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The material matrix (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;


    /**
     * @brief      { function_description }
     *
     * @param[in]  patch  The patch
     * @param[in]  u      { parameter_description }
     * @param[in]  z      { parameter_description }
     * @param[in]  out    The out
     *
     * @tparam     _mat   { description }
     * @tparam     _comp  { description }
     *
     * @return     { description_of_the_return_value }
     */
    template <enum Material _mat, bool _comp>
    typename std::enable_if<!(!_comp && _mat==Material::NH), gsMatrix<T>>::type _eval3D_dmatrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const;

    template <enum Material _mat, bool _comp>
    typename std::enable_if< (!_comp && _mat==Material::NH), gsMatrix<T>>::type _eval3D_dmatrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const;

    /// Same as \ref _eval3D_vector_impl for the Cauchy stress
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_CauchyVector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of jacobina determinant (detF)
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The stress vector (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_detF_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_detF_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_detF_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /// Same as \ref _eval3D_vector_impl for the Cauchy stress
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_CauchyVector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for SvK material (incompressible)
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The stress vector (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_stress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_stress_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_CauchyStress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /// Same as \ref _eval3D_vector_impl for the Cauchy stress
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_CauchyVector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for incompressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The stress vector (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_stress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_stress_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z)const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_CauchyStress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for compressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The stress vector (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_stress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_stress_C_impl(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_CauchyStress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for compressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The principal stresses (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_pstress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, void>::type _setRatio_impl(const gsFunctionSet<T> & Ratio);
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::MR, void>::type _setRatio_impl(const gsFunctionSet<T> & Ratio);

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, const function_ptr>::type _getRatio_impl() const;
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::MR, const function_ptr>::type _getRatio_impl() const;

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::OG, void>::type _setMu_impl(const index_t & i, const gsFunctionSet<T> & Mu_i);
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::OG, void>::type _setMu_impl(const index_t & i, const gsFunctionSet<T> & Mu_i);

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::OG, const function_ptr>::type _getMu_impl(const index_t & i) const;
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::OG, const function_ptr>::type _getMu_impl(const index_t & i) const;

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::OG, void>::type _setAlpha_impl(const index_t & i, const gsFunctionSet<T> & Alpha_i);
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::OG, void>::type _setAlpha_impl(const index_t & i, const gsFunctionSet<T> & Alpha_i);

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::OG, const function_ptr>::type _getAlpha_impl(const index_t & i) const;
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::OG, const function_ptr>::type _getAlpha_impl(const index_t & i) const;

protected:

    /**
     * @brief      Returns an entry of the material tensor C for incompressible materials
     *
     * @param[in]  i,j,k,l     The indices
     *
     * @return     C^{ijkl}
     */
    constexpr T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;

    /**
     * @brief      Returns an entry of the material tensor C for compressible materials without static condensation
     *
     * @param[in]  i,j,k,l     The indices
     * @param[in]  c           The deformation tensor
     * @param[in]  cinv        The inverse of the deformation tensor
     *
     * @return     C^{ijkl}
     */
    constexpr T _Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Returns an entry of the material tensor C for compressible materials, with static condensation included
     *
     * @param[in]  i,j,k,l     The indices
     * @param[in]  c           The deformation tensor
     * @param[in]  cinv        The inverse of the deformation tensor
     *
     * @return     C^{ijkl}
     */
    constexpr T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Returns an entry of the stress tensor S for incompressible materials
     *
     * @param[in]  i,j     The indices
     *
     * @return     S^{ij}
     */
    constexpr T _Sij    (const index_t i, const index_t j) const;

    /**
     * @brief      Returns an entry of the stress tensor S for compressible materials
     *
     * @param[in]  i,j     The indices
     * @param[in]  c       The deformation tensor
     * @param[in]  cinv    The inverse of the deformation tensor
     *
     * @return     S^{ij}
     */
    constexpr T _Sij    (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Returns an entry of the diagonal of the stress tensor S for incompressible materials
     *
     * @param[in]  i       The index
     *
     * @return     S^{ii}
     */
    constexpr T _Sii    (const index_t i) const;

    /**
     * @brief      Returns an entry of the diagonal of the stress tensor S for incompressible materials
     *
     * @param[in]  i       The index
     * @param[in]  c       The deformation tensor
     *
     * @return     S^{ii}
     */
    constexpr T _Sii    (const index_t i, const gsMatrix<T> & c) const;

private:
    ///////////////////////////////////////////////////////
    // Incompressible formulations                       //
    ///////////////////////////////////////////////////////
    /// Specialization for incompressible Cijkl(i,j,k,l) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH  && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::MR  && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Cijkl(i,j,k,l) for Extended NH materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Cijkl(i,j,k,l) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;


    ///////////////////////////////////////////////////////
    // Compressible formulations                         //
    ///////////////////////////////////////////////////////
    /**
     * @brief      Specialization for compressible Cijkl(i,j,k,l,c,cinv) for all materials implemented generally
     *
     * @note        The specialization is required because static condensation for spectrally implemented material models is performed before the transformation. For all other implementation, static condensation is performed using the Cijkl3D function.
     *
     */
    template<enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Spectral,   T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl(i,j,k,l,c,cinv) for all materials implemented not spectrally
    template<enum Implementation _imp>
    constexpr typename std::enable_if<!(_imp==Implementation::Spectral),T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for Extended NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;


    /// Specialization for incompressible Sij(i,j) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Sij(i,j) for Extended NH materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Sij(i,j) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Sij_impl(const index_t i, const index_t j) const;


    /// Specialization for compressible Sij(i,j,c,cinv) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for compressible Sij(i,j,c,cinv) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for compressible Sij(i,j,c,cinv) for Extended NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    constexpr typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

protected:
    /**
     * @brief      Provides the derivative of the incompressible strain energy density function w.r.t. component C_{ij} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C
     *
     * @return     The derivative of psi w.r.t. C_{ij}
     */
    constexpr T _dPsi   (const index_t i, const index_t j) const;

    /**
     * @brief      Provides the derivative of the compressible strain energy density function w.r.t. component C_{ij} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The derivative of psi w.r.t. C_{ij}
     */
    constexpr T _dPsi   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the derivative of the volumetric part of the compressible strain energy density function w.r.t. component C_{ij} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The derivative of psi w.r.t. C_{ij}
     */
    constexpr T _dPsi_vol(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the second (mixed) derivative of the incompressible strain energy density function w.r.t. components C_{ij} and C_{kl} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C_{ij}
     * @param[in]  k,l   The indices of C_{kl}
     *
     * @return     The second (mixed) derivative of psi w.r.t. C_{ij} and C_{kl}
     */
    constexpr T _d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l) const;

    /**
     * @brief      Provides the second (mixed) derivative of the compressible strain energy density function w.r.t. components C_{ij} and C_{kl} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C_{ij}
     * @param[in]  k,l   The indices of C_{kl}
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The second (mixed) derivative of psi w.r.t. C_{ij} and C_{kl}
     */
    constexpr T _d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the second (mixed) derivative of the volumetric part of the compressible strain energy density function w.r.t. components C_{ij} and C_{kl} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C_{ij}
     * @param[in]  k,l   The indices of C_{kl}
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The second (mixed) derivative of psi w.r.t. C_{ij} and C_{kl}
     */
    constexpr T _d2Psi_vol(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the derivative of the first invariant compressible materials w.r.t. component C_{ij} of the deformation tensor
     *
     * @note       \f$I_1=\trace{\mathbf{C}}\f$
     *
     * @param[in]  i,j   The indices of C
     *
     * @return     The derivative the first invariant w.r.t. C_{ij}
     */
    constexpr T _dI_1   (const index_t i, const index_t j) const;

    /**
     * @brief      Provides the derivative of the second invariant for compressible materials w.r.t. component C_{ij} of the deformation tensor
     *
     * @note       \f$I_1=\frac{1}{2}\left( \trace{\mathbf{C}^2} + (\trace{\mathbf{C}})^2 \right)\f$
     *
     * @param[in]  i,j   The indices of C
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The derivative the second invariant w.r.t. C_{ij}
     */
    constexpr T _dI_2   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;



private:
    /*
        Note: when adding a new implementation, make sure to add the condition in the 'other' implementation : !(_mat==xx || ...)
    */

    /// Implementation of _dPsi(i,j) for NH materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::NH, T>::type _dPsi_impl(const index_t i, const index_t j) const;

    /// Implementation of _dPsi(i,j) for MR materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::MR, T>::type _dPsi_impl(const index_t i, const index_t j) const;
    // other
    template<enum Material _mat>
    constexpr typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR), T>::type _dPsi_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};
    // add other

    /// Implementation of _dPsi(i,j,c,cinv) for NH materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::NH, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _dPsi(i,j,c,cinv) for MR materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::MR, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _dPsi(i,j,c,cinv) for Extended NH materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::NH_ext, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // other
    template<enum Material _mat>
    constexpr typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR || _mat==Material::NH_ext), T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Implementation of _d2Psi(i,j) for NH materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::NH, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Implementation of _d2Psi(i,j) for MR materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::MR, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    // other
    template<enum Material _mat>
    constexpr typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR), T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};

    /// Implementation of _d2Psi(i,j,c,cinv) for NH materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::NH, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _d2Psi(i,j,c,cinv) for MR materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::MR, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _d2Psi(i,j,c,cinv) for Extended NH materials
    template<enum Material _mat>
    constexpr typename std::enable_if<_mat==Material::NH_ext, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // other
    template<enum Material _mat>
    constexpr typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR || _mat==Material::NH_ext), T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};
protected:
    ///////////////////////////////////////////////////////
    // Stretch based formulation                         //
    ///////////////////////////////////////////////////////

    /**
     * @brief      First derivative of a strain energy density function \f$\Psi\f$ w.r.t. the stretch \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     *
     */
    constexpr T _dPsi_da   (const index_t a) const;

    /**
     * @brief      Second derivative of a strain energy density function \f$\Psi\f$ w.r.t. the stretches \f$\lambda_a\f$ and \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     *
     */
    constexpr T _d2Psi_dab (const index_t a, const index_t b) const;

    /**
     * @brief      First derivative of the volumetric part of a strain energy density function \f$\Psi_{vol}\f$ w.r.t. the stretch \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     */
    constexpr T _dPsi_da_vol(const index_t a) const;

    /**
     * @brief      Second derivative of the volumetric part of a strain energy density function \f$\Psi_{vol}\f$ w.r.t. the stretches \f$\lambda_a\f$ and \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     */
    constexpr T _d2Psi_dab_vol(const index_t a, const index_t b) const;

    /**
     * @brief      First derivative of the compressibilty function \f$\J\f$ w.r.t. the stretche \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     */
    constexpr T _dJ_da     (const index_t a) const;

    /**
     * @brief      First derivative of the compressibilty function \f$\J\f$ w.r.t. the stretches \f$\lambda_a\f$ and \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     */
    constexpr T _d2J_dab   (const index_t a, const index_t b) const;

    /**
     * @brief      Lagrange multiplier for incompressible materials
     */
    constexpr T _p()                                          const;

    /**
     * @brief      First derivative of the Lagrange multiplier for incompressible materials w.r.t. the stretch \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     */
    constexpr T _dp_da     (const index_t a) const;

    /**
     * @brief     Component \f$a\f$ of the stress
     *
     * @param[in]  a     The index a
     */
    constexpr T _Sa        (const index_t a) const;

    /**
     * @brief     First derivative of the \f$a^{\text{th}\f$ component of the stress w.r.t. the stretch \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     */
    constexpr T _dSa_db    (const index_t a, const index_t b) const;

    /**
     * @brief     The material matrix for stretch-based implementations
     *
     * @param[in]  a,b,c,d     The indices a,b,c,d
     */
    constexpr T _Cabcd     (const index_t a, const index_t b, const index_t c, const index_t d) const;

private:
    // ----------------------------------------------------------------------------------

    /// Specialization of _dPsi_da(a) for compressible NH materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::NH), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible NH materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::NH), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for compressible MR materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::MR), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible MR materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::MR), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for compressible OG materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::OG), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible OG materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::OG), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for compressible Extended NH materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::NH_ext), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible Extended NH materials (not implemented)
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::NH_ext), T>::type _dPsi_da_impl(const index_t a) const
    {GISMO_NO_IMPLEMENTATION};

    // other
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<
                            !(     _mat==Material::NH
                                || _mat==Material::MR
                                || _mat==Material::OG
                                || _mat==Material::NH_ext
                              )
                                                                , T>::type _dPsi_da_impl(const index_t a) const
    {GISMO_NO_IMPLEMENTATION}

    // ----------------------------------------------------------------------------------

    /// Specialization of _d2Psi_dab(a,b) for compressible NH materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::NH), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible NH materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::NH), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for compressible MR materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::MR), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible MR materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::MR), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for compressible OG materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::OG), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible OG materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::OG), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for compressible Extended NH materials
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<_com && (_mat==Material::NH_ext), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible Extended NH materials (not implemented)
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<!_com && (_mat==Material::NH_ext), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const
    {GISMO_NO_IMPLEMENTATION};

    // other
    template<enum Material _mat, bool _com>
    constexpr typename std::enable_if<
                            !(     _mat==Material::NH
                                || _mat==Material::MR
                                || _mat==Material::OG
                                || _mat==Material::NH_ext
                              )
                                                                , T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const
    {GISMO_NO_IMPLEMENTATION}

    // ----------------------------------------------------------------------------------

    /// Specialization of _Sa(a) for compressible materials
    template<bool _com>
    constexpr typename std::enable_if<_com , T>::type _Sa_impl(const index_t a) const;

    /// Specialization of _Sa(a) for incompressible materials
    template<bool _com>
    constexpr typename std::enable_if<!_com, T>::type _Sa_impl(const index_t a) const;

    // ----------------------------------------------------------------------------------

    /// Specialization of _dSa_db(a,b) for compressible materials
    template<bool _com>
    constexpr typename std::enable_if<_com , T>::type _dSa_db_impl(const index_t a, const index_t b) const;

    /// Specialization of _dSa_db(a,b) for incompressible materials
    template<bool _com>
    constexpr typename std::enable_if<!_com, T>::type _dSa_db_impl(const index_t a, const index_t b) const;

    // ----------------------------------------------------------------------------------

    /// Specialization of _Cabcd(a,b,c,d) for compressible materials
    template<bool _com>
    constexpr typename std::enable_if<_com , T>::type _Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;

    /// Specialization of _Cabcd(a,b,c,d) for incompressible materials
    template<bool _com>
    constexpr typename std::enable_if<!_com, T>::type _Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;

    // ----------------------------------------------------------------------------------

protected:
    using Base::m_thickness;
    using Base::m_pars;
    using Base::m_density;

    // Geometric data
    using Base::m_data;

    using Base::m_options;

private:
    static index_t delta(const index_t a, const index_t b)
    {
        return (a==b) ? 1 : 0;
    }

    static index_t idelta(const index_t a, const index_t b)
    {
        return (a!=b) ? 1 : 0;
    }
};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsMaterialMatrix
   */
  void pybind11_init_gsMaterialMatrixNH2i(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixNH2c(pybind11::module &m);

  void pybind11_init_gsMaterialMatrixNH3i(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixNH3c(pybind11::module &m);

  void pybind11_init_gsMaterialMatrixMR2i(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixMR2c(pybind11::module &m);

  void pybind11_init_gsMaterialMatrixMR3i(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixMR3c(pybind11::module &m);

  void pybind11_init_gsMaterialMatrixOG2i(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixOG2c(pybind11::module &m);

  void pybind11_init_gsMaterialMatrixOG3i(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixOG3c(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixNonlinear.hpp)
#endif
