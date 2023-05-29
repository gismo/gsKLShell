/** @file gsMaterialMatrixBase.h

    @brief Provides a base class for material matrices

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
#include <gsKLShell/gsMaterialMatrixUtils.h>

namespace gismo
{

/**
 * @brief      This class defines the base class for material matrices
 *
 * @tparam     T     Real type
 *
 * @ingroup    KLShell
 */
template <class T>
class gsMaterialMatrixBase
{
public:

    enum {Linear=0};

    /// Shared pointer for gsGeometry
    typedef memory::shared_ptr< gsMaterialMatrixBase > Ptr;

    /// Unique pointer for gsGeometry
    typedef memory::unique_ptr< gsMaterialMatrixBase > uPtr;

    gsMaterialMatrixBase()
    :
    m_defpatches(nullptr),
    m_patches(nullptr)
    {}
    
    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsMaterialMatrixBase, clone)

    /// Destructor
    virtual ~gsMaterialMatrixBase() {};

    /**
     * @brief      Specifies how the matrix is integrated
     *
     * Possible options to specify the integration of the matrix (see enum MatIntegration)
     *  - NotIntegrated     The matrix is not integrated, hence the values of eval3D_matrix provide the through-thickness value for each z coordinate
     *  - Integrated        The matrix is already integrated; the integrated values are returned for any z coordinate
     *  - Constant          The matrix is constant over the thickness, but not yet integrated
     *  - Linear            The matrix is linearly dependent of the thickness value, say a*z+b, hence the 0th moment w.r.t. z is a*t, the 1st moment is 1/2*b*t^2, etc.
     *
     * @return     Return MatIntegration enumeration
     */
    inline virtual enum MatIntegration isMatIntegrated() const {return MatIntegration::NotIntegrated; }

    /**
     * @brief      Specifies how the vector is integrated
     *
     * Possible options to specify the integration of the vector (see enum MatIntegration)
     *  - NotIntegrated     The vector is not integrated, hence the values of eval3D_matrix provide the through-thickness value for each z coordinate
     *  - Integrated        The vector is already integrated; the integrated values are returned for any z coordinate
     *  - Constant          The vector is constant over the thickness, but not yet integrated
     *  - Linear            The vector is linearly dependent of the thickness value, say a*z+b, hence the 0th moment w.r.t. z is a*t, the 1st moment is 1/2*b*t^2, etc.
     *
     * @return     Return MatIntegration enumeration
     */
    inline virtual enum MatIntegration isVecIntegrated() const {return MatIntegration::NotIntegrated; }

    /**
     * @brief      Returns the options
     *
     * @return     \ref gsOptionList
     */
    inline virtual gsOptionList & options()
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Sets the options
     *
     * @param[in]  opt   @ref gsOptionList
     */
    inline virtual void setOptions(gsOptionList opt)
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the density multiplied by the thickness of the shell (scalar)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    inline virtual void    density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the stretches in the shell (3x1 vector)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    inline virtual void    stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the directions of the stretches in the shell (3x1 vector per direction)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    inline virtual void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the thickness of the shell (scalar)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    inline virtual void  thickness_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the parameters of the shell
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    inline virtual void  parameters_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      todo
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */    
    inline virtual void  spec2cov_transform_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      todo
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */    
    inline virtual void  spec2con_transform_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      todo
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */    
    inline virtual void  cov2cart_transform_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      todo
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The resut
     */
    inline virtual void  con2cart_transform_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    inline virtual void  transform_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Computes the deformation tensor C = F'F
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    inline virtual void  deformation_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the matrix on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (matrix) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    inline virtual gsMatrix<T>  eval3D_matrix(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }

    inline virtual gsMatrix<T>  eval3D_matrix(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_matrix(patch,u,zmat,out);
    }

    inline virtual gsMatrix<T>  eval3D_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T>& u, const T z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the derivative of the matrix on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (matrix) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     *             every column has 9*3=27 entries. To obtain the derivatives of C, one can do the following
     *
     *             dC = eval.reshapeCol(k,9,3);
     *             dCd11 = dC.reshapeCol(0,3,3);
     *             dCd22 = dC.reshapeCol(1,3,3);
     *             dCd11 = dC.reshapeCol(2,3,3);
     *
     */
    inline virtual gsMatrix<T>  eval3D_dmatrix(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    {
        GISMO_NO_IMPLEMENTATION;
        // return gsMatrix<T>::Zero(27,u.cols()*z.rows());
    }

    inline virtual gsMatrix<T>  eval3D_dmatrix(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_dmatrix(patch,u,zmat,out);
    }

    // /**
    //  * @brief      { function_description }
    //  *
    //  * @param[in]  patch  The patch
    //  * @param[in]  u      { parameter_description }
    //  * @param[in]  z      { parameter_description }
    //  * @param[in]  out    The out
    //  *
    //  * @return     { description_of_the_return_value }
    //  */
    // inline virtual gsMatrix<T>  eval3D_dmatrix(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    // {
    //     GISMO_NO_IMPLEMENTATION;
    //     // return gsMatrix<T>::Zero(27,u.cols()*z.rows());
    // }

    // inline virtual gsMatrix<T>  eval3D_dmatrix(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    // {
    //     gsMatrix<T> zmat(1,1);
    //     zmat<<z;
    //     return eval3D_dmatrix(patch,u,zmat,out);
    // }

    /**
     * @brief      Evaluates the vector on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (vector) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    inline virtual gsMatrix<T>  eval3D_vector(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }
    inline virtual gsMatrix<T>  eval3D_vector(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_vector(patch,u,zmat,out);
    }
    inline virtual gsMatrix<T>  eval3D_vector_C(const gsMatrix<T> & C, const index_t patch, const gsVector<T>& u, const T z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the Cauchy Stress vector on \a patch on in-plane points \a u with height \a z
     *
     * note: the Cauchy stress vector is returned in the actual basis
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (vector) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    inline virtual gsMatrix<T>  eval3D_CauchyVector(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    {
        GISMO_NO_IMPLEMENTATION;
    }
    // = 0;
    /**
     * @brief      Evaluates the principal stress on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (principal stress) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    inline virtual gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }
    inline virtual gsMatrix<T> eval3D_CauchyPStress(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the principal strain on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (principal strain) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    inline virtual gsMatrix<T> eval3D_pstrain(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }

    /// See \ref gsMaterialMatrixBase for details
    inline virtual gsMatrix<T> eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }
    inline virtual gsMatrix<T> eval3D_strain(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_strain(patch,u,zmat,out);
    }

    /// See \ref gsMaterialMatrixBase for details
    inline virtual gsMatrix<T> eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }
    inline virtual gsMatrix<T> eval3D_stress(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_stress(patch,u,zmat,out);
    }

    inline virtual gsMatrix<T> eval3D_detF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }

    /// See \ref gsMaterialMatrixBase for details
    inline virtual gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }
    inline virtual gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_CauchyStress(patch,u,zmat,out);
    }

    /**
     * @brief      to do
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (principal strain) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    inline virtual gsMatrix<T> eval3D_tensionfield(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const
    { GISMO_NO_IMPLEMENTATION; }

    virtual void setYoungsModulus(const gsFunction<T> & YoungsModulus)
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsFunction<T> * getYoungsModulus()
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setPoissonsRatio(const gsFunction<T> & PoissonsRatio)
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsFunction<T> * getPoissonsRatio()
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setRatio(const gsFunction<T> & Ratio)
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsFunction<T> * getRatio()
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setMu(const index_t & i, const gsFunction<T> & Mu_i)
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsFunction<T> * getMu(const index_t & i)
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setAlpha(const index_t & i, const gsFunction<T> & Alpha_i)
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsFunction<T> * getAlpha(const index_t & i)
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setThickness(const gsFunction<T> & thickness)
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsFunction<T> * getThickness()
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setDensity(const gsFunction<T> & Density)
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsFunction<T> * getDensity()
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    inline virtual void setParameters(const std::vector<gsFunction<T>*> &pars)
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Gets the number of parameters
     *
     */
    inline virtual index_t numParameters() const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Prints info
     */
    inline virtual void info() const
    { GISMO_NO_IMPLEMENTATION; }

    inline virtual gsMatrix<T> S(const gsMatrix<T> & strain) const
    { GISMO_NO_IMPLEMENTATION; }

    inline virtual gsMatrix<T> C(const gsMatrix<T> & strain) const
    { GISMO_NO_IMPLEMENTATION; }
    
    virtual void setUndeformed(const gsFunctionSet<T> * undeformed) {m_patches = undeformed; }
    virtual void setDeformed(const gsFunctionSet<T> * deformed) {m_defpatches = deformed; }

    const gsFunctionSet<T> & getUndeformed() { return *m_patches; }
    const gsFunctionSet<T> & getDeformed() { return *m_defpatches; }

    bool hasUndeformed() const { return m_defpatches!=nullptr; }
    bool hasDeformed() const { return m_defpatches!=nullptr; }

protected:
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;
};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsMaterialMatrixBase
   */
  void pybind11_init_gsMaterialMatrixBase(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace
