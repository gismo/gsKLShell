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

#include <gsCore/gsFunctionSet.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>

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

    typedef typename gsFunctionSet<T>::Ptr function_ptr;

    enum {Linear=0};

    /// Shared pointer for gsGeometry
    typedef memory::shared_ptr< gsMaterialMatrixBase > Ptr;

    /// Unique pointer for gsGeometry
    typedef memory::unique_ptr< gsMaterialMatrixBase > uPtr;

    gsMaterialMatrixBase()
    :
    m_patches(nullptr),
    m_defpatches(nullptr),
    m_thickness(nullptr),
    m_density(nullptr)
    {
    }

    gsMaterialMatrixBase(   const gsFunctionSet<T> & mp,
                            const gsFunctionSet<T> & mp_def,
                            const gsFunctionSet<T> & thickness,
                            const gsFunctionSet<T> & Density)
    :
    m_patches(memory::make_shared(mp.clone().release())),
    m_defpatches(memory::make_shared(mp_def.clone().release())),
    m_thickness(memory::make_shared(thickness.clone().release())),
    m_density(memory::make_shared(Density.clone().release()))
    {
    }

    gsMaterialMatrixBase(   const gsFunctionSet<T> * mp,
                            const gsFunctionSet<T> * mp_def,
                            const gsFunctionSet<T> * thickness,
                            const gsFunctionSet<T> * Density)
    :
    m_patches(memory::make_shared_not_owned(mp)),
    m_defpatches(memory::make_shared_not_owned(mp_def)),
    m_thickness(memory::make_shared_not_owned(thickness)),
    m_density(memory::make_shared_not_owned(Density))
    {
    }

    gsMaterialMatrixBase(   const function_ptr & mp,
                            const function_ptr & mp_def,
                            const function_ptr & thickness,
                            const function_ptr & Density)
    :
    m_patches(mp),
    m_defpatches(mp_def),
    m_thickness(thickness),
    m_density(Density)
    {
    }

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsMaterialMatrixBase, clone)

    /// Destructor
    virtual ~gsMaterialMatrixBase() {};

    /// Copy constructor (makes deep copy)
    gsMaterialMatrixBase( const gsMaterialMatrixBase<T> & other )
    {
        operator=(other);
    }

    /// Move constructor
    gsMaterialMatrixBase( gsMaterialMatrixBase<T> && other )
    {
        operator=(give(other));
    }

    gsMaterialMatrixBase<T> & operator= ( const gsMaterialMatrixBase<T> & other )
    {
        if (this!=&other)
        {
            m_patches = other.m_patches;
            m_defpatches = other.m_defpatches;
            m_options = other.m_options;

            m_pars = other.m_pars;
            m_thickness = other.m_thickness;
            m_density = other.m_density;
        }
        return *this;
    }

    gsMaterialMatrixBase<T> & operator= ( gsMaterialMatrixBase<T> && other )
    {
        m_patches = give(other.m_patches);
        m_defpatches = give(other.m_defpatches);
        m_options = give(other.m_options);


        m_pars = give(other.m_pars);
        m_thickness = give(other.m_thickness);
        m_density = give(other.m_density);
        return *this;
    }

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
    virtual enum MatIntegration isMatIntegrated() const {return MatIntegration::NotIntegrated; }

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
    virtual enum MatIntegration isVecIntegrated() const {return MatIntegration::NotIntegrated; }

    /**
     * @brief      Returns the options
     *
     * @return     \ref gsOptionList
     */
    virtual void defaultOptions()
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Returns the options
     *
     * @return     \ref gsOptionList
     */
    virtual gsOptionList & options() { return m_options; }
    /**
     * @brief      Sets the options
     *
     * @param[in]  opt   @ref gsOptionList
     */
    virtual void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    /**
     * @brief      Evaluates the density multiplied by the thickness of the shell (scalar)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    virtual void    density_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the stretches in the shell (3x1 vector)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    GISMO_DEPRECATED
    virtual void    pstretch_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the directions of the stretches in the shell (3x1 vector per direction)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    GISMO_DEPRECATED
    virtual void    pstretchDir_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }
        /**
     * @brief      Evaluates the priciple stresses in the shell (3x1 vector)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    GISMO_DEPRECATED
    virtual void    pstress_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the directions of the priciple stresses in the shell (3x1 vector per direction)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    GISMO_DEPRECATED
    virtual void    pstressDir_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the thickness of the shell (scalar)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    virtual void  thickness_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the parameters of the shell
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    GISMO_DEPRECATED
    virtual void  parameters_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Constructs a transformation matrix that transforms a quantity (IN VOIGHT NOTATION) in the spectral basis to the (undeformed) convariant basis
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The result
     */
    GISMO_DEPRECATED
    virtual void  transform_into(const index_t /*patch*/, const gsMatrix<T>& /*u*/, gsMatrix<T>& /*result*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Provides the transformation from the spectral basis to the
     *             covariant basis
     *
     * @param[in]  patch  The patch to be evaluated on
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param      result  The result
     *
     * @return     Matrix with the result (deformation tensor) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T>  eval3D_spec2cov(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Provides the transformation from the spectral basis to the
     *             contravariant basis
     *
     * @param[in]  patch  The patch to be evaluated on
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     *
     * @return     Matrix with the result (deformation tensor) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T>  eval3D_spec2con(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Provides the transformation from the covariant basis to the
     *             local cartesian basis
     *
     * @param[in]  patch  The patch to be evaluated on
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param      result  The result
     *
     * @return     Matrix with the result (deformation tensor) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_cov2cart(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Provides the transformation from the contravariant basis to
     *             the local cartesian basis
     *
     * @param[in]  patch  The patch to be evaluated on
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     *
     * @return     Matrix with the result (deformation tensor) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_con2cart(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Provides the transformation of the principle stretches
     *
     * @param[in]  patch  The patch to be evaluated on
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     *
     * @return     Matrix with the result (deformation tensor) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_pstretchTransform(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
    { return eval3D_spec2con(patch,u,z); }

    /**
     * @brief      Provides the transformation of the principle stresses
     *
     * @param[in]  patch  The patch to be evaluated on
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     *
     * @return     Matrix with the result (deformation tensor) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_pstressTransform(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
    { return eval3D_spec2cov(patch,u,z); }

    /**
     * @brief      Evaluates deformation tensor on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (deformation tensor) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T>  eval3D_deformation(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the matrix on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @note       For the matrix outputs (MatrixA, MatrixB, MatrixC, MatrixD), each column holds a
     *             3x3 matrix in column-major order (9 entries) with M_ij = dS_i/dE_j: the row index is
     *             the stress (Voigt) component and the column index is the strain (engineering Voigt)
     *             component. The shell assemblers form the tangent as K = dE_test * M * dE_trial, which
     *             requires this orientation for non-symmetric tangents.
     *
     * @return     Matrix with the result (matrix) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T>  eval3D_matrix(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }

    virtual gsMatrix<T>  eval3D_matrix(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_matrix(patch,u,zmat,out);
    }

    virtual gsMatrix<T>  eval3D_matrix_C(const gsMatrix<T> & /*Cmat*/, const index_t /*patch*/, const gsVector<T>& /*u*/, const T /*z*/, enum MaterialOutput /*out*/) const
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
     *             every column has 27 entries. Entry d + 3*(i + 3*l) is dC_il/dx_d: the derivative of entry (i,l)
     *             of the 3x3 matrix returned by \ref eval3D_matrix with respect to x = (g_11, g_22, 2 g_12), the
     *             deformed metric in Voigt form with doubled shear entry. Hence eval.reshapeCol(k,3,9) has the
     *             derivative variable d as its row and the column-major matrix entry i+3l as its column.
     *             With respect to the engineering strain E = (E_11, E_22, 2E_12), dC/dE_d = 2 dC/dx_d.
     *             \ref gsMaterialMatrixLinear returns zeros.
     *
     */
    virtual gsMatrix<T>  eval3D_dmatrix(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    {
        GISMO_NO_IMPLEMENTATION;
        // return gsMatrix<T>::Zero(27,u.cols()*z.rows());
    }

    virtual gsMatrix<T>  eval3D_dmatrix(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_dmatrix(patch,u,zmat,out);
    }

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
    virtual gsMatrix<T>  eval3D_vector(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T>  eval3D_vector(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_vector(patch,u,zmat,out);
    }
    virtual gsMatrix<T>  eval3D_vector_C(const gsMatrix<T> & /*C*/, const index_t /*patch*/, const gsVector<T>& /*u*/, const T /*z*/, enum MaterialOutput /*out*/) const
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
    virtual gsMatrix<T>  eval3D_CauchyVector(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    {
        GISMO_NO_IMPLEMENTATION;
    }
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
    virtual gsMatrix<T> eval3D_pstress(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_CauchyPStress(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the principal stress directions on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (principal stress directions) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     *
     * @note       Layout: 9 rows, the column-major 3x3 matrix [d_min | d_max | N] of Cartesian unit vectors
     *             (entry (i,l) in row i+3*l), minor principal direction first (ascending, as \ref eval3D_pstress;
     *             Verhelst's S_p,1 is Zhang's sigma_2), N = normalise(G_1 x G_2) the unit normal of the reference
     *             (undeformed) configuration at the evaluation point, d_min and d_max Cartesian unit vectors in its
     *             tangent plane, and d_max = N x d_min. Every implementation builds the frame from the reference
     *             covariant basis.
     */
    virtual gsMatrix<T> eval3D_pstressDir(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the principal stretch on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (principal stretch) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_pstretch(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }
    /**
     * @brief      Evaluates the principal stretch directions on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (for classes with MatIntegration==Integrated, more details about \ref MaterialOutput can be found in \ref gsMaterialMatrixUtils)
     *
     * @return     Matrix with the result (principal stretch directions) ordered per z coordinate per point
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_pstretchDir(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
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
    virtual gsMatrix<T> eval3D_pstrain(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      { function_description }
     *
     * @param[in]  patch  The patch
     * @param[in]  u      { parameter_description }
     * @param[in]  z      { parameter_description }
     * @param[in]  out    The out
     *
     * @return     { description_of_the_return_value }
     */
    virtual gsMatrix<T> eval3D_strain(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_strain(const index_t patch, const gsVector<T>& u, const T & z) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_strain(patch,u,zmat);
    }

    /**
     * @brief      { function_description }
     *
     * @param[in]  patch  The patch
     * @param[in]  u      { parameter_description }
     * @param[in]  z      { parameter_description }
     * @param[in]  out    The out
     *
     * @return     { description_of_the_return_value }
     */
    virtual gsMatrix<T> eval3D_stress(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_stress(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_stress(patch,u,zmat,out);
    }

    /**
     * @brief      { function_description }
     *
     * @param[in]  patch  The patch
     * @param[in]  u      { parameter_description }
     * @param[in]  z      { parameter_description }
     * @param[in]  out    The out
     *
     * @return     { description_of_the_return_value }
     */
    virtual gsMatrix<T> eval3D_detF(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      { function_description }
     *
     * @param[in]  patch  The patch
     * @param[in]  u      { parameter_description }
     * @param[in]  z      { parameter_description }
     * @param[in]  out    The out
     *
     * @return     { description_of_the_return_value }
     */
    virtual gsMatrix<T> eval3D_CauchyStress(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
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
    virtual gsMatrix<T> eval3D_tensionfield(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_tensionfield(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_tensionfield(patch,u,zmat,out);
    }

    virtual gsMatrix<T> eval3D_theta(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_theta(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_theta(patch,u,zmat,out);
    }

    virtual gsMatrix<T> eval3D_gamma(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_gamma(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_gamma(patch,u,zmat,out);
    }

    /**
     * @brief      Evaluates the unit vector along the wrinkles (tension-field materials) on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (not used)
     *
     * @return     3 rows: the Cartesian unit vector tangent to the wrinkles on the deformed surface,
     *             normalise(-sin(theta) g_1 + cos(theta) g_2). The zero vector where the point is not wrinkled.
     *             Ordered per z coordinate per point (column j*u.cols()+k for point k, fibre j)
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_wrinkleDir(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_wrinkleDir(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_wrinkleDir(patch,u,zmat,out);
    }

    /**
     * @brief      Evaluates the wrinkling intensity (tension-field materials) on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (not used)
     *
     * @return     1 row: max(0,-E_min) with E_min the minor (ascending order, first) principal value of the total
     *             Green-Lagrange strain w.r.t. the reference metric (|E_2| of Zhang & Kiendl 2024, who number the
     *             major value first). Zero where the point is not wrinkled.
     *             Ordered per z coordinate per point (column j*u.cols()+k for point k, fibre j)
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_wrinkleIntensity(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_wrinkleIntensity(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_wrinkleIntensity(patch,u,zmat,out);
    }

    /**
     * @brief      Evaluates the parametrisation-independent wrinkle strain (tension-field materials) on \a patch on in-plane points \a u with height \a z
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     * @param[in]  z      The point through-thickness coorinate
     * @param[in]  out    (not used)
     *
     * @return     1 row: gamma * G^{ab} n_a n_b with n = (cos theta, sin theta) and G^{ab} the inverse of the
     *             2x2 covariant reference metric. NaN where the point is not wrinkled.
     *             Ordered per z coordinate per point (column j*u.cols()+k for point k, fibre j)
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_gammaHat(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_gammaHat(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_gammaHat(patch,u,zmat,out);
    }

    /**
     * @brief      Evaluates the strain energy W on \a patch on in-plane points \a u with height \a z
     *
     * W is the strain energy per unit reference volume at the fibre point z; it is not integrated through the thickness.
     * It is normalised such that S = dW/dE holds exactly, with S the second Piola-Kirchhoff stress in Voigt form
     * (S^11, S^22, S^12) as returned by \ref eval3D_vector and \ref eval3D_vector_C, and E the engineering Voigt strain
     * (E_11, E_22, 2E_12), E_ab = 1/2 (g_ab - G_ab). W(E=0) = 0.
     * For plane stress, W is evaluated at the same condensed state as the stress: C33 = 1/J0^2 for
     * incompressible laws and the C33 of the through-thickness Newton iteration for compressible laws.
     * The default raises GISMO_NO_IMPLEMENTATION; gsMaterialMatrixLinear, gsMaterialMatrixNonlinear (except
     * Material::SvK) and gsMaterialMatrixTFTSplit implement it.
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinates, a matrix of size nz x u.cols(): column k holds the nz
     *                    fibre heights of point k (entry z(j,k))
     * @param[in]  out    The output type. gsMaterialMatrixLinear selects the strain by it (membrane outputs such
     *                    as VectorN: mid-surface metric; bending outputs such as VectorM: change of curvature;
     *                    Generic, Energy, Strain, Stress or PStress: fibre metric G(z); any other value is an
     *                    error); gsMaterialMatrixNonlinear ignores it
     *
     * @return     1 row. Ordered per z coordinate per point (column j*u.cols()+k for point k, fibre j)
     *                  [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
     */
    virtual gsMatrix<T> eval3D_energy(const index_t /*patch*/, const gsMatrix<T>& /*u*/, const gsMatrix<T>& /*z*/, enum MaterialOutput /*out*/ = MaterialOutput::Generic) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual gsMatrix<T> eval3D_energy(const index_t patch, const gsVector<T>& u, const T & z, enum MaterialOutput out = MaterialOutput::Generic) const
    {
        gsMatrix<T> zmat(1,1);
        zmat<<z;
        return eval3D_energy(patch,u,zmat,out);
    }

    /**
     * @brief      Evaluates the strain energy W for a given deformed metric, at a single point
     *
     * Same definition and normalisation as \ref eval3D_energy (per unit reference volume at fibre z,
     * S = dW/dE, W(E=0) = 0, plane stress condensed state).
     *
     * @param[in]  Cmat   The deformed covariant metric, a 3x1 unhalved vector [g11, g22, g12]
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     * @param[in]  out    The output type. gsMaterialMatrixLinear selects the strain by it (membrane outputs such
     *                    as VectorN: mid-surface metric; bending outputs such as VectorM: change of curvature;
     *                    Generic, Energy, Strain, Stress or PStress: fibre metric G(z); any other value is an
     *                    error); gsMaterialMatrixNonlinear ignores it
     *
     * @return     1x1 matrix holding W
     */
    virtual gsMatrix<T> eval3D_energy_C(const gsMatrix<T> & /*Cmat*/, const index_t /*patch*/, const gsVector<T>& /*u*/, const T /*z*/, enum MaterialOutput /*out*/ = MaterialOutput::Generic) const
    { GISMO_NO_IMPLEMENTATION; }

    virtual void setYoungsModulus(const gsFunctionSet<T> & /*YoungsModulus*/)
    { GISMO_NO_IMPLEMENTATION; }
    virtual const function_ptr getYoungsModulus() const
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setPoissonsRatio(const gsFunctionSet<T> & /*PoissonsRatio*/)
    { GISMO_NO_IMPLEMENTATION; }
    virtual const function_ptr getPoissonsRatio() const
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setRatio(const gsFunctionSet<T> & /*Ratio*/)
    { GISMO_NO_IMPLEMENTATION; }
    virtual const function_ptr getRatio() const
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setMu(const index_t & /*i*/, const gsFunctionSet<T> & /*Mu_i*/)
    { GISMO_NO_IMPLEMENTATION; }
    virtual const function_ptr getMu(const index_t &) const
    { GISMO_NO_IMPLEMENTATION; }
    virtual void setAlpha(const index_t & /*i*/, const gsFunctionSet<T> & /*Alpha_i*/)
    { GISMO_NO_IMPLEMENTATION; }
    virtual const function_ptr getAlpha(const index_t & /*i*/) const
    { GISMO_NO_IMPLEMENTATION; }

    /// Sets the thickness
    virtual void setThickness(const function_ptr & thickness) { m_thickness = thickness; }

    /// Sets the thickness
    virtual void setThickness(const gsFunctionSet<T> & thickness)
    {
        function_ptr fun = memory::make_shared(thickness.clone().release());
        m_thickness = fun;
    }

    /// Returns true if a thickness is assigned
    virtual bool hasThickness() const { return m_thickness!=nullptr; }

    /// Gets the Density
    virtual const function_ptr getThickness() const {return m_thickness;}

    /// Sets the density
    virtual void setDensity(function_ptr Density) { m_density = Density; }
    /// Sets the density
    virtual void setDensity(const gsFunctionSet<T> & Density)
    {
        function_ptr fun = memory::make_shared(Density.clone().release());
        m_density = fun;
    }
    /// Returns true if a density is assigned
    virtual bool hasDensity() const { return m_density!=nullptr; }

    /// Gets the Density
    virtual const function_ptr getDensity()  const {return m_density;}

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual void setParameters(const std::vector<function_ptr> &pars)
    {
        m_pars = pars;
    }

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual void setParameter(const index_t i, const function_ptr &par)
    {
        m_pars[i] = par;
    }

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual void setParameters(const std::vector<gsFunctionSet<T> *> &pars)
    {
        m_pars.resize(pars.size());
        for (size_t k = 0; k!=pars.size(); k++)
            m_pars[k] = memory::make_shared_not_owned(pars[k]);
    }

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual void setParameter(const index_t i, const gsFunctionSet<T> &par)
    {
        if ((index_t)m_pars.size() < i+1)
            m_pars.resize(i+1);
        m_pars[i] = memory::make_shared(par.clone().release());
    }

    /**
     * @brief      Gets parameter i
     *
     * @param[in]  i     The parameter index
     *
     * @return     The parameter.
     */
    virtual const function_ptr getParameter(const index_t i)  const
    {
        GISMO_ASSERT(i < (index_t)m_pars.size(),"Parameter "<<i<<" is unavailable");
        return m_pars[i] ;
    }

    /**
     * @brief      Gets the number of parameters
     *
     */
    virtual index_t numParameters() const { return m_pars.size(); }

    /// See \ref gsMaterialMatrixBase for details
    virtual void resetParameters()
    {
        m_pars.clear();
        m_pars.resize(0);
    }

    /**
     * @brief      Prints info
     */
    virtual void info() const
    { GISMO_NO_IMPLEMENTATION; }

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os<<"gsMaterialMatrixBase (type not understood).\n";
        return os;
    }

    /// Returns this
    virtual const gsMaterialMatrixBase<T> * material() const { return this; }
    /// Returns this
    virtual gsMaterialMatrixBase<T> * material() { return this; }

    // TODO: Change from single letter to better description
    virtual gsMatrix<T> S(const gsMatrix<T> & /*strain*/) const
    { GISMO_NO_IMPLEMENTATION; }

    virtual gsMatrix<T> C(const gsMatrix<T> & /*strain*/) const
    { GISMO_NO_IMPLEMENTATION; }

    virtual void setUndeformed(const gsFunctionSet<T> * undeformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(undeformed);
        m_patches = f_ptr;
    }
    virtual void setDeformed(const gsFunctionSet<T> * deformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(deformed);
        m_defpatches = f_ptr;
    }

    virtual void setUndeformed(const function_ptr undeformed) {m_patches = undeformed; }
    virtual void setDeformed(const function_ptr deformed) {m_defpatches = deformed; }

    const function_ptr getUndeformed() const { return m_patches; }
    const function_ptr getDeformed()  const { return m_defpatches; }

    virtual bool hasUndeformed() const
    { return m_patches!=nullptr; }
    virtual bool hasDeformed() const { return m_defpatches!=nullptr; }

    virtual bool initialized() const
    { GISMO_NO_IMPLEMENTATION; }

protected:
    function_ptr m_patches;
    function_ptr m_defpatches;
    gsOptionList m_options;

    std::vector< function_ptr > m_pars;
    function_ptr m_thickness;
    function_ptr m_density;
};

template<class T>
std::ostream& operator<<(std::ostream& os, const gsMaterialMatrixBase<T>& mm)
{
    return mm.print(os);
}

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsMaterialMatrixBase
   */
  void pybind11_init_gsMaterialMatrixBase(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixBase.hpp)
#endif
