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
 * @ingroup    MaterialMatrix
 */
template <class T>
class gsMaterialMatrixBase
{
public:

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
    inline virtual gsOptionList & options() = 0;
    /**
     * @brief      Sets the options
     *
     * @param[in]  opt   @ref gsOptionList
     */
    inline virtual void setOptions(gsOptionList opt) = 0;

    /**
     * @brief      Evaluates the density multiplied by the thickness of the shell (scalar)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The resut
     */
    inline virtual void    density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;
    /**
     * @brief      Evaluates the stretches in the shell (3x1 vector)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The resut
     */
    inline virtual void    stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;
    /**
     * @brief      Evaluates the directions of the stretches in the shell (3x1 vector per direction)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The resut
     */
    inline virtual void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;
    /**
     * @brief      Evaluates the thickness of the shell (scalar)
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The resut
     */
    inline virtual void  thickness_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;

    /**
     * @brief      ???????????????
     *
     * @param[in]  patch   The patch to be evaluated on
     * @param[in]  u       The in-plane shell coordinates to be eveluated on
     * @param      result  The resut
     */
    inline virtual void  transform_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;

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
    inline virtual gsMatrix<T>  eval3D_matrix(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const = 0;
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
    inline virtual gsMatrix<T>  eval3D_vector(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const = 0;
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
    inline virtual gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T>& u, const gsMatrix<T>& z, enum MaterialOutput out) const = 0;

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    inline virtual void setParameters(const std::vector<gsFunction<T>*> &pars) =0;

    /**
     * @brief      Prints info
     */
    inline virtual void info() const = 0;

    void setDeformed(const gsFunctionSet<T> & deformed) {m_defpatches = &deformed; }

protected:
    const gsFunctionSet<T> * m_defpatches;
};

} // namespace
