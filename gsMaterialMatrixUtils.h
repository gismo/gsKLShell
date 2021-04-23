/** @file gsMaterialMatrixUtils.h

    @brief Provides material matrix utilities

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
#include <gsIO/gsOptionList.h>

namespace gismo
{


/**
 * @brief      This class describes if an object is integrated through-thickness or not.
 *
 * NotIntegrated: The object has to be integrated
 * Integrated: The object is integrated
 * Constant: The object is constant through thickness, but is not integrated
 * Linear: The object is linear through thickness, but is not integrated
 *
 * @ingroup    MaterialMatrix
 */
enum class MatIntegration : short_t
{
    NotIntegrated   = 0,
    Integrated      = 1,
    Constant        = 2,
    Linear          = 3,
};


/**
 * @brief      This class describes a material model.
 *
 * @ingroup    MaterialMatrix
 *
 */
enum class Material : short_t
{
    SvK = 0,
    NH  = 1,
    NH_ext  = 2,
    MR  = 3,
    OG  = 4
};


/**
 * @brief      This class describes the way material models are implemented.
 *
 * Composite: laminate material model
 * Analytical: The expressions for Cijkl and Sij have to be provided in closed form
 * Generalized: Uses a generalized way, where only the derivatives of psi have to be implemented
 * Spectral: Implementation based on derivatives of psi w.r.t. principal stretches
 *
 * @ingroup    MaterialMatrix
 */
enum class Implementation : short_t
{
    Composite   = 0,
    Analytical  = 1,
    Generalized = 2,
    Spectral    = 3
};


/**
 * @brief      This class describes the output type.
 *
 * Generic: unspecified output type
 * Density
 * VectorN: Membrane forces
 * VectorM: Bending moments
 * MatrixA: 0th order moment of the differential of the material matrix w.r.t. membrane strains
 * MatrixB: 1st order moment of the differential of the material matrix w.r.t. membrane strains
 * MatrixC: 0th order moment of the differential of the material matrix w.r.t. bending strains
 * MatrixD: 1st order moment of the differential of the material matrix w.r.t. bending strains
 * PStressN: membrane principal stress
 * PStressN: bending principal stress
 * Stretch: Principal stretch
 * StretchDir: Principal stretch directions
 *
 * @ingroup    MaterialMatrix
 *
 */
enum class MaterialOutput : short_t
{
    Generic = 0,
    Density = 1,
    VectorN = 2,
    VectorM = 3,
    MatrixA  = 4,
    MatrixB  = 5,
    MatrixC  = 6,
    MatrixD  = 7,
    PStressN = 8,
    PStressM = 9,
    Stretch = 10,       // ONLY ON MID-PLANE
    StretchDir = 11,    // ONLY ON MID-PLANE
    Transformation = 12,
};


/**
 * @brief      Encodes the material model and implementation
 *
 * @tparam     material        Material model
 * @tparam     implementation  The way it is implemented
 *
 * @ingroup    MaterialMatrix
 *
 */
template<enum Material material, enum Implementation implementation>
struct encodeMat_id {
  static const constexpr short_t id = 10*(short_t)implementation + (short_t)material;
};

/**
 * @brief      Decodes the material model and implementation
 *
 * @tparam     id    identifier (from encoder)
 *
 * @ingroup    MaterialMatrix
 *
 */
template<short_t id>
struct decodeMat_id {
  static const constexpr enum Material          material = (enum Material)(id%10);
  static const constexpr enum Implementation    implementation = (enum Implementation)((id/10)%10);
};

} // namespace

