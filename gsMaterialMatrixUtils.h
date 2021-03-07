/** @file gsMaterialMatrix.h

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

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

enum class MatIntegration : short_t
{
    NotIntegrated   = 0,
    Integrated      = 1,
    Constant        = 2,
    Linear          = 3,
};

enum class Material : short_t
{
    SvK = 0,
    NH  = 1,
    NH_ext  = 2,
    MR  = 3,
    OG  = 4
};

enum class Implementation : short_t
{
    Composite   = 0,
    Analytical  = 1,
    Generalized = 2,
    Spectral    = 3
};

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
    Stretch = 10,
    StretchDir = 11,
};

template<enum Material material, enum Implementation implementation>
struct encodeMat_id {
  static const constexpr short_t id = 10*(short_t)implementation + (short_t)material;
};

template<short_t id>
struct decodeMat_id {
  static const constexpr enum Material          material = (enum Material)(id%10);
  static const constexpr enum Implementation    implementation = (enum Implementation)((id/10)%10);
};

} // namespace

