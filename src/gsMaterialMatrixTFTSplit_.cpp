/** @file gsMaterialMatrixTFTSplit_.cpp

    @brief Explicit instantiation of gsMaterialMatrixTFTSplit

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
*/

#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixTFTSplit.h>
#include <gsKLShell/src/gsMaterialMatrixTFTSplit.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixTFTSplit<2,real_t,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFTSplit<2,real_t,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFTSplit<3,real_t,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFTSplit<3,real_t,false>;
}
