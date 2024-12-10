/** @file gsMaterialMatrixXml.h

    @brief TO DO

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): TO DO
*/

#pragma once

#include <gsIO/gsXml.h>

namespace gismo
{
namespace internal
{

template<class Object>
Object getMaterialMatrixFromXml ( gsXmlNode * node);

template<class Object>
gsXmlNode * putMaterialMatrixToXml ( Object const & obj, gsXmlTree & data);

}
}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixXml.hpp)
#endif
