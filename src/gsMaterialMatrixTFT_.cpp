#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixTFT.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.hpp>


namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<2,real_t,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<2,real_t,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<3,real_t,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<3,real_t,false>;
}

