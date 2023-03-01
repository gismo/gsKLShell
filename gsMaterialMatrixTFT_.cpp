#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixTFT.h>
#include <gsKLShell/gsMaterialMatrixTFT.hpp>


namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<2,real_t,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<2,real_t,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<3,real_t,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<3,real_t,false>;
}

