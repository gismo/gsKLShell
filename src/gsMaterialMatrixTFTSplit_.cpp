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
