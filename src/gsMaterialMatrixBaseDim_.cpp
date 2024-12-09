#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<3,real_t>;

}

