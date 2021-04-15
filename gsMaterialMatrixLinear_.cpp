#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixLinear.hpp>


namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixLinear<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixLinear<3,real_t>;
}

