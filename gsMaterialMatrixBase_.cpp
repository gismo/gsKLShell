#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixBaseDim.h>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixBase<real_t>;

  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<3,real_t>;
}

