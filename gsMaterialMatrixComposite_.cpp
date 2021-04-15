#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixComposite.h>
#include <gsKLShell/gsMaterialMatrixComposite.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixComposite<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixComposite<3,real_t>;
}

