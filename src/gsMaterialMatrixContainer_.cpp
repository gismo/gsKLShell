#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixContainer.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixContainer<real_t>;

namespace internal
{
    CLASS_TEMPLATE_INST gsXml< gsMaterialMatrixContainer<real_t>>;
}
}

