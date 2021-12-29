#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixTFT.h>
#include <gsKLShell/gsMaterialMatrixTFT.hpp>


namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixTFT<3,real_t>;
}

