#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixBase.h>

#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrix.hpp>

namespace gismo
{
  // Material matrix <dimension, real_t, material model, compressibility>
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,11,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,21,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,31,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,12,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,22,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,32,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,13,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,23,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,33,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,34,true>;
}

