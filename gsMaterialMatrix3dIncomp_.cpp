#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixBase.h>

#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrix.hpp>

namespace gismo
{
  // Material matrix <dimension, real_t, material model, compressibility>
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,11,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,21,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,31,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,12,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,22,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,32,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,13,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,23,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,33,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,34,false>;

}

