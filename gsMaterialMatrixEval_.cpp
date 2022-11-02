#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixBase.h>

#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixEval.hpp>

namespace gismo
{
  // Material matrix <dimension, real_t, material model, compressibility>
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Generic>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Density>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::VectorN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::VectorM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Stretch>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStress>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStressDir>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::StretchTransform>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStressTransform>;
}

