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
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStressN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStressM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStrainN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStrainM>;
  // CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStressG>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Stretch>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Transformation>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::CovTransform>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::ConTransform>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::TensionField>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Strain>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::StrainN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::StrainM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Thickness>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Parameters>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Deformation>;

    // Material matrix <dimension, real_t, material model, compressibility>
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Generic>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Density>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::VectorN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::VectorM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::MatrixA>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::MatrixB>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::MatrixC>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::MatrixD>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::PStressN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::PStressM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::PStrainN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::PStrainM>;
  // CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::PStressG>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Stretch>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::StretchDir>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Transformation>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::CovTransform>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::ConTransform>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::TensionField>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Strain>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::StrainN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::StrainM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Thickness>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Parameters>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEvalSingle<real_t,MaterialOutput::Deformation>;
}

