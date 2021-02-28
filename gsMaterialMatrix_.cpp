#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/gsMaterialMatrixBase.h>

#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixEval.hpp>

#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrix.hpp>

namespace gismo
{
  // Material matrix <dimension, real_t, material model, compressibility>

  CLASS_TEMPLATE_INST gsMaterialMatrixBase<real_t>;

  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Density>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::VectorN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::VectorM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStressN>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::PStressM>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::Stretch>;
  CLASS_TEMPLATE_INST gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir>;

  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,10,true>;
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

  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,10,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,11,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,21,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,31,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,12,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,22,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,32,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,13,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,23,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,33,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,34,false>;

  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,10,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,11,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,21,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,31,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,12,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,22,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,32,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,13,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,23,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,33,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,34,true>;

  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,10,false>;
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

