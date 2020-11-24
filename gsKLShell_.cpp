
#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellAssembler.hpp>

// #include <gsKLShell/gsThinPlateAssembler.h>
// #include <gsKLShell/gsThinPlateAssembler.hpp>

// #include <gsKLShell/gsMembraneAssembler.h>
// #include <gsKLShell/gsMembraneAssembler.hpp>

#include <gsKLShell/gsThinShellFunctions.h>
#include <gsKLShell/gsThinShellFunctions.hpp>

#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrix.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsShellStressFunction<real_t>;

  // Shell assembler <dimension, real_t, bending terms>

  CLASS_TEMPLATE_INST gsThinShellAssemblerBase<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<2,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,true>;

  // Material matrix <dimension, real_t, material model, compressibility>

  CLASS_TEMPLATE_INST gsMaterialMatrixBase<real_t>;

  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,0,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,2,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,12,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,22,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,3,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,13,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,23,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,14,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,5,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,15,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,25,true>;

  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,0,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,2,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,12,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,22,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,3,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,13,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,23,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,14,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,5,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,15,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,25,false>;

  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,0,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,2,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,12,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,22,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,3,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,13,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,23,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,14,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,5,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,15,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,25,true>;

  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,0,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,2,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,12,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,22,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,3,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,13,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,23,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,14,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,5,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,15,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,25,false>;

}

