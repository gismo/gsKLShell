
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
  CLASS_TEMPLATE_INST gsMaterialMatrix<real_t>;
  CLASS_TEMPLATE_INST gsShellStressFunction<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<2,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,true>;
  // CLASS_TEMPLATE_INST gsThinPlateAssembler<real_t>;
  // CLASS_TEMPLATE_INST gsMembraneAssembler<real_t>;
}
