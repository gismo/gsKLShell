#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellAssembler.hpp>

#include <gsKLShell/gsThinShellFunctions.h>
#include <gsKLShell/gsThinShellFunctions.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsShellStressFunction<real_t>;

  // Shell assembler <dimension, real_t, bending terms>

  CLASS_TEMPLATE_INST gsThinShellAssemblerBase<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<2,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,true>;
}

