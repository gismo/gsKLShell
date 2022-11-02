#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsThinShellAssemblerDWR.hpp>

namespace gismo
{
  // CLASS_TEMPLATE_INST gsThinShellDWRFunction<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWRBase<real_t>;
  // CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<2,real_t,false>;
  // CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true>;
}

