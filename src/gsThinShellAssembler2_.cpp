#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/src/gsThinShellUtils.h>

// Included OUTSIDE the guard on purpose: it is the include that (transitively)
// pulls in the generated gsCore/gsConfigExt.h, which is where
// gsPhaseFieldFracture_ENABLED is defined. Same pattern as
// gsMaterialMatrix3D_.cpp:1-7.
#include <gsKLShell/src/gsThinShellAssembler2.h>

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsKLShell/src/gsThinShellAssembler2.hpp>

namespace gismo
{
  // Shell assembler <dimension, real_t, bending terms>
  // Mirrors gsThinShellAssembler_.cpp:20-22.
  CLASS_TEMPLATE_INST gsThinShellAssembler2<2,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler2<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler2<3,real_t,true>;
}

#endif // gsPhaseFieldFracture_ENABLED
