#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsThinShellAssemblerDWR.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWRBase<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<2,real_t,false,GoalFunction::DisplacementNorm>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<2,real_t,false,GoalFunction::Displacement>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<2,real_t,false,GoalFunction::MembraneStrain>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<2,real_t,false,GoalFunction::MembraneStress>;

  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,false,GoalFunction::DisplacementNorm>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,false,GoalFunction::Displacement>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,false,GoalFunction::MembraneStrain>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,false,GoalFunction::MembraneStress>;

  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::DisplacementNorm>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Displacement>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStrain>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStress>;
}

