#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsThinShellAssemblerDWR.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWRBase<real_t>;

  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Displacement,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Stretch,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStrain,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStrain,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStress,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStress,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneForce,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Modal,9>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Buckling,9>;

  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Displacement,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Stretch,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStrain,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStrain,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStress,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStress,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneForce,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Modal,0>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Buckling,0>;

  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Displacement,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Stretch,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStrain,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStrain,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStress,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStress,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneForce,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Modal,1>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Buckling,1>;

  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Displacement,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Stretch,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStrain,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStrain,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStress,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStress,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneForce,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Modal,2>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Buckling,2>;
}

