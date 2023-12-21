#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/src/gsThinShellUtils.h>

#include <gsKLShell/src/gsThinShellAssemblerDWR.h>
#include <gsKLShell/src/gsThinShellAssemblerDWR.hpp>

#include <gsKLShell/src/gsThinShellFunctions.h>
#include <gsKLShell/src/gsThinShellFunctions.hpp>

#include <gsKLShell/src/gsFunctionSum.h>


namespace gismo
{
  CLASS_TEMPLATE_INST gsFunctionSum<real_t>;

  CLASS_TEMPLATE_INST gsShellStressFunction<real_t>;

  // Shell assembler <dimension, real_t, bending terms>

  CLASS_TEMPLATE_INST gsThinShellAssemblerDWRBase<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<2,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssemblerDWR<3,real_t,true>;

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_enum_GoalFunction(py::module &m)
  {
    py::enum_<GoalFunction>(m, "GoalFunction")
        .value("Displacement" , GoalFunction::Displacement )
        .value("Stretch" , GoalFunction::Stretch )
        .value("PStrain", GoalFunction::PStrain)
        .value("PStress", GoalFunction::PStress)
        .value("MembraneStrain", GoalFunction::MembraneStrain)
        .value("MembraneStress" , GoalFunction::MembraneStress )
        .value("MembraneForce", GoalFunction::MembraneForce)
        .value("FlexuralStrain", GoalFunction::FlexuralStrain)
        .value("FlexuralStress" , GoalFunction::FlexuralStress )
        .value("FlexuralMoment", GoalFunction::FlexuralMoment)
        .value("Modal" , GoalFunction::Modal )
        .value("Buckling"   , GoalFunction::Buckling   )
        .export_values();
  }

  void pybind11_init_gsThinShellAssemblerDWRBase(py::module &m)
  {
    using Class = gsThinShellAssemblerDWRBase<real_t>;
    py::class_<Class>(m, "gsThinShellAssemblerDWRBase")

    // Member functions
    // .def("numDofs", &Class::numDofs, "Returns the number of degrees of freedom of the system")
    // .def("setSpaceBasis", &Class::setSpaceBasis, "Sets the basis used for discretization (but not for quadrature)")

    .def("assembleL", &Class::assembleL, "Assemble on the lower basis")
    .def("assembleH", &Class::assembleH, "Assemble on the higher basis")

    .def("assembleMatrixL", static_cast<ThinShellAssemblerStatus (Class::*)(                             )> (&Class::assembleMatrixL),
          "Assembles the linear matrix")
    .def("assembleMatrixL", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &)> (&Class::assembleMatrixL),
          "Assembles the nonlinear matrix")
    .def("assembleMatrixH", static_cast<ThinShellAssemblerStatus (Class::*)(                             )> (&Class::assembleMatrixH),
          "Assembles the linear matrix")
    .def("assembleMatrixH", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &)> (&Class::assembleMatrixH),
          "Assembles the nonlinear matrix")

    .def("assemblePrimalL", static_cast<ThinShellAssemblerStatus (Class::*)(                             )> (&Class::assemblePrimalL),
          "Assembles the linear primal")
    .def("assemblePrimalL", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &)> (&Class::assemblePrimalL),
          "Assembles the nonlinear primal")
    .def("assemblePrimalH", static_cast<ThinShellAssemblerStatus (Class::*)(                             )> (&Class::assemblePrimalH),
          "Assembles the linear primal")
    .def("assemblePrimalH", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &)> (&Class::assemblePrimalH),
          "Assembles the nonlinear primal")

    .def("assembleDualL", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &)> (&Class::assembleDualL),
          "Assembles the linear dual in the interior")
    .def("assembleDualH", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &)> (&Class::assembleDualH),
          "Assembles the linear dual in the interior")

    .def("assembleDualL", static_cast<ThinShellAssemblerStatus (Class::*)(const typename gsThinShellAssemblerDWRBase<real_t>::bContainer  &, const gsMultiPatch<real_t> &)> (&Class::assembleDualL),
          "Assembles the linear dual on a set of boundaries")
    .def("assembleDualH", static_cast<ThinShellAssemblerStatus (Class::*)(const typename gsThinShellAssemblerDWRBase<real_t>::bContainer  &, const gsMultiPatch<real_t> &)> (&Class::assembleDualH),
          "Assembles the linear dual on a set of boundaries")

    .def("assembleDualL", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualL),
          "Assembles the linear dual on a set of points")
    .def("assembleDualH", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualH),
          "Assembles the linear dual on a set of points")

    .def("assembleDualL", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualL),
          "Assembles the nonlinear dual in the interior")
    .def("assembleDualH", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualH),
          "Assembles the nonlinear dual in the interior")

    .def("assembleDualL", static_cast<ThinShellAssemblerStatus (Class::*)(const typename gsThinShellAssemblerDWRBase<real_t>::bContainer  &, const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualL),
          "Assembles the nonlinear dual on a set of boundaries")
    .def("assembleDualH", static_cast<ThinShellAssemblerStatus (Class::*)(const typename gsThinShellAssemblerDWRBase<real_t>::bContainer  &, const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualH),
          "Assembles the nonlinear dual on a set of boundaries")

    .def("assembleDualL", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t> &, const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualL),
          "Assembles the nonlinear dual on a set of points")
    .def("assembleDualH", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t> &, const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &)> (&Class::assembleDualH),
          "Assembles the nonlinear dual on a set of points")


    .def("assembleMassL", &Class::assembleMassL, "Assembles the mass matrix on the lower basis",
          py::arg("lumped") = false)
    .def("assembleMassH", &Class::assembleMassH, "Assembles the mass matrix on the higher basis",
          py::arg("lumped") = false)

    .def("matrixL", &Class::matrixL, "Returns the assembled matrix"         )
    .def("matrixH", &Class::matrixH, "Returns the assembled matrix"         )
    .def("massL", &Class::massL, "Returns the assembled mass matrix"         )
    .def("massH", &Class::massH, "Returns the assembled mass matrix"         )

    .def("primalL"   , &Class::primalL   , "Returns the assembled primal")
    .def("primalH"   , &Class::primalH   , "Returns the assembled primal")
    .def("dualL"   , &Class::dualL   , "Returns the assembled dual")
    .def("dualH"   , &Class::dualH   , "Returns the assembled dual")

    .def("constructSolutionL", static_cast<void                 (Class::*)(const gsMatrix<real_t> &, gsMultiPatch<real_t> &) > (&Class::constructSolutionL),
          "Constructs the solution as a gsMultiPatch")
    .def("constructSolutionL", static_cast<gsMultiPatch<real_t> (Class::*)(const gsMatrix<real_t> &                        ) > (&Class::constructSolutionL),
          "Constructs the solution as a gsMultiPatch")
    .def("constructSolutionH", static_cast<void                 (Class::*)(const gsMatrix<real_t> &, gsMultiPatch<real_t> &) > (&Class::constructSolutionH),
          "Constructs the solution as a gsMultiPatch")
    .def("constructSolutionH", static_cast<gsMultiPatch<real_t> (Class::*)(const gsMatrix<real_t> &                        ) > (&Class::constructSolutionH),
          "Constructs the solution as a gsMultiPatch")


    .def("constructDisplacementL", static_cast<gsMultiPatch<real_t> (Class::*)(const gsMatrix<real_t> &                        ) > (&Class::constructDisplacementL),
          "Constructs the displacements as a gsMultiPatch")
    .def("constructDisplacementH", static_cast<gsMultiPatch<real_t> (Class::*)(const gsMatrix<real_t> &                        ) > (&Class::constructDisplacementH),
          "Constructs the displacements as a gsMultiPatch")

    .def("computeError", static_cast<real_t (Class::*)(const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &, bool, std::string, unsigned, bool, bool)> (&Class::computeError),
          "Computes the error",
          py::arg("dualL"),
          py::arg("dualH"),
          py::arg("withLoads") = false,
          py::arg("filename") = std::string(),
          py::arg("np") = 1000,
          py::arg("parametric") = false,
          py::arg("mesh") = false
          )
    .def("computeError", static_cast<real_t (Class::*)(const gsMultiPatch<real_t> &,const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &, bool, std::string, unsigned, bool, bool)> (&Class::computeError),
          "Computes the error",
          py::arg("dualL"),
          py::arg("dualH"),
          py::arg("deformed"),
          py::arg("withLoads") = false,
          py::arg("filename") = std::string(),
          py::arg("np") = 1000,
          py::arg("parametric") = false,
          py::arg("mesh") = false
          )

    .def("computeSquaredError", static_cast<real_t (Class::*)(const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &, bool, std::string, unsigned, bool, bool)> (&Class::computeSquaredError),
          "Computes the error",
          py::arg("dualL"),
          py::arg("dualH"),
          py::arg("withLoads") = false,
          py::arg("filename") = std::string(),
          py::arg("np") = 1000,
          py::arg("parametric") = false,
          py::arg("mesh") = false
          )
    .def("computeSquaredError", static_cast<real_t (Class::*)(const gsMultiPatch<real_t> &,const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &, bool, std::string, unsigned, bool, bool)> (&Class::computeSquaredError),
          "Computes the error",
          py::arg("dualL"),
          py::arg("dualH"),
          py::arg("deformed"),
          py::arg("withLoads") = false,
          py::arg("filename") = std::string(),
          py::arg("np") = 1000,
          py::arg("parametric") = false,
          py::arg("mesh") = false
          )

    .def("computeErrorEig", static_cast<real_t (Class::*)(const real_t, const real_t, const real_t,
                                                          const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &,
                                                          const gsMultiPatch<real_t> &,
                                                          std::string, unsigned, bool, bool)> (&Class::computeErrorEig),
          "Computes the error",
          py::arg("evPrimalL"),
          py::arg("evDualL"),
          py::arg("evDualH"),
          py::arg("dualL"),
          py::arg("dualH"),
          py::arg("primal"),
          py::arg("filename") = std::string(),
          py::arg("np") = 1000,
          py::arg("parametric") = false,
          py::arg("mesh") = false
          )
    .def("computeErrorEig", static_cast<real_t (Class::*)(const real_t, const real_t, const real_t,
                                                          const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &,
                                                          const gsMultiPatch<real_t> &, const gsMultiPatch<real_t> &,
                                                          std::string, unsigned, bool, bool)> (&Class::computeErrorEig),
          "Computes the error",
          py::arg("evPrimalL"),
          py::arg("evDualL"),
          py::arg("evDualH"),
          py::arg("dualL"),
          py::arg("dualH"),
          py::arg("primal"),
          py::arg("deformed"),
          py::arg("filename") = std::string(),
          py::arg("np") = 1000,
          py::arg("parametric") = false,
          py::arg("mesh") = false
          )

    .def("computeGoal", static_cast<real_t (Class::*)(const gsMultiPatch<real_t> &)> (&Class::computeGoal), "Computes the goal functional")
    .def("computeGoal", static_cast<real_t (Class::*)(const typename gsThinShellAssemblerDWRBase<real_t>::bContainer  &, const gsMultiPatch<real_t> &)> (&Class::computeGoal), "Computes the goal functional")
    .def("computeGoal", static_cast<real_t (Class::*)(const gsMatrix<real_t> &, const gsMultiPatch<real_t> &)> (&Class::computeGoal), "Computes the goal functional")

    .def("error", &Class::error , "Returns the error")
    .def("errors", &Class::errors , "Returns the errors on the elements")
    .def("absErrors", &Class::absErrors , "Returns the absolute errors on the elements")

    .def("setGoal", &Class::setGoal , "Sets the goal functional")
    ;
  }

  void pybind11_init_gsThinShellAssemblerDWR2(py::module &m)
  {
    using Base1 = gsThinShellAssembler<2,real_t,false>;
    using Base2 = gsThinShellAssemblerDWRBase<real_t>;
    using Class = gsThinShellAssemblerDWR<2,real_t,false>;
    py::class_<Class,Base1,Base2>(m, "gsThinShellAssemblerDWR2")

      // Constructors
      .def(py::init<gsMultiPatch<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsBoundaryConditions<real_t> &,
                    gsFunction<real_t> &,
                    gsMaterialMatrixBase<real_t> * >())

      // Member functions
      ;
  }

  void pybind11_init_gsThinShellAssemblerDWR3nb(py::module &m)
  {
    using Base1 = gsThinShellAssembler<3,real_t,false>;
    using Base2 = gsThinShellAssemblerDWRBase<real_t>;
    using Class = gsThinShellAssemblerDWR<3,real_t,false>;
    py::class_<Class,Base1,Base2>(m, "gsThinShellAssemblerDWR3nb")

      // Constructors
      .def(py::init<gsMultiPatch<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsBoundaryConditions<real_t> &,
                    gsFunction<real_t> &,
                    gsMaterialMatrixBase<real_t> * >())

      // Member functions
      ;
  }

  void pybind11_init_gsThinShellAssemblerDWR3(py::module &m)
  {
    using Base1 = gsThinShellAssembler<3,real_t,true>;
    using Base2 = gsThinShellAssemblerDWRBase<real_t>;
    using Class = gsThinShellAssemblerDWR<3,real_t,true>;
    py::class_<Class,Base1,Base2>(m, "gsThinShellAssemblerDWR3")

      // Constructors
      .def(py::init<gsMultiPatch<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsBoundaryConditions<real_t> &,
                    gsFunction<real_t> &,
                    gsMaterialMatrixBase<real_t> * >())

      // Member functions
      ;
  }

#endif
}

