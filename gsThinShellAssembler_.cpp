#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellAssembler.hpp>

#include <gsKLShell/gsThinShellFunctions.h>
#include <gsKLShell/gsThinShellFunctions.hpp>

#include <gsKLShell/gsFunctionSum.h>


namespace gismo
{
  CLASS_TEMPLATE_INST gsFunctionSum<real_t>;

  CLASS_TEMPLATE_INST gsShellStressFunction<real_t>;

  // Shell assembler <dimension, real_t, bending terms>

  CLASS_TEMPLATE_INST gsThinShellAssemblerBase<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<2,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,true>;

#ifdef GISMO_BUILD_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsThinShellAssemblerBase(py::module &m)
  {
    using Class = gsThinShellAssemblerBase<real_t>;
    py::class_<Class>(m, "gsThinShellAssemblerBase")

    // Member functions
    .def("numDofs", &Class::numDofs, "Returns the number of degrees of freedom of the system")
    .def("setSpaceBasis", &Class::setSpaceBasis, "Sets the basis used for discretization (but not for quadrature)")

    .def("assemble", static_cast<void (Class::*)()> (&Class::assemble),
          "Assembles the linear system")
    .def("assemble", static_cast<void (Class::*)(const gsFunctionSet<real_t> &, bool)> (&Class::assemble),
          "Assembles the nonlinear system (matrix optional)")//,
          // py::arg("Matrix") = false)
    .def("assemble", static_cast<void (Class::*)(const gsMatrix<real_t> &     , bool)> (&Class::assemble),
          "Assembles the nonlinear system (matrix optional)")//,
          // py::arg("Matrix") = false)

    .def("assembleMatrix", static_cast<void (Class::*)(const gsFunctionSet<real_t> &)> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix")
    .def("assembleMatrix", static_cast<void (Class::*)(const gsMatrix<real_t> &     )> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix")

    .def("assembleMatrix", static_cast<void (Class::*)(const gsFunctionSet<real_t> &, const gsFunctionSet<real_t> &, gsMatrix<real_t> & update)> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix using the Mixed Integration Point (MIP) method")
    .def("assembleMatrix", static_cast<void (Class::*)(const gsMatrix<real_t>      &, const gsMatrix<real_t>      &                           )> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix using the Mixed Integration Point (MIP) method")

    .def("assembleVector", static_cast<void (Class::*)(const gsFunctionSet<real_t> &)> (&Class::assembleVector),
          "Assembles the nonlinear vector")
    .def("assembleVector", static_cast<void (Class::*)(const gsMatrix<real_t> &     )> (&Class::assembleVector),
          "Assembles the nonlinear vector")

    .def("assembleMass", &Class::assembleMass, "Assembles the mass matrix",
          py::arg("lumped") = false)

    .def("matrix", &Class::matrix, "Returns the assembled matrix"         )
    .def("rhs"   , &Class::rhs   , "Returns the assembled right-hand side")

    .def("constructSolution", static_cast<void                 (Class::*)(const gsMatrix<real_t> &, gsMultiPatch<real_t> &) const> (&Class::constructSolution),
          "Constructs the solution as a gsMultiPatch")
    .def("constructSolution", static_cast<gsMultiPatch<real_t> (Class::*)(const gsMatrix<real_t> &                        ) const> (&Class::constructSolution),
          "Constructs the solution as a gsMultiPatch")

    .def("constructDisplacement", static_cast<void                 (Class::*)(const gsMatrix<real_t> &, gsMultiPatch<real_t> &) const> (&Class::constructDisplacement),
          "Constructs the displacements as a gsMultiPatch")
    .def("constructDisplacement", static_cast<gsMultiPatch<real_t> (Class::*)(const gsMatrix<real_t> &                        ) const> (&Class::constructDisplacement),
          "Constructs the displacements as a gsMultiPatch")
    ;
  }

  void pybind11_init_gsThinShellAssembler2(py::module &m)
  {
    using Base = gsThinShellAssemblerBase<real_t>;
    using Class = gsThinShellAssembler<2,real_t,false>;
    py::class_<Class,Base>(m, "gsThinShellAssembler2")

      // Constructors
      .def(py::init<gsMultiPatch<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsBoundaryConditions<real_t> &,
                    gsFunction<real_t> &,
                    gsMaterialMatrixBase<real_t> * >())

      // Member functions
      ;
  }

  void pybind11_init_gsThinShellAssembler3nb(py::module &m)
  {
    using Base = gsThinShellAssemblerBase<real_t>;
    using Class = gsThinShellAssembler<3,real_t,false>;
    py::class_<Class,Base>(m, "gsThinShellAssembler3nb")

      // Constructors
      .def(py::init<gsMultiPatch<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsBoundaryConditions<real_t> &,
                    gsFunction<real_t> &,
                    gsMaterialMatrixBase<real_t> * >())

      // Member functions
      ;
  }

  void pybind11_init_gsThinShellAssembler3(py::module &m)
  {
    using Base = gsThinShellAssemblerBase<real_t>;
    using Class = gsThinShellAssembler<3,real_t,true>;
    py::class_<Class,Base>(m, "gsThinShellAssembler3")

      // Constructors
      .def(py::init<gsMultiPatch<real_t>&,
                    gsMultiBasis<real_t>&,
                    gsBoundaryConditions<real_t> &,
                    gsFunction<real_t> &,
                    gsMaterialMatrixBase<real_t> * >())

      // Member functions
      ;
  }

#endif
}

