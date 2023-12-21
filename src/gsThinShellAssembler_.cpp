#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/src/gsThinShellUtils.h>

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsThinShellAssembler.hpp>

#include <gsKLShell/src/gsThinShellFunctions.h>
#include <gsKLShell/src/gsThinShellFunctions.hpp>

#include <gsKLShell/src/gsFunctionSum.h>


namespace gismo
{
  CLASS_TEMPLATE_INST gsFunctionSum<real_t>;

  CLASS_TEMPLATE_INST gsShellStressFunction<real_t>;

  // Shell assembler <dimension, real_t, bending terms>

  CLASS_TEMPLATE_INST gsThinShellAssemblerBase<real_t>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<2,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,false>;
  CLASS_TEMPLATE_INST gsThinShellAssembler<3,real_t,true>;

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_enum_gsThinShellAssemblerStatus(py::module &m)
  {
      py::enum_<ThinShellAssemblerStatus>(m, "assemblerStatus")
          .value("Success"   , ThinShellAssemblerStatus::Success )
          .value("AssemblyError"     , ThinShellAssemblerStatus::AssemblyError )
          .value("DimensionError", ThinShellAssemblerStatus::DimensionError)
          .export_values();
  }

  void pybind11_init_gsThinShellAssemblerBase(py::module &m)
  {
    using Class = gsThinShellAssemblerBase<real_t>;
    py::class_<Class>(m, "gsThinShellAssemblerBase")

    // Member functions
    .def("numDofs", &Class::numDofs, "Returns the number of degrees of freedom of the system")
    .def("setSpaceBasis", &Class::setSpaceBasis, "Sets the basis used for discretization (but not for quadrature)")

    .def("assemble", static_cast<ThinShellAssemblerStatus (Class::*)()> (&Class::assemble),
          "Assembles the linear system")
    .def("assemble", static_cast<ThinShellAssemblerStatus (Class::*)(const gsFunctionSet<real_t> &, const bool, const bool)> (&Class::assemble),
          "Assembles the nonlinear system (matrix optional)",
          py::arg("deformed"),
          py::arg("matrix") = false,
          py::arg("homogenize") = true)
    .def("assemble", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t> &     , const bool, const bool)> (&Class::assemble),
          "Assembles the nonlinear system (matrix optional)",
          py::arg("solVector"),
          py::arg("matrix") = false,
          py::arg("homogenize") = true)

    .def("assembleMatrix", static_cast<ThinShellAssemblerStatus (Class::*)(const gsFunctionSet<real_t> &)> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix")
    .def("assembleMatrix", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t> &     )> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix")

    .def("assembleMatrix", static_cast<ThinShellAssemblerStatus (Class::*)(const gsFunctionSet<real_t> &, const gsFunctionSet<real_t> &, gsMatrix<real_t> & update)> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix using the Mixed Integration Point (MIP) method")
    .def("assembleMatrix", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t>      &, const gsMatrix<real_t>      &                           )> (&Class::assembleMatrix),
          "Assembles the nonlinear matrix using the Mixed Integration Point (MIP) method")

    .def("assembleVector", static_cast<ThinShellAssemblerStatus (Class::*)(const gsFunctionSet<real_t> &, const bool )> (&Class::assembleVector),
          "Assembles the nonlinear vector",
          py::arg("deformed"),
          py::arg("homogenize") = true)
    .def("assembleVector", static_cast<ThinShellAssemblerStatus (Class::*)(const gsMatrix<real_t> &     , const bool)> (&Class::assembleVector),
          "Assembles the nonlinear vector",
          py::arg("solVector"),
          py::arg("homogenize") = true)
    
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
    .def("constructSolutionVector", &Class::constructSolutionVector, "Returns the solution vector based on a multipatch solution")

    .def("getArea"   , &Class::getArea   , "Returns the area of a geometry")


    .def("setGeometry"   , &Class::setGeometry   , "Sets the geometry")
    .def("setBasis"   , &Class::setBasis   , "Sets the basis")
    .def("setSpaceBasis"   , &Class::setSpaceBasis   , "Sets the basis for the space")

    .def("setPointLoads"      , &Class::setPointLoads       , "Sets point loads")
    .def("updateBCs"          , &Class::updateBCs           , "update BCs")
    .def("homogenizeDirichlet", &Class::homogenizeDirichlet , "Homogenize dirichlet BCs")

    .def("options"   , &Class::options   , "Access options")
    .def("addWeakC0"   , &Class::addWeakC0   , "Adds interfaces for weak C0 coupling")
    .def("addWeakC1"   , &Class::addWeakC1   , "Adds interfaces for weak C1 coupling")

    .def("geometry"   , static_cast<const gsMultiPatch<real_t> & (Class::*)() const> (&Class::geometry), "Gets the geometry")
    .def("geometry"   , static_cast<      gsMultiPatch<real_t> & (Class::*)()      > (&Class::geometry), "Gets the geometry")
    .def("basis"      , static_cast<const gsMultiBasis<real_t> & (Class::*)() const> (&Class::basis)   , "Gets the basis")
    .def("basis"      , static_cast<      gsMultiBasis<real_t> & (Class::*)()      > (&Class::basis)   , "Gets the basis")
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

