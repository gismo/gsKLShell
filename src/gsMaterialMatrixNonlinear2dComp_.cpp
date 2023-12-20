#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.hpp>

namespace gismo
{
  // Material matrix <dimension, real_t, material model, compressibility>

  // NH
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,11,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,21,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,31,true>;

  // NH_ext
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,12,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,22,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,32,true>;

  // MR
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,13,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,23,true>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,33,true>;

  // OG
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,34,true>;

  #ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsMaterialMatrixNH2c(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<2,real_t>;
      using Class = gsMaterialMatrix<2,real_t,11,true>;
      py::class_<Class,Base>(m, "gsMaterialNH2c")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

    void pybind11_init_gsMaterialMatrixMR2c(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<2,real_t>;
      using Class = gsMaterialMatrix<2,real_t,13,true>;
      py::class_<Class,Base>(m, "gsMaterialMR2c")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

    void pybind11_init_gsMaterialMatrixOG2c(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<2,real_t>;
      using Class = gsMaterialMatrix<2,real_t,34,true>;
      py::class_<Class,Base>(m, "gsMaterialOG2c")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

  #endif
}

