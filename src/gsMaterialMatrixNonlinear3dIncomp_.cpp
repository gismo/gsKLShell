#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.hpp>

namespace gismo
{
  // Material matrix <dimension, real_t, material model, compressibility>

  // NH
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,11,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,21,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,31,false>;

  // NH_ext does not exist for incomp models
  // CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,12,false>;
  // CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,22,false>;
  // CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,32,false>;

  // MR
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,13,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,23,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,33,false>;

  // OG
  CLASS_TEMPLATE_INST gsMaterialMatrix<3,real_t,34,false>;

  #ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsMaterialMatrixNH3i(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<3,real_t>;
      using Class = gsMaterialMatrix<3,real_t,11,false>;
      py::class_<Class,Base>(m, "gsMaterialNH3i")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

    void pybind11_init_gsMaterialMatrixMR3i(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<3,real_t>;
      using Class = gsMaterialMatrix<3,real_t,13,false>;
      py::class_<Class,Base>(m, "gsMaterialMR3i")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

    void pybind11_init_gsMaterialMatrixOG3i(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<3,real_t>;
      using Class = gsMaterialMatrix<3,real_t,34,false>;
      py::class_<Class,Base>(m, "gsMaterialOG3i")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

  #endif

}

