#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.hpp>

namespace gismo
{
  // Material matrix <dimension, real_t, material model, compressibility>

  // NH
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,11,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,21,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,31,false>;

  // NH_ext does not exist for incomp models
  // CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,12,false>;
  // CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,22,false>;
  // CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,32,false>;

  // MR
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,13,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,23,false>;
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,33,false>;

  // OG
  CLASS_TEMPLATE_INST gsMaterialMatrix<2,real_t,34,false>;

  #ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsMaterialMatrixNH2i(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<2,real_t>;
      using Class = gsMaterialMatrix<2,real_t,11,false>;
      py::class_<Class,Base>(m, "gsMaterialNH2i")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

    void pybind11_init_gsMaterialMatrixMR2i(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<2,real_t>;
      using Class = gsMaterialMatrix<2,real_t,13,false>;
      py::class_<Class,Base>(m, "gsMaterialMR2i")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

    void pybind11_init_gsMaterialMatrixOG2i(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<2,real_t>;
      using Class = gsMaterialMatrix<2,real_t,34,false>;
      py::class_<Class,Base>(m, "gsMaterialOG2i")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())
      ;
    }

  #endif

}

