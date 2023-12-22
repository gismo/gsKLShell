#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.hpp>

#include <gsKLShell/src/gsMaterialMatrixXml.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixLinear<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixLinear<3,real_t>;

  #ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsMaterialMatrixLinear2(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<2,real_t>;
      using Class = gsMaterialMatrixLinear<2,real_t>;
      py::class_<Class,Base>(m, "gsMaterialMatrixLinear2")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())

      // Member functions
      .def("setYoungsModulus",  &Class::setYoungsModulus,   "Sets the Young's Modulus")
      .def("setPoissonsRatio",  &Class::setPoissonsRatio,   "Sets the Poisson's Ratio")
      ;
    }

    void pybind11_init_gsMaterialMatrixLinear3(py::module &m)
    {
      using Base = gsMaterialMatrixBaseDim<3,real_t>;
      using Class = gsMaterialMatrixLinear<3,real_t>;
      py::class_<Class,Base>(m, "gsMaterialMatrixLinear3")

      // Constructors
      .def(py::init<gsFunctionSet<real_t>&,gsFunction<real_t>&>())

      // Member functions
      .def("setYoungsModulus",  &Class::setYoungsModulus,   "Sets the Young's Modulus")
      .def("setPoissonsRatio",  &Class::setPoissonsRatio,   "Sets the Poisson's Ratio")
      ;
    }

  #endif
}

