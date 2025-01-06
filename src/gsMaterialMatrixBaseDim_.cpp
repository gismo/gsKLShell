#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<3,real_t>;

  #ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsMaterialMatrixBaseDim2(py::module &m)
  {
    using Base  = gsMaterialMatrixBase<real_t>;
    using Class = gsMaterialMatrixBaseDim<2,real_t>;
    py::class_<Class,Base>(m, "gsMaterialMatrixBaseDim2")
    ;
  }

  void pybind11_init_gsMaterialMatrixBaseDim3(py::module &m)
  {
    using Base  = gsMaterialMatrixBase<real_t>;
    using Class = gsMaterialMatrixBaseDim<3,real_t>;
    py::class_<Class,Base>(m, "gsMaterialMatrixBaseDim3")
    ;
  }

  #endif

}

