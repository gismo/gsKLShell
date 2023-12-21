#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixBase.hpp>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMaterialMatrixBase<real_t>;

  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrixBaseDim<3,real_t>;

  // namespace internal
  // {
  //   CLASS_TEMPLATE_INST gsXml<gsMaterialMatrixBase<real_t>>;
  // }

  #ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsMaterialMatrixBase(py::module &m)
  {
    using Class = gsMaterialMatrixBase<real_t>;
    py::class_<Class>(m, "gsMaterialMatrixBase")

    // Member functions
    .def("setOptions", &Class::setOptions, "Sets the options")
    .def("density_into", &Class::density_into, "Computes the density into a matrix")
    .def("pstretch_into", &Class::pstretch_into, "Computes the stretches into a matrix")
    .def("pstretchDir_into", &Class::pstretchDir_into, "Computes the stretch directions into a matrix")
    .def("thickness_into", &Class::thickness_into, "Computes the thickness into a matrix")
    .def("transform_into", &Class::transform_into, "Computes the stretch transformation into a matrix")
    .def("eval3D_matrix", static_cast<gsMatrix<real_t> (Class::*)(   const index_t, const gsMatrix<real_t> &, const gsMatrix<real_t> &, enum MaterialOutput    ) const> (&Class::eval3D_matrix), "Evaluates the material matrix")
    .def("eval3D_matrix", static_cast<gsMatrix<real_t> (Class::*)(   const index_t, const gsVector<real_t> &, const          real_t  &, enum MaterialOutput    ) const> (&Class::eval3D_matrix), "Evaluates the material matrix")
    .def("eval3D_vector", static_cast<gsMatrix<real_t> (Class::*)(   const index_t, const gsMatrix<real_t> &, const gsMatrix<real_t> &, enum MaterialOutput    ) const> (&Class::eval3D_vector), "Evaluates the material matrix")
    .def("eval3D_vector", static_cast<gsMatrix<real_t> (Class::*)(   const index_t, const gsVector<real_t> &, const          real_t  &, enum MaterialOutput    ) const> (&Class::eval3D_vector), "Evaluates the material matrix")
    .def("eval3D_pstress", &Class::eval3D_pstress, "Evaluates the principal stress vector")

    .def("setYoungsModulus",  &Class::setYoungsModulus,   "Sets the Young's Modulus")
    .def("setPoissonsRatio",  &Class::setPoissonsRatio,   "Sets the Poisson's Ratio")
    .def("setRatio"        ,  &Class::setRatio        ,   "Sets the Ratio for MR model")
    .def("setMu"           ,  &Class::setMu           ,   "Sets the Mu_i for OG model")
    .def("setAlpha"        ,  &Class::setAlpha        ,   "Sets the Alpha_i for OG model")
    .def("setDensity"      ,  static_cast<void (Class::*)(   const gsFunctionSet<real_t> &    ) > (&Class::setDensity)      ,   "Sets the Density")
    .def("setThickness"    ,  static_cast<void (Class::*)(   const gsFunctionSet<real_t> &    ) > (&Class::setThickness)    ,   "Sets the Thickness")
    ;
  }

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

