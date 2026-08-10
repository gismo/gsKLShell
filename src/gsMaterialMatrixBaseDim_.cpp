#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.hpp>

namespace gismo
{
  // ----------------------------------------------------------------------------
  // Diagnostics counter (gsKLShell issue #28). Defined here so it is compiled
  // exactly once into libgismo, giving a single exported symbol shared across
  // the shared-library boundary. Declared GISMO_EXPORT in gsMaterialMatrixBaseDim.h
  // (see there for the -fvisibility=hidden rationale). Not thread-exact; single-
  // threaded benchmarking only.
  // ----------------------------------------------------------------------------
  namespace {
      size_t & _gsMaterialMatrixComputePointsCounter()
      {
          static size_t counter = 0;
          return counter;
      }
      // Misses = full executions that passed the same-input guard (issue #28).
      size_t & _gsMaterialMatrixComputePointsMissCounter()
      {
          static size_t counter = 0;
          return counter;
      }
  }
  size_t gsMaterialMatrixComputePointsCalls()
  { return _gsMaterialMatrixComputePointsCounter(); }
  void gsMaterialMatrixResetComputePointsCalls()
  { _gsMaterialMatrixComputePointsCounter() = 0; }
  void gsMaterialMatrixIncrementComputePointsCalls()
  { ++_gsMaterialMatrixComputePointsCounter(); }

  size_t gsMaterialMatrixComputePointsMisses()
  { return _gsMaterialMatrixComputePointsMissCounter(); }
  void gsMaterialMatrixResetComputePointsMisses()
  { _gsMaterialMatrixComputePointsMissCounter() = 0; }
  void gsMaterialMatrixIncrementComputePointsMisses()
  { ++_gsMaterialMatrixComputePointsMissCounter(); }

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

