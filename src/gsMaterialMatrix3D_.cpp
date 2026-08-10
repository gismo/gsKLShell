#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrix3D.h>

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsKLShell/src/gsMaterialMatrix3D.hpp>

namespace gismo
{
  // --------------------------------------------------------------------------
  // Cross-output cache diagnostics counters (gsKLShell issue #28, constitutive
  // side). Defined here so each is compiled exactly once into libgismo, giving a
  // single exported symbol shared across the shared-library boundary. Declared
  // GISMO_EXPORT in gsMaterialMatrix3D.h (see there for the -fvisibility=hidden
  // rationale). Not thread-exact; single-threaded benchmarking only.
  // --------------------------------------------------------------------------
  namespace {
      size_t & _gsMaterialMatrix3DSweepCounter()
      {
          static size_t counter = 0;
          return counter;
      }
      size_t & _gsMaterialMatrix3DHitCounter()
      {
          static size_t counter = 0;
          return counter;
      }
  }
  size_t gsMaterialMatrix3DSweeps()          { return _gsMaterialMatrix3DSweepCounter(); }
  void   gsMaterialMatrix3DResetSweeps()     { _gsMaterialMatrix3DSweepCounter() = 0; }
  void   gsMaterialMatrix3DIncrementSweeps() { ++_gsMaterialMatrix3DSweepCounter(); }

  size_t gsMaterialMatrix3DHits()            { return _gsMaterialMatrix3DHitCounter(); }
  void   gsMaterialMatrix3DResetHits()       { _gsMaterialMatrix3DHitCounter() = 0; }
  void   gsMaterialMatrix3DIncrementHits()   { ++_gsMaterialMatrix3DHitCounter(); }

  CLASS_TEMPLATE_INST gsMaterialMatrix3D<2,real_t>;
  CLASS_TEMPLATE_INST gsMaterialMatrix3D<3,real_t>;
}

#endif // gsPhaseFieldFracture_ENABLED
