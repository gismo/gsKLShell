#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>

#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsThinShellAssemblerDWR.h>


namespace gismo
{

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsKLShell(py::module &m)
  {
    gismo::pybind11_init_gsMaterialMatrixBase( m );

    gismo::pybind11_init_gsMaterialMatrixBaseDim2( m );
    gismo::pybind11_init_gsMaterialMatrixBaseDim3( m );

    gismo::pybind11_init_gsMaterialMatrixLinear2( m );
    gismo::pybind11_init_gsMaterialMatrixLinear3( m );

    gismo::pybind11_init_gsMaterialMatrixNH2i( m );
    gismo::pybind11_init_gsMaterialMatrixNH2c( m );

    gismo::pybind11_init_gsMaterialMatrixNH3i( m );
    gismo::pybind11_init_gsMaterialMatrixNH3c( m );

    gismo::pybind11_init_gsMaterialMatrixMR2i( m );
    gismo::pybind11_init_gsMaterialMatrixMR2c( m );

    gismo::pybind11_init_gsMaterialMatrixMR3i( m );
    gismo::pybind11_init_gsMaterialMatrixMR3c( m );

    gismo::pybind11_init_gsMaterialMatrixOG2i( m );
    gismo::pybind11_init_gsMaterialMatrixOG2c( m );

    gismo::pybind11_init_gsMaterialMatrixOG3i( m );
    gismo::pybind11_init_gsMaterialMatrixOG3c( m );

    gismo::pybind11_init_gsThinShellAssemblerBase( m );
    gismo::pybind11_enum_gsThinShellAssemblerStatus( m );

    gismo::pybind11_init_gsThinShellAssembler2( m );
    gismo::pybind11_init_gsThinShellAssembler3( m );
    gismo::pybind11_init_gsThinShellAssembler3nb( m );

    gismo::pybind11_enum_GoalFunction( m );

    gismo::pybind11_init_gsThinShellAssemblerDWRBase( m );

    gismo::pybind11_init_gsThinShellAssemblerDWR2( m );
    gismo::pybind11_init_gsThinShellAssemblerDWR3( m );
    gismo::pybind11_init_gsThinShellAssemblerDWR3nb( m );
  }

#endif
}

