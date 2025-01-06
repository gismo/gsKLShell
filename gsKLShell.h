#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixComposite.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsKLShell/src/gsMaterialMatrixXml.h>

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsThinShellAssemblerDWR.h>

#include <gsKLShell/src/gsThinShellFunctions.h>
#include <gsKLShell/src/gsThinShellUtils.h>

#include <gsKLShell/src/getMaterialMatrix.h>

namespace gismo
{
#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsKLShell(py::module &m);

#endif
}

