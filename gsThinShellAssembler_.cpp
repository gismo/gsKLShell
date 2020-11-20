
#include <gsCore/gsTemplateTools.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellAssembler.hpp>

// #include <gsKLShell/gsThinPlateAssembler.h>
// #include <gsKLShell/gsThinPlateAssembler.hpp>

// #include <gsKLShell/gsMembraneAssembler.h>
// #include <gsKLShell/gsMembraneAssembler.hpp>

#include <gsKLShell/gsThinShellFunctions.h>
#include <gsKLShell/gsThinShellFunctions.hpp>

#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrix.hpp>

namespace gismo
{
	CLASS_TEMPLATE_INST gsMaterialMatrix<real_t>;
	CLASS_TEMPLATE_INST gsShellStressFunction<real_t>;
	CLASS_TEMPLATE_INST gsThinShellAssembler<real_t>;
	// CLASS_TEMPLATE_INST gsThinPlateAssembler<real_t>;
	// CLASS_TEMPLATE_INST gsMembraneAssembler<real_t>;
}
