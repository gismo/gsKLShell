
#include <gsCore/gsTemplateTools.h>

// #include <gsThinShell2/gsThinShellFunctions.h>
// #include <gsThinShell2/gsThinShellFunctions.hpp>

#include <gsThinShell2/gsThinShellAssembler.h>
#include <gsThinShell2/gsThinShellAssembler.hpp>

#include <gsThinShell2/gsMaterialMatrix.h>
#include <gsThinShell2/gsMaterialMatrix.hpp>

namespace gismo
{
    	CLASS_TEMPLATE_INST gsThinShellAssembler<real_t>;
    	CLASS_TEMPLATE_INST gsMaterialMatrix<real_t>;
}
