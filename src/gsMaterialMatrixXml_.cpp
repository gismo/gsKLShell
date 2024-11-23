#include <gsCore/gsTemplateTools.h>

#include <gsKLShell/src/gsMaterialMatrixXml.h>
#include <gsKLShell/src/gsMaterialMatrixXml.hpp>

namespace gismo
{

// Explicit instantiation

namespace internal
{
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrixBase<real_t>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrixContainer<real_t>>;
} // end namespace internal


} // end namespace gismo
