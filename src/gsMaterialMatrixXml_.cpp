#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXmlUtils.h>

#include <gsKLShell/src/gsMaterialMatrixXml.hpp>
#include <gsKLShell/src/gsMaterialMatrixBase.hpp>
#include <gsKLShell/src/gsMaterialMatrixLinear.hpp>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.hpp>

namespace gismo
{

// Explicit instantiation

namespace internal
{
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrixBase<real_t>>;

    ////////////////////////////////////////////////////////////////
    // 2D compressible
    // NH
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,11,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,21,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,31,true>>;

    // NH_ext
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,12,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,22,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,32,true>>;

    // MR
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,13,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,23,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,33,true>>;

    // OG
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,34,true>>;
    ////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////
    // 3D compressible
    // NH
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,11,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,21,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,31,true>>;

    // NH_ext
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,12,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,22,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,32,true>>;

    // MR
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,13,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,23,true>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,33,true>>;

    // OG
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,34,true>>;
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // 2D incompressible
    // NH
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,11,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,21,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,31,false>>;

    // NH_ext
    // CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,12,false>>;
    // CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,22,false>>;
    // CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,32,false>>;

    // MR
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,13,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,23,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,33,false>>;

    // OG
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<2,real_t,34,false>>;
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // 3D incompressible
    // NH
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,11,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,21,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,31,false>>;

    // NH_ext
    // CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,12,false>>;
    // CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,22,false>>;
    // CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,32,false>>;

    // MR
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,13,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,23,false>>;
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,33,false>>;

    // OG
    CLASS_TEMPLATE_INST gsXml<gsMaterialMatrix<3,real_t,34,false>>;
    ////////////////////////////////////////////////////////////////


} // end namespace internal


} // end namespace gismo
