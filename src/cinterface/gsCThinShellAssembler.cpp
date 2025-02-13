// #ifdef gsCInterface_ENABLED // This does not work....?

#include <gismo.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCTypes.h>

// Original code:
#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellFunctions.h>

#include <gsKLShell/gsThinShellUtils.h>

// C bindings
#include <gsKLShell/cinterface/gsCMaterialMatrix.h>
#include <gsKLShell/cinterface/gsCThinShellAssembler.h>

#define RICAST_MM  reinterpret_cast<gismo::gsMaterialMatrixBase<double> *>
#define RICAST_TSA  reinterpret_cast<gismo::gsThinShellAssemblerBase<double> *>
#define RICAST_CMM  reinterpret_cast<gsCMaterialMatrixBase *>
#define RICAST_CTSA reinterpret_cast<gsCThinShellAssemblerBase *>

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT void gsThinShellAssemblerBase_delete(gsCThinShellAssemblerBase * ptr)
{ delete reinterpret_cast<gismo::gsThinShellAssemblerBase<double>*>(ptr); }

// GISMO_EXPORT void gsThinShellAssemblerBase_print(gsCThinShellAssemblerBase * ptr)
// { reinterpret_cast<gismo::gsThinShellAssemblerBase<double>*>(ptr)->print(gsInfo); }


GISMO_EXPORT gsCThinShellAssemblerBase * gsThinShellAssembler2_create( gsCMultiPatch * mp_ptr,
                                                                        gsCMultiBasis * mb_ptr,
                                                                        gsCBoundaryConditions * bc_ptr,
                                                                        gsCFunctionSet * force_ptr,
                                                                        gsCMaterialMatrixBase * mm_ptr)
{
    auto * mp   = RICAST_MP(mp_ptr);
    auto * mb   = RICAST_MB(mb_ptr);
    auto * bc   = RICAST_BC(bc_ptr);
    auto * force= RICAST_F(force_ptr);
    auto * mm   = RICAST_MM(mm_ptr);
    return (RICAST_CTSA(new gismo::gsThinShellAssembler<2,double,false>(*mp,*mb,*bc,*force,mm)));
}

GISMO_EXPORT gsCThinShellAssemblerBase * gsThinShellAssembler3_bending_create( gsCMultiPatch * mp_ptr,
                                                                                gsCMultiBasis * mb_ptr,
                                                                                gsCBoundaryConditions * bc_ptr,
                                                                                gsCFunctionSet * force_ptr,
                                                                                gsCMaterialMatrixBase * mm_ptr)
{
    auto * mp   = RICAST_MP(mp_ptr);
    auto * mb   = RICAST_MB(mb_ptr);
    auto * bc   = RICAST_BC(bc_ptr);
    auto * force= RICAST_F(force_ptr);
    auto * mm   = RICAST_MM(mm_ptr);
    return (RICAST_CTSA(new gismo::gsThinShellAssembler<3,double,true>(*mp,*mb,*bc,*force,mm)));
}

GISMO_EXPORT gsCThinShellAssemblerBase * gsThinShellAssembler3_nobending_create( gsCMultiPatch * mp_ptr,
                                                                                  gsCMultiBasis * mb_ptr,
                                                                                  gsCBoundaryConditions * bc_ptr,
                                                                                  gsCFunctionSet * force_ptr,
                                                                                  gsCMaterialMatrixBase * mm_ptr)
{
    auto * mp   = RICAST_MP(mp_ptr);
    auto * mb   = RICAST_MB(mb_ptr);
    auto * bc   = RICAST_BC(bc_ptr);
    auto * force= RICAST_F(force_ptr);
    auto * mm   = RICAST_MM(mm_ptr);
    return (RICAST_CTSA(new gismo::gsThinShellAssembler<3,double,false>(*mp,*mb,*bc,*force,mm)));
}

GISMO_EXPORT int gsThinShellAssembler_assemble( gsCThinShellAssemblerBase * assembler_ptr)
{
    auto * assembler = RICAST_TSA(assembler_ptr);
    auto status = assembler->assemble();
    if (status==gismo::ThinShellAssemblerStatus::Success)
        return 0;
    else if (status==gismo::ThinShellAssemblerStatus::AssemblyError)
        return 1;
    else if (status==gismo::ThinShellAssemblerStatus::DimensionError)
        return 2;
    else
        return 99;
}

GISMO_EXPORT int gsThinShellAssembler_assembleMatrix( gsCThinShellAssemblerBase * assembler_ptr, gsCFunctionSet * deformed_ptr)
{
    auto * assembler = RICAST_TSA(assembler_ptr);
    auto * deformed = RICAST_F(deformed_ptr);
    auto status = assembler->assembleMatrix(*deformed);
    if (status==gismo::ThinShellAssemblerStatus::Success)
        return 0;
    else if (status==gismo::ThinShellAssemblerStatus::AssemblyError)
        return 1;
    else if (status==gismo::ThinShellAssemblerStatus::DimensionError)
        return 2;
    else
        return 99;
}

GISMO_EXPORT int gsThinShellAssembler_assembleVector( gsCThinShellAssemblerBase * assembler_ptr, gsCFunctionSet * deformed_ptr)
{
    auto * assembler = RICAST_TSA(assembler_ptr);
    auto * deformed = RICAST_F(deformed_ptr);
    auto status = assembler->assembleVector(*deformed);
    if (status==gismo::ThinShellAssemblerStatus::Success)
        return 0;
    else if (status==gismo::ThinShellAssemblerStatus::AssemblyError)
        return 1;
    else if (status==gismo::ThinShellAssemblerStatus::DimensionError)
        return 2;
    else
        return 99;
}

GISMO_EXPORT void gsThinShellAssembler_matrix_into( gsCThinShellAssemblerBase * assembler_ptr, gsCSparseMatrix * mat_ptr )
{
    auto * assembler = RICAST_TSA(assembler_ptr);
    *RICAST_SM( mat_ptr ) = assembler->matrix(); // THIS LINE IS PROBLEMATIC
    // mat_ptr = RICAST_CSM(new gismo::gsSparseMatrix<double>(assembler->matrix()));
}

GISMO_EXPORT void gsThinShellAssembler_rhs_into( gsCThinShellAssemblerBase * assembler_ptr, gsCMatrix * rhs_ptr )
{
    *RICAST_M(rhs_ptr) = RICAST_TSA(assembler_ptr)->rhs();
}

GISMO_EXPORT gsCMultiPatch * gsThinShellAssembler_constructSolution( gsCThinShellAssemblerBase * assembler_ptr, gsCMatrix * solVector_ptr)
{
    auto * assembler = RICAST_TSA(assembler_ptr);
    auto * solVector = RICAST_M(solVector_ptr);
    return RICAST_CMP(new gismo::gsMultiPatch<double>(assembler->constructSolution(*solVector)));
}

#ifdef __cplusplus
}
#endif

// #endif // gsCInterface_ENABLED