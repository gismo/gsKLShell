#ifdef __cplusplus
extern "C"
{
#endif

    struct gsCThinShellAssemblerBase; // Opaque type that we use as a handle
    typedef struct gsCThinShellAssemblerBase gsCThinShellAssemblerBase;

    GISMO_EXPORT void gsThinShellAssemblerBase_delete(gsCThinShellAssemblerBase * ptr);

    // GISMO_EXPORT void gsThinShellAssemblerBase_print(gsCThinShellAssemblerBase * ptr);

    GISMO_EXPORT gsCThinShellAssemblerBase * gsThinShellAssembler2_create( gsCMultiPatch * mp_ptr,
                                                                            gsCMultiBasis * t_ptr,
                                                                            gsCBoundaryConditions * bc_ptr,
                                                                            gsCFunctionSet * f_ptr,
                                                                            gsCMaterialMatrixBase * m_ptr);

    GISMO_EXPORT gsCThinShellAssemblerBase * gsThinShellAssembler3_bending_create( gsCMultiPatch * mp_ptr,
                                                                                    gsCMultiBasis * t_ptr,
                                                                                    gsCBoundaryConditions * bc_ptr,
                                                                                    gsCFunctionSet * f_ptr,
                                                                                    gsCMaterialMatrixBase * m_ptr);


    GISMO_EXPORT int gsThinShellAssembler_assemble( gsCThinShellAssemblerBase * assembler_ptr);

    GISMO_EXPORT int gsThinShellAssembler_assembleMatrix( gsCThinShellAssemblerBase * assembler_ptr, gsCFunctionSet * deformed_ptr);

    GISMO_EXPORT int gsThinShellAssembler_assembleVector( gsCThinShellAssemblerBase * assembler_ptr, gsCFunctionSet * deformed_ptr);

    GISMO_EXPORT void gsThinShellAssembler_matrix_into( gsCThinShellAssemblerBase * assembler_ptr, gsCSparseMatrix * mat_ptr );

    GISMO_EXPORT void gsThinShellAssembler_rhs_into( gsCThinShellAssemblerBase * assembler_ptr, gsCMatrix * rhs_ptr );

    GISMO_EXPORT gsCMultiPatch * gsThinShellAssembler_constructSolution( gsCThinShellAssemblerBase * assembler_ptr, gsCMatrix * solVector_ptr);


#ifdef __cplusplus
}
#endif