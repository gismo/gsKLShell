#ifdef __cplusplus
extern "C"
{
#endif

    /* BASE CLASS */
    struct gsCMaterialMatrixBase; // Opaque type that we use as a handle
    typedef struct gsCMaterialMatrixBase gsCMaterialMatrixBase;

    GISMO_EXPORT void gsMaterialMatrixBase_delete(gsCMaterialMatrixBase * ptr);

    GISMO_EXPORT void gsMaterialMatrixBase_print(gsCMaterialMatrixBase * ptr);

    /* DERIVED CLASSES */
    // gsMaterialMatrixLinear
#   define gsMaterialMatrixLinear_print gsMaterialMatrixBase_print
#   define gsMaterialMatrixLinear_delete gsMaterialMatrixBase_delete

    GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear2_create( gsCFunctionSet * mp_ptr,
                                                                         gsCFunctionSet * t_ptr,
                                                                         gsCFunctionSet * E_ptr,
                                                                         gsCFunctionSet * nu_ptr,
                                                                         gsCFunctionSet * rho_ptr);

    GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear2const_create( gsCFunctionSet * mp_ptr,
                                                                         double t, double E, double nu, double rho);

    GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear3_create( gsCFunctionSet * mp_ptr,
                                                                         gsCFunctionSet * t_ptr,
                                                                         gsCFunctionSet * E_ptr,
                                                                         gsCFunctionSet * nu_ptr,
                                                                         gsCFunctionSet * rho_ptr);

    GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear3const_create( gsCFunctionSet * mp_ptr,
                                                                         double t, double E, double nu, double rho);


#ifdef __cplusplus
}
#endif