// #ifdef gsCInterface_ENABLED // This does not work....?

#include <gismo.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCTypes.h>

// Original code:
#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>

// C bindings
#include <gsKLShell/cinterface/gsCMaterialMatrix.h>

#define RICAST_CMM reinterpret_cast<gsCMaterialMatrixBase *>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT void gsMaterialMatrixBase_delete(gsCMaterialMatrixBase * ptr)
{ delete reinterpret_cast<gsMaterialMatrixBase<double>*>(ptr); }

GISMO_EXPORT void gsMaterialMatrixBase_print(gsCMaterialMatrixBase * ptr)
{ reinterpret_cast<gsMaterialMatrixBase<double>*>(ptr)->print(gsInfo); }

GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear2_create(gsCFunctionSet * mp_ptr,
                                                                   gsCFunctionSet * t_ptr,
                                                                   gsCFunctionSet * E_ptr,
                                                                   gsCFunctionSet * nu_ptr,
                                                                   gsCFunctionSet * rho_ptr)
{
    auto * mp  = RICAST_F(mp_ptr);
    auto * t   = RICAST_F(t_ptr);
    auto * E   = RICAST_F(E_ptr);
    auto * nu  = RICAST_F(nu_ptr);
    auto * rho = RICAST_F(rho_ptr);
    return (RICAST_CMM(new gsMaterialMatrixLinear<2,double>(*mp,*t,*E,*nu,*rho)));
}

GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear2const_create( gsCFunctionSet * mp_ptr,
                                                                     double t, double E, double nu, double rho)
{
    auto * mp  = RICAST_F(mp_ptr);
    gsConstantFunction<double> tfun(t,2);
    gsConstantFunction<double> Efun(t,2);
    gsConstantFunction<double> nfun(t,2);
    gsConstantFunction<double> rfun(t,2);
    return (RICAST_CMM(new gsMaterialMatrixLinear<2,double>(*mp,tfun,Efun,nfun,rfun)));
}


GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear3_create(gsCFunctionSet * mp_ptr,
                                                                   gsCFunctionSet * t_ptr,
                                                                   gsCFunctionSet * E_ptr,
                                                                   gsCFunctionSet * nu_ptr,
                                                                   gsCFunctionSet * rho_ptr)
{
    auto * mp  = RICAST_F(mp_ptr);
    auto * t   = RICAST_F(t_ptr);
    auto * E   = RICAST_F(E_ptr);
    auto * nu  = RICAST_F(nu_ptr);
    auto * rho = RICAST_F(rho_ptr);
    return (RICAST_CMM(new gsMaterialMatrixLinear<3,double>(*mp,*t,*E,*nu,*rho)));
}

GISMO_EXPORT gsCMaterialMatrixBase * gsMaterialMatrixLinear3const_create( gsCFunctionSet * mp_ptr,
                                                                     double t, double E, double nu, double rho)
{
    auto * mp  = RICAST_F(mp_ptr);
    gsConstantFunction<double> tfun(t,3);
    gsConstantFunction<double> Efun(t,3);
    gsConstantFunction<double> nfun(t,3);
    gsConstantFunction<double> rfun(t,3);
    return (RICAST_CMM(new gsMaterialMatrixLinear<3,double>(*mp,tfun,Efun,nfun,rfun)));
}

#ifdef __cplusplus
}
#endif

// #endif // gsCInterface_ENABLED