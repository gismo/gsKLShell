/** @file gsMaterialMatrix_test.cpp

    @brief Code for the arc-length method of a shell based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#include <gsKLShell/src/gsThinShellUtils.h>

using namespace gismo;

int main (int argc, char** argv)
{

    int material = 0;
	int impl = 1;
	int Compressibility = 0;
	gsCmdLine cmd("Thin shell plate example.");

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "I", "Impl", "Implementation",  impl );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsMultiPatch<> mp, mp_def;
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);

    mp_def = mp;
    mp_def.patch(0).coefs().col(0) *= 2;
    // mp_def.patch(0).coefs() *= 2;
    mp_def.patch(0).coefs()(3,0) *= 2;

    real_t thickness = 1e-2;
    real_t E_modulus = 1.0;
    real_t PoissonRatio = (Compressibility==1 || material==0) ? 0.45 : 0.5;
    real_t Ratio = 2;
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsConstantFunction<> ratio(Ratio,3);

    std::vector<gsFunctionSet<>*> parameters;
    if (material==0 || material==1 || material==2)
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
    }
    else if (material==3)
    {
        parameters.resize(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
    }
    else
        GISMO_ERROR("Material unknown");

    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);


    gsVector<> pt(2);
    pt.col(0)<<0.125,0.375;

    gsMatrix<> pts(2,6);
    pts.col(0)<<0.125,0.375;
    pts.col(1)<<0.375,0.125;
    pts.col(2)<<0.125,0.25;
    pts.col(3)<<0.25,0.125;
    pts.col(4)<<0.25,0.25;
    pts.col(5)<<0.5,0.5;

    gsMatrix<> result;

    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;

    geometryMap map = A.getMap(mp); // the last map counts
    geometryMap def = A.getMap(mp_def); // the last map counts

    gsExprEvaluator<> ev(A);


    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////Integrals///////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"______________________________Integrals________________________________\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> mmA_i(materialMatrix,&mp_def);
    variable mmAi = A.getCoeff(mmA_i);
    gsInfo<<"matrix A = \n"<<ev.eval(mmAi,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixB> mmB_i(materialMatrix,&mp_def);
    variable mmBi = A.getCoeff(mmB_i);
    gsInfo<<"matrix B = \n"<<ev.eval(mmBi,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixC> mmC_i(materialMatrix,&mp_def);
    variable mmCi = A.getCoeff(mmC_i);
    gsInfo<<"matrix C = \n"<<ev.eval(mmCi,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixD> mmD_i(materialMatrix,&mp_def);
    variable mmDi = A.getCoeff(mmD_i);
    gsInfo<<"matrix D = \n"<<ev.eval(mmDi,pt)<<"\n";


    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> S0_i(materialMatrix,&mp_def);
    variable S0i = A.getCoeff(S0_i);
    gsInfo<<"Vector N = \n"<<ev.eval(S0i,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> S1_i(materialMatrix,&mp_def);
    variable S1i = A.getCoeff(S1_i);
    gsInfo<<"Vector M = \n"<<ev.eval(S1i,pt)<<"\n";

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////Point values////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"______________________________Point values_____________________________\n";
    gsMatrix<> Z0(1,1);
    Z0.setZero();
    gsMatrix<> Zf(1,1);
    Zf.setOnes();

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> mmA_p(materialMatrix,&mp_def,Z0);
    variable mmAp = A.getCoeff(mmA_p);
    gsInfo<<"matrix A = \n"<<ev.eval(mmAp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> mmB_p(materialMatrix,&mp_def,Z0);
    variable mmBp = A.getCoeff(mmB_p);
    gsInfo<<"matrix B = \n"<<ev.eval(mmBp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> mmC_p(materialMatrix,&mp_def,Z0);
    variable mmCp = A.getCoeff(mmC_p);
    gsInfo<<"matrix C = \n"<<ev.eval(mmCp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> mmD_p(materialMatrix,&mp_def,Z0);
    variable mmDp = A.getCoeff(mmD_p);
    gsInfo<<"matrix D = \n"<<ev.eval(mmDp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> S0_p(materialMatrix,&mp_def,Z0);
    variable S0p = A.getCoeff(S0_p);
    gsInfo<<"Vector N = \n"<<ev.eval(S0p,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> S1_p(materialMatrix,&mp_def,Zf);
    variable S1p = A.getCoeff(S1_p);
    gsInfo<<"Vector M = \n"<<ev.eval(S1p,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::Generic> St_p(materialMatrix,&mp_def,Zf);
    variable Stp = A.getCoeff(St_p);
    gsInfo<<"Vector total = \n"<<ev.eval(Stp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::PStressN> P0_p(materialMatrix,&mp_def,Z0);
    variable P0p = A.getCoeff(P0_p);
    gsInfo<<"Pstress N = \n"<<ev.eval(P0p,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::PStressM> P1_p(materialMatrix,&mp_def,Zf);
    variable P1p = A.getCoeff(P1_p);
    gsInfo<<"Pstress M = \n"<<ev.eval(P1p,pt)<<"\n";


    gsMaterialMatrixEval<real_t,MaterialOutput::Stretch> lambda_p(materialMatrix,&mp_def,Z0);
    variable lambdap = A.getCoeff(lambda_p);
    gsInfo<<"Stretch = \n"<<ev.eval(lambdap,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir> lambdaDir_p(materialMatrix,&mp_def,Z0);
    variable lambdaDirp = A.getCoeff(lambdaDir_p);
    gsInfo<<"Stretch dirs = \n"<<ev.eval(lambdaDirp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::Spec2CovTransform> trans_p1(materialMatrix,&mp_def,Z0);
    variable transp1 = A.getCoeff(trans_p1);
    gsInfo<<"Spec2CovTransform = \n"<<ev.eval(transp1,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::Spec2ConTransform> trans_p2(materialMatrix,&mp_def,Z0);
    variable transp2 = A.getCoeff(trans_p2);
    gsInfo<<"Spec2ConTransform = \n"<<ev.eval(transp2,pt)<<"\n";

	return 0;
}