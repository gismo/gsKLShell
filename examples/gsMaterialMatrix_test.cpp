/** @file gsMaterialMatrix_test.cpp

    @brief Code for the arc-length method of a shell based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsThinShellUtils.h>

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
    gsReadFile<>("deformed_plate_lin_T=" + std::to_string(1.0) + ".xml",mp_def);



    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);


    // gsMatrix<> ones;
    // // translate second patch up
    // ones = gsMatrix<>::Ones(mp.patch(0).coefs().rows(),1); // patch 1
    // // mp_def.patch(0).coefs().col(1) += ones; // patch 1
    // mp_def.patch(0).scale(4,0); // patch 1
    // mp_def.patch(0).scale(2,1); // patch 1
    // mp_def.patch(0).rotate(3.141592/4.); // patch 1

    // mp.embed(3);
    // mp_def.embed(3);


    real_t thickness = 1.0;
    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.499;
    real_t Ratio = 2;
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsConstantFunction<> ratio(Ratio,3);


    std::vector<gsFunction<>*> parameters(3);
    parameters[0] = &E;
    parameters[1] = &nu;
    parameters[2] = &ratio;

    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);


    gsVector<> pt(2);
    pt.setConstant(0.25);
    pt<<0.125,0.375;

    gsMatrix<> result;

    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;

    geometryMap map = A.getMap(mp); // the last map counts

    gsExprEvaluator<> ev(A);


    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////Integrals///////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"______________________________Integrals________________________________\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> mmA_i(materialMatrix,mp_def);
    variable mmAi = A.getCoeff(mmA_i);
    gsInfo<<"matrix A = \n"<<ev.eval(reshape(mmAi,3,3),pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixB> mmB_i(materialMatrix,mp_def);
    variable mmBi = A.getCoeff(mmB_i);
    gsInfo<<"matrix B = \n"<<ev.eval(reshape(mmBi,3,3),pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixC> mmC_i(materialMatrix,mp_def);
    variable mmCi = A.getCoeff(mmC_i);
    gsInfo<<"matrix C = \n"<<ev.eval(reshape(mmCi,3,3),pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixD> mmD_i(materialMatrix,mp_def);
    variable mmDi = A.getCoeff(mmD_i);
    gsInfo<<"matrix D = \n"<<ev.eval(reshape(mmDi,3,3),pt)<<"\n";



    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> S0_i(materialMatrix,mp_def);
    variable S0i = A.getCoeff(S0_i);
    gsInfo<<"Vector N = \n"<<ev.eval(reshape(S0i,3,1),pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> S1_i(materialMatrix,mp_def);
    variable S1i = A.getCoeff(S1_i);
    gsInfo<<"Vector M = \n"<<ev.eval(reshape(S1i,3,1),pt)<<"\n";


    gsMaterialMatrixIntegrate<real_t,MaterialOutput::PStressN> P0_i(materialMatrix,mp_def);
    variable P0i = A.getCoeff(P0_i);
    gsInfo<<"Pstress N = \n"<<ev.eval(reshape(P0i,2,1),pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::PStressM> P1_i(materialMatrix,mp_def);
    variable P1i = A.getCoeff(P1_i);
    gsInfo<<"Pstress M = \n"<<ev.eval(reshape(P1i,2,1),pt)<<"\n";


    gsMaterialMatrixIntegrate<real_t,MaterialOutput::Stretch> lambda_i(materialMatrix,mp_def);
    variable lambdai = A.getCoeff(lambda_i);
    gsInfo<<"Stretch = \n"<<ev.eval(reshape(lambdai,3,1),pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::StretchDir> lambdaDir_i(materialMatrix,mp_def);
    variable lambdaDiri = A.getCoeff(lambdaDir_i);
    gsInfo<<"Stretch dirs = \n"<<ev.eval(reshape(lambdaDiri,3,3),pt)<<"\n";

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////Point values////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"______________________________Point values_____________________________\n";
    gsMatrix<> Z(1,1);
    Z.setZero();

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> mmA_p(materialMatrix,mp_def,Z);
    variable mmAp = A.getCoeff(mmA_p);
    gsInfo<<"matrix A = \n"<<ev.eval(reshape(mmAp,3,3),pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> mmB_p(materialMatrix,mp_def,Z);
    variable mmBp = A.getCoeff(mmB_p);
    gsInfo<<"matrix B = \n"<<ev.eval(reshape(mmBp,3,3),pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> mmC_p(materialMatrix,mp_def,Z);
    variable mmCp = A.getCoeff(mmC_p);
    gsInfo<<"matrix C = \n"<<ev.eval(reshape(mmCp,3,3),pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> mmD_p(materialMatrix,mp_def,Z);
    variable mmDp = A.getCoeff(mmD_p);
    gsInfo<<"matrix D = \n"<<ev.eval(reshape(mmDp,3,3),pt)<<"\n";



    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> S0_p(materialMatrix,mp_def,Z);
    variable S0p = A.getCoeff(S0_p);
    gsInfo<<"Vector N = \n"<<ev.eval(reshape(S0i,3,1),pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> S1_p(materialMatrix,mp_def,Z);
    variable S1p = A.getCoeff(S1_p);
    gsInfo<<"Vector M = \n"<<ev.eval(reshape(S1p,3,1),pt)<<"\n";

    /// FIX THIS ONE
    gsMaterialMatrixEval<real_t,MaterialOutput::Generic> St_p(materialMatrix,mp_def,Z);
    variable Stp = A.getCoeff(St_p);
    gsInfo<<"Vector total = \n"<<ev.eval(reshape(Stp,3,1),pt)<<"\n";


    gsMaterialMatrixEval<real_t,MaterialOutput::PStressN> P0_p(materialMatrix,mp_def,Z);
    variable P0p = A.getCoeff(P0_p);
    gsInfo<<"Pstress N = \n"<<ev.eval(reshape(P0p,2,1),pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::PStressM> P1_p(materialMatrix,mp_def,Z);
    variable P1p = A.getCoeff(P1_p);
    gsInfo<<"Pstress M = \n"<<ev.eval(reshape(P1p,2,1),pt)<<"\n";


    gsMaterialMatrixEval<real_t,MaterialOutput::Stretch> lambda_p(materialMatrix,mp_def,Z);
    variable lambdap = A.getCoeff(lambda_p);
    gsInfo<<"Stretch = \n"<<ev.eval(reshape(lambdap,3,1),pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir> lambdaDir_p(materialMatrix,mp_def,Z);
    variable lambdaDirp = A.getCoeff(lambdaDir_p);
    gsInfo<<"Stretch dirs = \n"<<ev.eval(reshape(lambdaDirp,3,3),pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::Transformation> trans_p(materialMatrix,mp_def,Z);
    variable transp = A.getCoeff(trans_p);
    gsInfo<<"Transformation = \n"<<ev.eval(reshape(transp,3,3),pt)<<"\n";


    gsInfo<<"PStressN (transformed) = \n"<<ev.eval(reshape(transp,3,3) * S0p,pt)<<"\n";
    gsInfo<<"PStressM (transformed) = \n"<<ev.eval(reshape(transp,3,3) * S1p,pt)<<"\n";
    gsInfo<<"PStressM (transformed) = \n"<<ev.eval(reshape(transp,3,3) * Stp,pt)<<"\n";



    //////////////////////////////////////////////////////////////////////////////////////////////




    gsWriteParaview(mp,"mp",1000);
    gsWriteParaview(mp_def,"mp_def",1000);

	return 0;
}