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

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> mmA_i(materialMatrix,mp_def);
    variable mmAi = A.getCoeff(mmA_i);
    gsInfo<<"matrix A = \n"<<ev.eval(mmAi,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixB> mmB_i(materialMatrix,mp_def);
    variable mmBi = A.getCoeff(mmB_i);
    gsInfo<<"matrix B = \n"<<ev.eval(mmBi,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixC> mmC_i(materialMatrix,mp_def);
    variable mmCi = A.getCoeff(mmC_i);
    gsInfo<<"matrix C = \n"<<ev.eval(mmCi,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixD> mmD_i(materialMatrix,mp_def);
    variable mmDi = A.getCoeff(mmD_i);
    gsInfo<<"matrix D = \n"<<ev.eval(mmDi,pt)<<"\n";


    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> S0_i(materialMatrix,mp_def);
    variable S0i = A.getCoeff(S0_i);
    gsInfo<<"Vector N = \n"<<ev.eval(S0i,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> S1_i(materialMatrix,mp_def);
    variable S1i = A.getCoeff(S1_i);
    gsInfo<<"Vector M = \n"<<ev.eval(S1i,pt)<<"\n";


    gsMaterialMatrixIntegrate<real_t,MaterialOutput::PStressN> P0_i(materialMatrix,mp_def);
    variable P0i = A.getCoeff(P0_i);
    gsInfo<<"Pstress N = \n"<<ev.eval(P0i,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::PStressM> P1_i(materialMatrix,mp_def);
    variable P1i = A.getCoeff(P1_i);
    gsInfo<<"Pstress M = \n"<<ev.eval(P1i,pt)<<"\n";


    gsMaterialMatrixIntegrate<real_t,MaterialOutput::Stretch> lambda_i(materialMatrix,mp_def);
    variable lambdai = A.getCoeff(lambda_i);
    gsInfo<<"Stretch = \n"<<ev.eval(lambdai,pt)<<"\n";

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::StretchDir> lambdaDir_i(materialMatrix,mp_def);
    variable lambdaDiri = A.getCoeff(lambdaDir_i);
    gsInfo<<"Stretch dirs = \n"<<ev.eval(lambdaDiri,pt)<<"\n";

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////Point values////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"______________________________Point values_____________________________\n";
    gsMatrix<> Z(1,1);
    Z.setZero();

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> mmA_p(materialMatrix,mp_def,Z);
    variable mmAp = A.getCoeff(mmA_p);
    gsInfo<<"matrix A = \n"<<ev.eval(mmAp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> mmB_p(materialMatrix,mp_def,Z);
    variable mmBp = A.getCoeff(mmB_p);
    gsInfo<<"matrix B = \n"<<ev.eval(mmBp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> mmC_p(materialMatrix,mp_def,Z);
    variable mmCp = A.getCoeff(mmC_p);
    gsInfo<<"matrix C = \n"<<ev.eval(mmCp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> mmD_p(materialMatrix,mp_def,Z);
    variable mmDp = A.getCoeff(mmD_p);
    gsInfo<<"matrix D = \n"<<ev.eval(mmDp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> S0_p(materialMatrix,mp_def,Z);
    variable S0p = A.getCoeff(S0_p);
    gsInfo<<"Vector N = \n"<<ev.eval(S0p,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> S1_p(materialMatrix,mp_def,Z);
    variable S1p = A.getCoeff(S1_p);
    gsInfo<<"Vector M = \n"<<ev.eval(S1p,pt)<<"\n";

    /// FIX THIS ONE
    gsMaterialMatrixEval<real_t,MaterialOutput::Generic> St_p(materialMatrix,mp_def,Z);
    variable Stp = A.getCoeff(St_p);
    gsInfo<<"Vector total = \n"<<ev.eval(Stp,pt)<<"\n";


    gsMaterialMatrixEval<real_t,MaterialOutput::PStressN> P0_p(materialMatrix,mp_def,Z);
    variable P0p = A.getCoeff(P0_p);
    gsInfo<<"Pstress N = \n"<<ev.eval(P0p,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::PStressM> P1_p(materialMatrix,mp_def,Z);
    variable P1p = A.getCoeff(P1_p);
    gsInfo<<"Pstress M = \n"<<ev.eval(P1p,pt)<<"\n";


    gsMaterialMatrixEval<real_t,MaterialOutput::Stretch> lambda_p(materialMatrix,mp_def,Z);
    variable lambdap = A.getCoeff(lambda_p);
    gsInfo<<"Stretch = \n"<<ev.eval(lambdap,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir> lambdaDir_p(materialMatrix,mp_def,Z);
    variable lambdaDirp = A.getCoeff(lambdaDir_p);
    gsInfo<<"Stretch dirs = \n"<<ev.eval(lambdaDirp,pt)<<"\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::Transformation> trans_p(materialMatrix,mp_def,Z);
    variable transp = A.getCoeff(trans_p);
    gsInfo<<"Transformation = \n"<<ev.eval(transp,pt)<<"\n";

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    variable m2 = A.getCoeff(mult2t);

    auto Em = 0.5 * ( flat(jac(def).tr()*jac(def)) - flat(jac(map).tr()* jac(map)) );
    auto Cm = ( flat(jac(def).tr()*jac(def)) ) * reshape(m2,3,3);
    gsInfo<<"Em = \n"<<ev.eval(Em,pt)<<"\n";
    gsInfo<<"Sm = \n"<<ev.eval(Em * reshape(mmAp,3,3),pt)<<"\n";
    gsInfo<<"Em = \n"<<ev.eval(Cm * reshape(transp,3,3).tr(),pt)<<"\n";

    gsInfo<<"PStressN (transformed) = \n"<<ev.eval(reshape(transp,3,3) * S0p,pt)<<"\n";
    gsInfo<<"PStressM (transformed) = \n"<<ev.eval(reshape(transp,3,3) * S1p,pt)<<"\n";
    gsInfo<<"PStressM (transformed) = \n"<<ev.eval(reshape(transp,3,3) * Stp,pt)<<"\n";



    //////////////////////////////////////////////////////////////////////////////////////////////

    gsMatrix<> Tmats;
    gsMatrix<> Stretches;
    gsMatrix<> Stresses;
    gsMatrix<> PStresses;
    trans_p.eval_into(pts,Tmats);
    lambda_p.eval_into(pts,Stretches);
    S0_p.eval_into(pts,Stresses);
    P0_p.eval_into(pts,PStresses);
    for (index_t k=0; k!=pts.cols(); k++)
    {
        gsMatrix<> Tmat = Tmats.reshapeCol(k,3,3);
        gsMatrix<> Stretch = Stretches.reshapeCol(k,3,1);
        gsMatrix<> Stress = Stresses.reshapeCol(k,3,1);
        gsMatrix<> PStress = PStresses.reshapeCol(k,3,1);

        gsDebugVar((PStress - Tmat * Stress).norm());
        gsDebugVar((PStress.transpose() - Stress.transpose() * Tmat.transpose()).norm());
    }


    gsWriteParaview(mp,"mp",1000);
    gsWriteParaview(mp_def,"mp_def",1000);

	return 0;
}