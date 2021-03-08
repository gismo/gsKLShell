/** @file gsThinShell_test2.cpp

    @brief Example testing and debugging thin shell solver. Based on gsThinShell_test

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/getMaterialMatrix.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool Compressibility = 0;
    index_t material = 0;
    bool verbose = false;
    std::string fn;
    bool membrane = false;


    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;
    real_t Ratio = 7.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addReal( "R", "Ratio", "Mooney Rivlin Ratio",  Ratio );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addSwitch( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.embed(3);
    mp.addAutoBoundaries();
    E_modulus = 1e0;
    thickness = 1e0;
    PoissonRatio = 0.3;

    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);


    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    gsConstantFunction<> displx(0.1,3);
    gsConstantFunction<> disply(0.25,3);

    gsFunctionExpr<> neuDataFun1;
    gsConstantFunction<> neuData(neu,3);
    real_t pressure = 0.0;

    // neu <<0.1,0,0;
    // neuData.setValue(neu,3);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,0);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,0);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,0);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,0);

    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,1);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,1);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,1);

    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,2);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,2);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,2);

    tmp << 0,0,-1;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsConstantFunction<> pressFun(pressure,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(Ratio,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));

    gsConstantFunction<> alpha1(2.0,3);
    gsConstantFunction<> mu1(7.0*mu/8.0,3);
    gsConstantFunction<> alpha2(-2.0,3);
    gsConstantFunction<> mu2(-mu/8.0,3);
    // gsConstantFunction<> alpha3()
    // gsConstantFunction<> mu3()

    // gsMaterialMatrix materialMatrixNonlinear(mp,mp_def,t,E,nu,rho);
    std::vector<gsFunction<>*> parameters(3);
    parameters[0] = &E;
    parameters[1] = &nu;
    parameters[2] = &ratio;
    gsMaterialMatrixBase<real_t>* materialMatrixNonlinear;
    gsMaterialMatrixBase<real_t>* materialMatrixLinear;
    gsMaterialMatrixBase<real_t>* materialMatrixComposite;

    materialMatrixLinear = new gsMaterialMatrixLinear<3,real_t>(mp,mp_def,t,parameters,rho);

    // Linear anisotropic material model
    real_t pi = math::atan(1)*4;

    index_t kmax = 2;
    gsVector<> E11(kmax), E22(kmax), G12(kmax), nu12(kmax), nu21(kmax), thick(kmax), phi(kmax);
    E11.setZero(); E22.setZero(); G12.setZero(); nu12.setZero(); nu21.setZero(); thick.setZero(); phi.setZero();
    for (index_t k=0; k != kmax; ++k)
    {
        E11.at(k) = E22.at(k) = E_modulus;
        nu12.at(k) = nu21.at(k) = PoissonRatio;
        G12.at(k) = 0.5 * E_modulus / (1+PoissonRatio);
        thick.at(k) = thickness/kmax;
        phi.at(kmax) = 0.* k / kmax * pi/2.0;
    }
    std::vector<gsFunction<>*> compParameters(6);

    gsConstantFunction<> E11fun(E11,3);
    gsConstantFunction<> E22fun(E22,3);
    gsConstantFunction<> G12fun(G12,3);
    gsConstantFunction<> nu12fun(nu12,3);
    gsConstantFunction<> nu21fun(nu21,3);
    gsConstantFunction<> thickfun(thick,3);
    gsConstantFunction<> phifun(phi,3);
    compParameters[0] = &E11fun;
    compParameters[1] = &E22fun;
    compParameters[2] = &G12fun;
    compParameters[3] = &nu12fun;
    compParameters[4] = &nu21fun;
    compParameters[5] = &phifun;
    materialMatrixComposite = new gsMaterialMatrixComposite<3,real_t>(mp,mp_def,thickfun,compParameters);

    gsOptionList optionsNL;
    optionsNL.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
    optionsNL.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
    optionsNL.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",3);
    materialMatrixNonlinear = getMaterialMatrix<3,real_t>(mp,mp_def,t,parameters,rho,optionsNL);


    gsOptionList optionsL;
    optionsL.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    materialMatrixLinear = getMaterialMatrix<3,real_t>(mp,mp_def,t,parameters,rho,optionsL);

    gsOptionList optionsC;
    optionsC.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    optionsC.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",0);
    materialMatrixComposite = getMaterialMatrix<3,real_t>(mp,mp_def,t,compParameters,rho,optionsC);

    std::vector<gsFunction<>*> parameters2(6);
    if (material==4)
    {
        parameters2[0] = &E;
        parameters2[1] = &nu;
        parameters2[2] = &mu1;
        parameters2[3] = &alpha1;

        parameters2[4] = &mu2;
        parameters2[5] = &alpha2;

        // parameters[6] = ;
        // parameters[7] = ;
        materialMatrixNonlinear->setParameters(parameters2);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    if(membrane)
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,bc,force,materialMatrixNonlinear);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrixNonlinear);

    assembler->setPointLoads(pLoads);
    if (pressure!= 0.0)
        assembler->setPressure(pressFun);

    assembler->assemble();

    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();

    gsDebugVar(assembler->matrix().toDense());
    gsDebugVar(assembler->rhs().transpose());

    // Solve linear problem
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    gsDebugVar(solVector.transpose());

    mp_def = assembler->constructSolution(solVector);

    gsMatrix<> pts(2,3);
    pts.col(0).setConstant(0.25);
    pts.col(1).setConstant(0.50);
    pts.col(2).setConstant(0.75);

    gsMatrix<> result;


    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> materialMatrixEvalA(materialMatrixNonlinear);
    materialMatrixEvalA.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> materialMatrixEvalB(materialMatrixNonlinear);
    materialMatrixEvalB.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> materialMatrixEvalC(materialMatrixNonlinear);
    materialMatrixEvalC.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> materialMatrixEvalD(materialMatrixNonlinear);
    materialMatrixEvalD.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> materialVectorEvalN(materialMatrixNonlinear);
    materialVectorEvalN.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> materialVectorEvalM(materialMatrixNonlinear);
    materialVectorEvalM.eval_into(pts,result);
    gsDebugVar(result);

    gsDebug<<"======================================================\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> materialMatrixEvalA_lin(materialMatrixLinear);
    materialMatrixEvalA_lin.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> materialMatrixEvalB_lin(materialMatrixLinear);
    materialMatrixEvalB_lin.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> materialMatrixEvalC_lin(materialMatrixLinear);
    materialMatrixEvalC_lin.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> materialMatrixEvalD_lin(materialMatrixLinear);
    materialMatrixEvalD_lin.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> materialVectorEvalN_lin(materialMatrixLinear);
    materialVectorEvalN_lin.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> materialVectorEvalM_lin(materialMatrixLinear);
    materialVectorEvalM_lin.eval_into(pts,result);
    gsDebugVar(result);

    gsDebug<<"======================================================\n";

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> materialMatrixEvalA_com(materialMatrixComposite);
    materialMatrixEvalA_com.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> materialMatrixEvalB_com(materialMatrixComposite);
    materialMatrixEvalB_com.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> materialMatrixEvalC_com(materialMatrixComposite);
    materialMatrixEvalC_com.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> materialMatrixEvalD_com(materialMatrixComposite);
    materialMatrixEvalD_com.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> materialVectorEvalN_com(materialMatrixComposite);
    materialVectorEvalN_com.eval_into(pts,result);
    gsDebugVar(result);

    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> materialVectorEvalM_com(materialMatrixComposite);
    materialVectorEvalM_com.eval_into(pts,result);
    gsDebugVar(result);

delete assembler;
    return 0;

}// end main

