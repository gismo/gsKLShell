/** @file example_shell3D.cpp

    @brief Examples for the surface thin shells including the shell obstacle course

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsThinShellUtils.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Number of adaptive refinement loops
    index_t RefineLoopMax;
    // Flag for refinemet criterion
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    index_t refCriterion;
    // Parameter for computing adaptive refinement threshold
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    real_t refParameter;  // ...specified below with the examples


    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t goal = 1;
    bool nonlinear = false;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t thickness = 1.0;

    refCriterion = GARU;
    refParameter = 0.85;
    RefineLoopMax = 1;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt("R", "refine", "Maximum number of adaptive refinement steps to perform",
        RefineLoopMax);
    cmd.addReal( "T", "thickness", "thickness",  thickness );
    cmd.addInt( "g", "goal", "Goal function to use", goal );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsVector<> pts(2);
    pts.setConstant(0.25);

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    gsMultiPatch<> mp_ex;
    gsReadFile<>("deformed_plate_T=" + std::to_string(thickness) + ".xml",mp_ex);
    gsMultiBasis<> basisR(mp_ex);

    // Unit square
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);
    E_modulus = 1.0;
    // thickness = 1.0;

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    gsMultiBasis<> dbasis(mp);
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    // // Cast all patches of the mp object to THB splines
    // gsTHBSpline<2,real_t> thb;
    // for (index_t k=0; k!=mp.nPatches(); ++k)
    // {
    //     gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
    //     thb = gsTHBSpline<2,real_t>(*geo);
    //     mp.patch(k) = thb;
    // }

    mp_def = mp;

    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH = basisL;
    basisH.degreeElevate(1);

    // gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n";
    // gsInfo<<"Basis Dual:   "<<basisH.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    real_t load = 1.0;
    real_t D = E_modulus * math::pow(thickness,3) / ( 12 * ( 1- math::pow(PoissonRatio,2) ) );

    gsFunctionExpr<> u_ex( "0","0","w:= 0; for (u := 1; u < 100; u += 2) { for (v := 1; v < 100; v += 2) { w += -16.0 * " + std::to_string(load) + " / ( pi^6*" + std::to_string(D) + " ) * 1 / (v * u * ( v^2 + u^2 )^2 ) * sin( v * pi * x) * sin(u * pi * y) } }",3);
    gsFunctionExpr<> z_ex( "0","0","w:= 0; for (u := 1; u < 100; u += 2) { for (v := 1; v < 100; v += 2) { w += 16.0 * 1 / ( pi^6*" + std::to_string(D) + " ) * 1 / (v * u * ( v^2 + u^2 )^2 ) * sin( v * pi * x) * sin(u * pi * y) } }",3);

    for (index_t i=0; i!=3; ++i)
    {
        bc.addCondition(boundary::north,condition_type::dirichlet, 0, i );
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, i );
        bc.addCondition(boundary::south,condition_type::dirichlet, 0, i );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, i );
    }
    tmp << 0,0,-load;
    //! [Refinement]

    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);

    gsThinShellAssembler<3, real_t, true > shellAssemblerL(mp,basisL,bc,force,materialMatrix);
    gsThinShellAssemblerBase<real_t>*assemblerL = &shellAssemblerL;

    gsThinShellAssembler<3, real_t, true > shellAssemblerH(mp,basisH,bc,force,materialMatrix);
    gsThinShellAssemblerBase<real_t>*assemblerH = &shellAssemblerH;

    gsThinShellAssemblerDWRBase<real_t> * DWR;
    if (goal==1)
        DWR = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::DisplacementNorm>(mp,basisL,basisH,bc,force,materialMatrix);
    else if (goal==2)
        DWR = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Displacement>(mp,basisL,basisH,bc,force,materialMatrix);
    else if (goal==3)
        DWR = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStrain>(mp,basisL,basisH,bc,force,materialMatrix);
    else if (goal==4)
        DWR = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStress>(mp,basisL,basisH,bc,force,materialMatrix);
    else if (goal==5)
        DWR = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneForce>(mp,basisL,basisH,bc,force,materialMatrix);
    else if (goal==6)
        DWR = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStress>(mp,basisL,basisH,bc,force,materialMatrix);
    else
        GISMO_ERROR("Goal function unknown");

    gsSparseSolver<>::LU solver;
    gsVector<> solVector, updateVector;
    gsMultiPatch<> primalL,dualL,dualH;

    gsMatrix<> points(2,0);
    // points.col(0).setConstant(0.25);
    // points.col(1).setConstant(0.50);
    // points.col(2).setConstant(0.75);

    gsInfo << "Assembling primal... "<< std::flush;
    DWR->assembleMatrixL();
    DWR->assemblePrimalL();
    gsInfo << "done\n";

    gsInfo << "Solving primal, size ="<<DWR->matrixL().rows()<<","<<DWR->matrixL().cols()<<"... "<< "\n";
    solver.compute(DWR->matrixL());
    solVector = solver.solve(DWR->primalL());
    DWR->constructMultiPatchL(solVector,primalL);
    DWR->constructSolutionL(solVector,mp_def);
    index_t itMax = 100;
    real_t tol = 1e-14;
    real_t residual = DWR->primalL().norm();
    real_t residual0 = residual;
    real_t residualOld = residual;
    for (index_t it = 0; it != itMax; ++it)
    {
        DWR->assembleMatrixL(mp_def);
        DWR->assemblePrimalL(mp_def);
        solver.compute(DWR->matrixL());
        updateVector = solver.solve(DWR->primalL());
        solVector += updateVector;

        residual = DWR->primalL().norm();

        gsInfo<<"Iteration: "<< it
               <<", residue: "<< residual
               <<", update norm: "<<updateVector.norm()
               <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
               <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
               <<"\n";

        residualOld = residual;
        DWR->constructSolutionL(solVector,mp_def);
        if (updateVector.norm() < tol)
            break;
    }

    DWR->constructMultiPatchL(solVector,primalL);

    gsInfo << "done.\n";

    gsInfo << "Assembling dual vector (L)... "<< std::flush;
    DWR->assembleDualL(primalL);
    DWR->assembleDualL(points,primalL);
    gsInfo << "done.\n";

    gsInfo << "Solving dual (low), size = "<<DWR->matrixL().rows()<<","<<DWR->matrixL().cols()<<"... "<< std::flush;
    solVector = solver.solve(DWR->dualL());

    gsDebugVar(DWR->dualL());
    gsDebugVar(solVector);

    DWR->constructMultiPatchL(solVector,dualL);
    gsInfo << "done.\n";
    // gsInfo << "done." << " --> ";
    // gsInfo <<"Dual L error: \t"<<evL.integral(((dual_exL - zL_sol).norm()*meas(mapL)))<<"\n";

    gsInfo << "Assembling dual matrix (H)... "<< std::flush;
    DWR->assembleMatrixH(mp_def);
    gsInfo << "done.\n";

    gsInfo << "Assembling dual vector (H)... "<< std::flush;
    DWR->assembleDualH(primalL);
    DWR->assembleDualH(points,primalL);
    gsInfo << "done.\n";

    gsInfo << "Solving dual (high), size = "<<DWR->matrixH().rows()<<","<<DWR->matrixH().cols()<<"... "<< std::flush;
    solver.compute(DWR->matrixH());
    solVector = solver.solve(DWR->dualH());
    DWR->constructMultiPatchH(solVector,dualH);
    gsInfo << "done.\n";

    gsDebugVar(DWR->dualH());
    gsDebugVar(solVector);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real_t approx, exact = 0;

    gsThinShellAssemblerDWRBase<real_t> * DWR2;
    if (goal==1)
        DWR2 = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::DisplacementNorm>(mp,basisR,basisR,bc,force,materialMatrix);
    else if (goal==2)
        DWR2 = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::Displacement>(mp,basisR,basisR,bc,force,materialMatrix);
    else if (goal==3)
        DWR2 = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStrain>(mp,basisR,basisR,bc,force,materialMatrix);
    else if (goal==4)
        DWR2 = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneStress>(mp,basisR,basisR,bc,force,materialMatrix);
    else if (goal==5)
        DWR2 = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembraneForce>(mp,basisR,basisR,bc,force,materialMatrix);
    else if (goal==6)
        DWR2 = new gsThinShellAssemblerDWR<3,real_t,true,GoalFunction::MembranePStress>(mp,basisR,basisR,bc,force,materialMatrix);
    else
        GISMO_ERROR("Goal function unknown");

    exact += DWR2->computeGoal(mp_ex);
    exact += DWR2->computeGoal(points,mp_ex);
    exact -= DWR->computeGoal(mp_def);
    exact -= DWR->computeGoal(points,mp_def);
    approx = DWR->computeError(dualL,dualH,mp_def);

    gsInfo<<"approx = "<<approx<<"\n";
    gsInfo<<"Exact = "<<exact<<"\n";
    gsInfo<<"Efficiency = "<<approx/exact<<"\n";

    if (plot)
    {
        gsField<> fieldDL(mp, dualL);
        gsField<> fieldDH(mp, dualH);

        gsField<> fieldPL(mp, primalL);


        gsWriteParaview<>( fieldDL, "dualL", 1000);
        gsWriteParaview<>( fieldDH, "dualH", 1000);
        gsWriteParaview<>( fieldPL, "primalL", 1000);
    }


    return EXIT_SUCCESS;

}// end main
