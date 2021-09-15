/** @file example_shell3D.cpp

    @brief Examples for the surface thin shells including the shell obstacle course

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <iostream>
#include <fstream>
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
    bool write = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t goal = 1;
    index_t component = 9;
    bool nonlinear = false;
    bool loop = false;
    bool adaptive = false;
    std::string fn;

    // real_t E_modulus = 1;
    // real_t PoissonRatio = 0.3;
    // real_t thickness = 0.01;

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
    cmd.addInt( "g", "goal", "Goal function to use", goal );
    cmd.addInt( "C", "comp", "Component", component );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write convergence to file", write);
    cmd.addSwitch("loop", "Uniform Refinement loop", loop);
    cmd.addSwitch("adaptive", "Adaptive Refinemenct loop", adaptive);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsVector<> pts(2);
    pts.setConstant(0.25);

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    // Unit square
    // real_t L = 2;
    real_t L = 1;
    real_t B = 1;
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.patch(0).coefs().col(0) *= L;
    mp.patch(0).coefs().col(1) *= B;
    mp.embed(3);
    mp.addAutoBoundaries();

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    if (!loop)
    {
        for (index_t r =0; r < numRefine; ++r)
            mp.uniformRefine();
        numRefine = 0;
    }

    // Cast all patches of the mp object to THB splines
    if (adaptive)
    {
        gsTHBSpline<2,real_t> thb;
        for (index_t k=0; k!=mp.nPatches(); ++k)
        {
            gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
            thb = gsTHBSpline<2,real_t>(*geo);
            mp.patch(k) = thb;
        }
    }

    gsMultiPatch<> mp_ex = mp;
    mp_ex.degreeElevate(2);
    mp_ex.uniformRefine();
    gsMultiBasis<> basisR(mp_ex);

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    // real_t load = 1e-5;
    real_t load = 1.0;

    // real_t D = E_modulus * math::pow(thickness,3) / ( 12 * ( 1- math::pow(PoissonRatio,2) ) );
    // gsFunctionExpr<> exact( "x","y","w:= 0; for (u := 1; u < 100; u += 2) { for (v := 1; v < 100; v += 2) { w += -16.0 * " + std::to_string(load) + " / ( pi^6*" + std::to_string(D) + " ) * 1 / (v * u * ( v^2 + u^2 )^2 ) * sin( v * pi * x) * sin(u * pi * y) } }",2);

    for (index_t i=0; i!=3; ++i)
    {
        bc.addCondition(boundary::north,condition_type::dirichlet, 0, i );
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, i );
        bc.addCondition(boundary::south,condition_type::dirichlet, 0, i );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, i );
    }

    bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x
    bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x
    bc.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x
    bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x

    real_t a = 1e0;
    real_t nu = PoissonRatio;
    real_t E = E_modulus;
    real_t t = thickness;

    char buffer[2000];
    sprintf(buffer,"%e^2 / (%e^2-1) * 8*%e*(y - 1)^2*x*(x - 1)*(-1/2 + x)*%e*y^2*(((%e + 7)*y^2 + (-%e - 7)*y + 3/2 + %e/2)*x^4 + ((-2*%e - 14)*y^2 + (2*%e + 14)*y - %e - 3)*x^3 + (6*y^4 - 12*y^3 + (%e + 13)*y^2 + (-%e - 7)*y + 3/2 + %e/2)*x^2 - 6*y^2*(y - 1)^2*x + y^2*(y - 1)^2)*%e^2",
        a,nu,E,t,nu,nu,nu,nu,nu,nu,nu,nu,nu,a);
    std::string fx = buffer;


    sprintf(buffer,"%e^2 / (%e^2-1) * 48*%e*(y - 1)*x^2*(x - 1)^2*%e*(-1/2 + y)*y*(((%e/6 + 7/6)*x^2 + (-%e/6 - 7/6)*x + %e/12 + 1/4)*y^4 + ((-%e/3 - 7/3)*x^2 + (%e/3 + 7/3)*x - %e/6 - 1/2)*y^3 + (x^4 - 2*x^3 + (%e/6 + 13/6)*x^2 + (-%e/6 - 7/6)*x + %e/12 + 1/4)*y^2 - x^2*(x - 1)^2*y + x^2*(x - 1)^2/6)*%e^2",
    a,nu,E,t,nu,nu,nu,nu,nu,nu,nu,nu,nu,a);
    std::string fy = buffer;

    // sprintf(buffer,"%e / (3*(%e^2-1)) * 864*%e*%e*((y - 1)^2*y^2*(-1/2 + y)^2*(y^2 - y + 1/6)*%e^2*x^12 - 6*(y - 1)^2*y^2*(-1/2 + y)^2*(y^2 - y + 1/6)*%e^2*x^11 + (22*(y - 1)^2*y^2*(y^6 - 3*y^5 + 75/8*y^4 - 55/4*y^3 + 393/44*y^2 - 225/88*y + 45/176)*%e^2*x^10)/9 - (110*(y - 1)^2*y^2*(y^6 - 3*y^5 + 39/8*y^4 - 19/4*y^3 + 225/88*y^2 - 15/22*y + 3/44)*%e^2*x^9)/9 + (y - 1)^2*(y^8 - 4*y^7 + 1117/36*y^6 - 949/12*y^5 + 1747/18*y^4 - 2411/36*y^3 + 245/9*y^2 - 25/4*y + 5/8)*y^2*%e^2*x^8 - 4*(y - 1)^2*(y^8 - 4*y^7 + 457/36*y^6 - 289/12*y^5 + 1741/72*y^4 - 116/9*y^3 + 67/18*y^2 - 5/8*y + 1/16)*y^2*%e^2*x^7 + (77*(y - 1)^2*(y^8 - 4*y^7 + 1952/231*y^6 - 874/77*y^5 + 100/11*y^4 - 908/231*y^3 + 62/77*y^2 - 5/77*y + 1/154)*y^2*%e^2*x^6)/12 - (21*(y - 1)^4*(y^4 - 2*y^3 + 361/189*y^2 - 172/189*y + 41/189)*y^4*%e^2*x^5)/4 + (-1/144*%e^2 + 55/24*%e^2*y^12 - 55/4*%e^2*y^11 + 839/24*%e^2*y^10 - 195/4*%e^2*y^9 + 2905/72*%e^2*y^8 - 725/36*%e^2*y^7 + 145/24*%e^2*y^6 - 41/36*%e^2*y^5 + 5/36*%e^2*y^4)*x^4 + (1/72*%e^2 - 1/2*%e^2*y^12 + 3*%e^2*y^11 - 15/2*%e^2*y^10 + 10*%e^2*y^9 - 15/2*%e^2*y^8 + 3*%e^2*y^7 - 1/2*%e^2*y^6)*x^3 + (-1/48*%e^2 + 1/24*%e^2*y^12 - 1/4*%e^2*y^11 + 5/8*%e^2*y^10 - 5/6*%e^2*y^9 + 5/8*%e^2*y^8 - 1/4*%e^2*y^7 + 1/24*%e^2*y^6 - 1/12*%e^2*y^2 + 1/12*%e^2*y)*x^2 + %e^2*(y^2 - y + 1/6)*x/12 - %e^2*(y^4 - 2*y^3 + 3*y^2 - 2*y + 1/3)/144)*%e"
    // ,a,nu,E,t,a,a,a,a,a,a,a,a,t,a,a,a,a,a,a,a,a,a,t,a,a,a,a,a,a,a,t,a,a,a,a,a,a,a,t,t,t,t,a);
    sprintf(buffer,"-6*%e*%e^3*%e*(x^4 - 2*x^3 + 12*(-1/2 + y)^2*x^2 + (-12*y^2 + 12*y - 2)*x + y^4 - 2*y^3 + 3*y^2 - 2*y + 1/3)/(3*%e^2 - 3)",a,t,E,nu);
    std::string fz = buffer;

    std::string ux = "0";
    std::string uy = "0";
    // sprintf(buffer,"%e*sin(pi*x)*sin(pi*y)",a);
    // // sprintf(buffer,"%e*x*(x - 1)*y*(y - 1)",a);
    sprintf(buffer,"%e*x^2*(x - 1)^2*y^2*(y - 1)^2",a);
    std::string uz = buffer;

    // tmp << 0,0,-load;
    //! [Refinement]

    // gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> force(fx,fy,fz,3);
    gsFunctionExpr<> thick(std::to_string(thickness), 3);
    gsFunctionExpr<> Emod(std::to_string(E_modulus),3);
    gsFunctionExpr<> Pois(std::to_string(PoissonRatio),3);

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &Emod;
    parameters[1] = &Pois;
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,thick,parameters,options);

    gsSparseSolver<>::LU solver;

    gsMatrix<> points(2,0);
    // points.col(0).setConstant(0.25);
    // points.col(1).setConstant(0.50);
    // points.col(2).setConstant(0.75);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    gsFunctionExpr<> exact(ux,uy,uz,3);
    gsField<> u_ex(mp,exact);

    gsWriteParaview(u_ex,"exact");

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(basisR);
    geometryMap G   = A.getMap(mp);

    space u = A.getSpace(basisR, 3);
    u.setup(bc,dirichlet::interpolation,0);
    auto function = A.getCoeff(exact, G);
    A.initSystem();
    A.assemble(u*u.tr()*meas(G),u * function*meas(G));
    solver.compute(A.matrix());
    gsMatrix<> result = solver.solve(A.rhs());

    solution u_sol = A.getSolution(u, result);

    gsMatrix<> cc;
    for ( size_t k =0; k!=mp_ex.nPatches(); ++k) // Deform the geometry
    {
        // // extract deformed geometry
        u_sol.extract(cc, k);
        mp_ex.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }

    gsExprEvaluator<> ev(A);
    gsDebugVar(ev.integral((u_sol-function).sqNorm()));


    // gsThinShellAssembler<3,real_t,true> assembler(mp,basisR,bc,force,materialMatrix);
    // gsMultiPatch<> mp_ex2 = mp;
    // assembler.projectL2_into(exact,mp_ex2);

    // gsMatrix<> coefs;
    // gsQuasiInterpolate<real_t>::localIntpl(basisR.basis(0), exact, coefs);
    // gsMultiPatch<> mp_ex;
    // mp_ex.addPatch(*basisR.basis(0).makeGeometry(give(coefs)));

    gsField<> ump_ex(mp,mp_ex);
    gsWriteParaview(ump_ex,"mp_exact");

    gsThinShellAssemblerDWRBase<real_t> * DWR2;

    // gsMultiPatch<> mp_ex;
    // // gsReadFile<>("deformations/deformed_plate_r7e5.xml",mp_ex);
    // gsReadFile<>("deformed_plate_lin_T=" + std::to_string(thickness) + ".xml",mp_ex);
    // gsMultiBasis<> basisR(mp_ex);

    DWR2 = new gsThinShellAssemblerDWR<3,real_t,true>(mp,basisR,basisR,bc,force,materialMatrix);
    if (goal==1)
        DWR2->setGoal(GoalFunction::Displacement,component);
    else if (goal==2)
        DWR2->setGoal(GoalFunction::Stretch,component);
    else if (goal==3)
        DWR2->setGoal(GoalFunction::MembraneStrain,component);
    else if (goal==4)
        DWR2->setGoal(GoalFunction::MembranePStrain,component);
    else if (goal==5)
        DWR2->setGoal(GoalFunction::MembraneStress,component);
    else if (goal==6)
        DWR2->setGoal(GoalFunction::MembranePStress,component);
    else if (goal==7)
        DWR2->setGoal(GoalFunction::MembraneForce,component);
    else if (goal==8)
        DWR2->setGoal(GoalFunction::FlexuralStrain,component);
    else if (goal==9)
        DWR2->setGoal(GoalFunction::FlexuralStress,component);
    else if (goal==10)
        DWR2->setGoal(GoalFunction::FlexuralMoment,component);
    else
        GISMO_ERROR("Goal function unknown");

    real_t exactGoal = 0;
    exactGoal += DWR2->computeGoal(mp_ex);
    exactGoal += DWR2->computeGoal(points,mp_ex);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<real_t> exacts(numRefine+1);
    std::vector<real_t> approxs(numRefine+1);
    std::vector<real_t> efficiencies(numRefine+1);
    std::vector<real_t> numGoal(numRefine+1);
    std::vector<real_t> estGoal(numRefine+1);
    std::vector<real_t> exGoal(numRefine+1);

    gsVector<> solVector, updateVector;
    gsMultiPatch<> primalL,dualL,dualH;

    gsThinShellAssemblerDWRBase<real_t> * DWR;

    gsParaviewCollection collection("solution");
    for (index_t r=0; r!=numRefine+1; r++)
    {
        if (adaptive)
        {
            // [Mark elements for refinement]
            std::vector<index_t> bools(mp.basis(0).numElements());
            std::vector<bool> refVec(mp.basis(0).numElements());

            //////// RANDOMLY
            std::srand(std::time(nullptr)); // use current time as seed for random generator
            std::generate(bools.begin(), bools.end(), rand);
            for (index_t k = 0; k!=bools.size(); k++)
                refVec[k] = static_cast<bool>(std::round(static_cast<real_t>(bools[k]) / ( RAND_MAX+1u )));

            gsRefineMarkedElements(mp,refVec,0);
            gsMultiPatch<> mp_def = mp;
        }

        // -----------------------------------------------------------------------------------------
        // ----------------------------Prepare basis------------------------------------------------
        // -----------------------------------------------------------------------------------------
        // Set deformed multipatch
        mp_def = mp;

        // Set bases
        gsMultiBasis<> basisL(mp);
        gsMultiBasis<> basisH = basisL;
        basisH.degreeElevate(1);

        gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n";
        gsInfo<<"Basis Dual:   "<<basisH.basis(0)<<"\n";

        DWR = new gsThinShellAssemblerDWR<3,real_t,true>(mp,basisL,basisH,bc,force,materialMatrix);
        if (goal==1)
            DWR->setGoal(GoalFunction::Displacement,component);
        else if (goal==2)
            DWR->setGoal(GoalFunction::Stretch,component);
        else if (goal==3)
            DWR->setGoal(GoalFunction::MembraneStrain,component);
        else if (goal==4)
            DWR->setGoal(GoalFunction::MembranePStrain,component);
        else if (goal==5)
            DWR->setGoal(GoalFunction::MembraneStress,component);
        else if (goal==6)
            DWR->setGoal(GoalFunction::MembranePStress,component);
        else if (goal==7)
            DWR->setGoal(GoalFunction::MembraneForce,component);
        else if (goal==8)
            DWR->setGoal(GoalFunction::FlexuralStrain,component);
        else if (goal==9)
            DWR->setGoal(GoalFunction::FlexuralStress,component);
        else if (goal==10)
            DWR->setGoal(GoalFunction::FlexuralMoment,component);
        else
            GISMO_ERROR("Goal function unknown");

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
        gsVector<> rhsL;
        DWR->assembleDualL(primalL);
        rhsL = DWR->dualL();
        DWR->assembleDualL(points,primalL);
        rhsL += DWR->dualL();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (L), size = "<<DWR->matrixL().rows()<<","<<DWR->matrixL().cols()<<"... "<< std::flush;
        solVector = solver.solve(rhsL);

        DWR->constructMultiPatchL(solVector,dualL);
        gsInfo << "done.\n";
        // gsInfo << "done." << " --> ";
        // gsInfo <<"Dual L error: \t"<<evL.integral(((dual_exL - zL_sol).norm()*meas(mapL)))<<"\n";

        gsInfo << "Assembling dual matrix (H)... "<< std::flush;
        DWR->assembleMatrixH(mp_def);
        gsInfo << "done.\n";

        gsInfo << "Assembling dual vector (H)... "<< std::flush;
        gsVector<> rhsH;
        DWR->assembleDualH(primalL);
        rhsH = DWR->dualH();
        DWR->assembleDualH(points,primalL);
        rhsH += DWR->dualH();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (H), size = "<<DWR->matrixH().rows()<<","<<DWR->matrixH().cols()<<"... "<< std::flush;
        solver.compute(DWR->matrixH());
        solVector = solver.solve(rhsH);
        DWR->constructMultiPatchH(solVector,dualH);
        gsInfo << "done.\n";

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        exacts[r] = 0;
        numGoal[r] = DWR->computeGoal(mp_def)+DWR->computeGoal(points,mp_def);
        exGoal[r] = exactGoal;

        exacts[r] += exactGoal;
        exacts[r] -= numGoal[r];
        approxs[r] = DWR->computeError(dualL,dualH,mp_def);

        estGoal[r] = numGoal[r]+approxs[r];

        efficiencies[r] = approxs[r]/exacts[r];

        if (!adaptive)
            mp.uniformRefine();
    }

    if (plot) collection.save();

    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    gsInfo<<"Ref.\tApprox    \tExact     \tEfficiency\tNumGoal   \tEstGoal   \texGoal    \n";
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    for(index_t r=0; r!=numRefine+1; r++)
    {
        gsInfo  <<std::setw(4 )<<std::left<<r<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<approxs[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<exacts[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<efficiencies[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<numGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<estGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<exGoal[r]<<"\n";
    }
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";


    if (write)
    {
        std::string filename;
        filename = "example_shell3D_DWR_NL_r" + std::to_string(numRefine) + "_e" + std::to_string(numElevate) + "_g" + std::to_string(goal) + "_C" + std::to_string(component);
        filename = filename + ".csv";
        std::ofstream file_out;
        file_out.open (filename);

        file_out<<"Ref,Approx,Exact,Efficiency,NumGoal,EstGoal,exGoal\n";
        for(index_t r=0; r!=numRefine+1; r++)
        {
            file_out<<r<<","<<approxs[r]<<","<<exacts[r]<<","<<efficiencies[r]<<","<<numGoal[r]<<","<<estGoal[r]<<","<<exGoal[r]<<"\n";
        }

        file_out.close();
    }

    if (plot)
    {
        gsField<> fieldDL(mp, dualL);
        gsField<> fieldDH(mp, dualH);

        gsField<> fieldPL(mp, primalL);


        gsWriteParaview<>( fieldDL, "dualL", 1000);
        gsWriteParaview<>( fieldDH, "dualH", 1000);
        gsWriteParaview<>( fieldPL, "primalL", 1000);
    }


    delete DWR;
    delete DWR2;
    delete materialMatrix;
    return EXIT_SUCCESS;

}// end main
