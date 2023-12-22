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

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsThinShellAssemblerDWR.h>
#include <gsKLShell/src/gsThinShellUtils.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>

#include <gsAssembler/gsAdaptiveMeshing.h>
#include <gsAssembler/gsAdaptiveMeshingUtils.h>


//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Flag for refinemet criterion
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    index_t markstrat = 2;
    // Parameter for computing adaptive refinement threshold
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    real_t adaptRefParam = 0.9;


    //! [Parse command line]
    bool plot = false;
    bool write = false;
    index_t numRefine  = 1;
    index_t numRefineIni = 0;
    index_t numElevate = 1;
    index_t goal = 1;
    index_t component = 9;
    bool nonlinear = false;
    bool last = false;
    std::string fn;

    real_t E_modulus = 1;
    real_t PoissonRatio = 0.3;
    real_t thickness = 0.01;

    // real_t E_modulus = 1.0;
    // real_t PoissonRatio = 0.3;
    // real_t thickness = 1.0;

    int testCase = 0;

    int adaptivity = 0;

    std::string mesherOptionsFile("options/shell_mesher_options.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt("R", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", numRefineIni);
    cmd.addInt( "r", "refine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt("t", "testcase",
                "Test case: 0: Beam - pinned-pinned, 1: Beam - fixed-fixed, 2: beam - fixed-free, 3: plate - fully pinned, 4: plate - fully fixed, 5: circle - fully pinned, 6: 5: circle - fully fixed",
               testCase);

    cmd.addInt("A", "adaptivity", "Adaptivity scheme: 0) uniform refinement, 1) adaptive refinement, 2) adaptive refinement and coarsening", adaptivity);

    cmd.addInt( "g", "goal", "Goal function to use", goal );
    cmd.addInt( "c", "comp", "Component", component );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write convergence to file", write);
    cmd.addSwitch("last", "Only last refinement", last);
    cmd.addString( "O", "mesherOpt", "Input XML file for mesher options", mesherOptionsFile );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    real_t L = 1;
    real_t B = 1;
    if (testCase==0 || testCase==1)
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= L;
        mp.patch(0).coefs().col(1) *= B;
        mp.embed(3);
        mp.addAutoBoundaries();

        E_modulus = 1;
        thickness = 0.01;
        PoissonRatio = 0.3;

    }
    else if (testCase==2)
    {
        L=3;
        B=1;
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= L;
        mp.patch(0).coefs().col(1) *= B;
        mp.patch(0).coefs().col(1).array() -= B/2;
        mp.patch(0).coefs()(0,0) = B;
        mp.embed(3);
        mp.addAutoBoundaries();

        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0.3;
    }
    else
        GISMO_ERROR("Testcase "<<testCase<<" unknown");

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    gsMultiPatch<> mp_ex;
    gsMultiBasis<> basisR;
    if (testCase==0)
    {
        mp_ex = mp;
        mp_ex.degreeElevate(2);
        for (index_t r =0; r < std::max(std::min(numRefine,5),2); ++r)
            mp_ex.uniformRefine();
        basisR = gsMultiBasis<>(mp_ex);
    }

    // h-refine
    if (last)
    {
        for (index_t r =0; r < numRefine; ++r)
            mp.uniformRefine();
        numRefine = 0;
    }

    // Cast all patches of the mp object to THB splines
    if (adaptivity!=0)
    {
        gsTHBSpline<2,real_t> thb;
        for (index_t k=0; k!=mp.nPatches(); ++k)
        {
            gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
            thb = gsTHBSpline<2,real_t>(*geo);
            gsMatrix<> bbox = geo->support();
            for (index_t i = 0; i< numRefineIni; ++i)
                thb.refineElements(thb.basis().asElements(bbox));

            mp.patch(k) = thb;
        }
    }
    else
    {
        for (index_t i = 0; i< numRefineIni; ++i)
            mp.uniformRefine();
    }

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    std::string fx, fy, fz;
    fx = fy = fz = "0";

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    // real_t load = 1e-5;
    // gsFunctionExpr<> exact( "x","y","w:= 0; for (u := 1; u < 100; u += 2) { for (v := 1; v < 100; v += 2) { w += -16.0 * " + std::to_string(load) + " / ( pi^6*" + std::to_string(D) + " ) * 1 / (v * u * ( v^2 + u^2 )^2 ) * sin( v * pi * x) * sin(u * pi * y) } }",2);

    real_t ampl;

    gsMatrix<> points;
    typename std::vector<patchSide> bContainer;
    gsFunctionExpr<> neuData;
    if (testCase == 0)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north,condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::south,condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
        }

        bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 0 - x

        ampl = 1e2;
        char buffer[2000];
        // sprintf(buffer,"-%e^3*%e*%e*pi^4*sin(pi*x)*sin(pi*y)/(3*%e^2 - 3)",thickness,E_modulus,ampl,PoissonRatio);
        // // sprintf(buffer,"-2*%e^3*%e*%e/(3*%e^2 - 3)",thickness,E_modulus,ampl,PoissonRatio);
        sprintf(buffer,"-6*%e*%e^3*%e*(x^4 - 2*x^3 + 12*(-1/2 + y)^2*x^2 + (-12*y^2 + 12*y - 2)*x + y^4 - 2*y^3 + 3*y^2 - 2*y + 1/3)/(3*%e^2 - 3)",ampl,thickness,E_modulus,PoissonRatio);
        fz = buffer;

        points.resize(2,0);
        // bContainer.push_back(patchSide(0,4));
    }
    else if (testCase == 1)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north,condition_type::dirichlet, 0, i );
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, i );
            bc.addCondition(boundary::south,condition_type::dirichlet, 0, i );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, i );
        }

        gsVector<> point(2); point<< 0.5, 0.5 ;
        gsVector<> load (3); load << 0.0, 0.0, -1e-6 ;
        pLoads.addLoad(point, load, 0 );

        points.resize(2,1);
        points.col(0) = point;
    }
    else if (testCase == 2)
    {
        std::string nx = "0";
        std::string ny = "0";
        std::string nz = "0";

        real_t sigmax = 1e-1;
        char buffer[2000];
        sprintf(buffer,"%e ( 1 - y/%e)",sigmax,B/2);
        sprintf(buffer,"%e ( - y/%e)",sigmax,B/2);
        nx = buffer;

        neuData = gsFunctionExpr<> (nx,ny,nz,3);

        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        GISMO_ERROR("This does not work, since the computeErrorElwise does not take into account the neumann conditions since we miss a function gsExprEvaluator::integralBdrBcElwise");
        // bc.addCondition(boundary::east, condition_type::neumann, &neuData );
        // gsField<> field(mp.patch(0),neuData);
        // gsWriteParaview(field,"field");
    }
    else
        GISMO_ERROR("Test case "<<testCase<<" unknown");

    gsSparseSolver<>::LU solver;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Find reference solution
    if (testCase == 0)
    {
        char buffer[2000];
        std::string ux = "0";
        std::string uy = "0";
        // sprintf(buffer,"%e*sin(pi*x)*sin(pi*y)",ampl);
        // // sprintf(buffer,"%e*x*(x - 1)*y*(y - 1)",ampl);
        sprintf(buffer,"%e*x^2*(x - 1)^2*y^2*(y - 1)^2",ampl);
        std::string uz = buffer;

        ///////////////////////////////////////////////////////////////////////////////////////////////
        gsFunctionExpr<> exact(ux,uy,uz,3);
        gsField<> u_ex(mp,exact);

        gsWriteParaview(u_ex,"exact");

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(basisR);

        geometryMap G   = A.getMap(mp);

        space u = A.getSpace(basisR, 3);
        u.setup(bc,dirichlet::interpolation,0);
        auto    function = A.getCoeff(exact, G);
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

    }
    else if (testCase == 1)
    {
        if (gsFileManager::fileExists("deformed_plate.xml"))
            gsReadFile<>("deformed_plate.xml",mp_ex);
        else
            mp_ex = mp;
        basisR = gsMultiBasis<>(mp_ex);
    }
    else if (testCase == 2)
    {
        if (gsFileManager::fileExists("corner_piece.xml"))
            gsReadFile<>("corner_piece.xml",mp_ex);
        else
            mp_ex = mp;
        basisR = gsMultiBasis<>(mp_ex);
    }
    else
        GISMO_ERROR("Test case "<<testCase<<" unknown");

    gsField<> ump_ex(mp,mp_ex);
    gsWriteParaview(ump_ex,"mp_exact");

    //! [Refinement]
    gsFunctionExpr<> force(fx,fy,fz,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);

    // Compute exact goal
    gsThinShellAssemblerDWR<3,real_t,true> DWR2(mp,basisR,basisR,bc,force,materialMatrix);
    if (goal==1)
        DWR2.setGoal(GoalFunction::Displacement,component);
    else if (goal==2)
        DWR2.setGoal(GoalFunction::Stretch,component);
    else if (goal==3)
        DWR2.setGoal(GoalFunction::MembraneStrain,component);
    else if (goal==4)
        DWR2.setGoal(GoalFunction::PStrain,component);
    else if (goal==5)
        DWR2.setGoal(GoalFunction::MembraneStress,component);
    else if (goal==6)
        DWR2.setGoal(GoalFunction::PStress,component);
    else if (goal==7)
        DWR2.setGoal(GoalFunction::MembraneForce,component);
    else if (goal==8)
        DWR2.setGoal(GoalFunction::FlexuralStrain,component);
    else if (goal==9)
        DWR2.setGoal(GoalFunction::FlexuralStress,component);
    else if (goal==10)
        DWR2.setGoal(GoalFunction::FlexuralMoment,component);
    else
        GISMO_ERROR("Goal function unknown");


    real_t exactGoal = 0;
    exactGoal += DWR2.computeGoal(mp_ex);
    exactGoal += DWR2.computeGoal(bContainer,mp_ex);
    exactGoal += DWR2.computeGoal(points,mp_ex);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<real_t> exacts(numRefine+1);
    std::vector<real_t> approxs(numRefine+1);
    std::vector<real_t> efficiencies(numRefine+1);
    std::vector<real_t> numGoal(numRefine+1);
    std::vector<real_t> estGoal(numRefine+1);
    std::vector<real_t> exGoal(numRefine+1);
    std::vector<real_t> DoFs(numRefine+1);
    std::vector<real_t> hmins(numRefine+1);

    gsVector<> solVector, updateVector;
    gsMultiPatch<> primalL,dualL,dualH;

    gsParaviewCollection collection("solution");
    gsParaviewCollection indicator("indicator");
    gsParaviewCollection errors("error_elem_ref");

    std::vector<real_t> elErrors;

    gsAdaptiveMeshing<real_t> mesher;
    if (adaptivity!=0)
    {
        gsFileData<> fd_mesher(mesherOptionsFile);
        gsOptionList mesherOpts;
        fd_mesher.getFirst<gsOptionList>(mesherOpts);

        mesher = gsAdaptiveMeshing<real_t>(mp);
        mesher.options() = mesherOpts;
        mesher.getOptions();
    }

    real_t cummImprovement, cummImprovementOld;
    for (index_t r=0; r!=numRefine+1; r++)
    {

        // -----------------------------------------------------------------------------------------
        // ----------------------------Prepare basis------------------------------------------------
        // -----------------------------------------------------------------------------------------

        // Set bases
        gsMultiBasis<> basisL(mp);
        gsMultiBasis<> basisH = basisL;
        basisH.degreeElevate(1);

        gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n";
        gsInfo<<"Basis Dual:   "<<basisH.basis(0)<<"\n";

        gsThinShellAssemblerDWRBase<real_t> * DWR = new gsThinShellAssemblerDWR<3,real_t,true>(mp,basisL,basisH,bc,force,materialMatrix);
        if (goal==1)
            DWR->setGoal(GoalFunction::Displacement,component);
        else if (goal==2)
            DWR->setGoal(GoalFunction::Stretch,component);
        else if (goal==3)
            DWR->setGoal(GoalFunction::MembraneStrain,component);
        else if (goal==4)
            DWR->setGoal(GoalFunction::PStrain,component);
        else if (goal==5)
            DWR->setGoal(GoalFunction::MembraneStress,component);
        else if (goal==6)
            DWR->setGoal(GoalFunction::PStress,component);
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

        DWR->setPointLoads(pLoads);

        gsInfo << "Assembling primal... "<< std::flush;
        DWR->assembleMatrixL();
        DWR->assemblePrimalL();
        gsInfo << "done\n";

        gsInfo << "Solving primal, size ="<<DWR->matrixL().rows()<<","<<DWR->matrixL().cols()<<"... "<< std::flush;
        solver.compute(DWR->matrixL());
        solVector = solver.solve(DWR->primalL());
        DWR->constructMultiPatchL(solVector,primalL);
        DWR->constructSolutionL(solVector,mp_def);
        if (testCase==2) gsWrite(mp_def,"corner");
        gsInfo << "done.\n";


        gsInfo << "Assembling dual vector (L)... "<< std::flush;
        gsVector<> rhsL(DWR->numDofsL());
        rhsL.setZero();
        DWR->assembleDualL(primalL);
        rhsL += DWR->dualL();
        DWR->assembleDualL(bContainer,primalL);
        rhsL += DWR->dualL();
        DWR->assembleDualL(points,primalL);
        rhsL += DWR->dualL();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (L), size = "<<DWR->matrixL().rows()<<","<<DWR->matrixL().cols()<<"... "<< std::flush;
        solVector = solver.solve(rhsL);
        DWR->constructMultiPatchL(solVector,dualL);
        gsInfo << "done.\n";

        gsInfo << "Assembling dual matrix (H)... "<< std::flush;
        DWR->assembleMatrixH();
        gsInfo << "done.\n";

        gsInfo << "Assembling dual vector (H)... "<< std::flush;
        gsVector<> rhsH(DWR->numDofsH());
        rhsH.setZero();
        DWR->assembleDualH(primalL);
        rhsH = DWR->dualH();
        DWR->assembleDualH(bContainer,primalL);
        rhsH += DWR->dualH();
        DWR->assembleDualH(points,primalL);
        rhsH += DWR->dualH();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (H), size = "<<DWR->matrixH().rows()<<","<<DWR->matrixH().cols()<<"... "<< std::flush;

        solver.compute(DWR->matrixH());
        solVector = solver.solve(rhsH);
        DWR->constructMultiPatchH(solVector,dualH);
        gsInfo << "done.\n";

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (plot)
        {
            gsField<> Def(mp,primalL, true);
            std::string fileName = "solution" + util::to_string(r);
            gsWriteParaview<>(Def, fileName, 5000, true);
            fileName = "solution" + util::to_string(r) + "0";
            collection.addTimestep(fileName,r,".vts");
            collection.addTimestep(fileName,r,"_mesh.vtp");
        }

        real_t numError;
        if (plot)
        {
            std::string fileName = "indicator" + util::to_string(r);
            numError = DWR->computeError(dualL,dualH,false,fileName,5000,false,false);
            fileName = "indicator" + util::to_string(r) + "0";
            indicator.addTimestep(fileName,r,".vts");
        }
        else
            numError = DWR->computeError(dualL,dualH);

        typename gsBasis<real_t>::domainIter domIt = basisL.basis(0).makeDomainIterator();
        real_t diam = domIt->getMinCellLength();
        for (; domIt->good(); domIt->next())
            diam = domIt->getMinCellLength() < diam ? domIt->getMinCellLength() : diam;

        exacts[r] = 0;
        numGoal[r] = DWR->computeGoal(mp_def)+DWR->computeGoal(bContainer,mp_def)+DWR->computeGoal(points,mp_def);
        exGoal[r] = exactGoal;
        // DoFs[r] = DWR->numDofsL();
        DoFs[r] = basisL.basis(0).numElements();
        // DoFs[r] = basisL.basis(0).size();

        hmins[r] = diam;

        exacts[r] += exactGoal;
        exacts[r] -= numGoal[r];
        approxs[r] = numError;

        estGoal[r] = numGoal[r]+approxs[r];

        efficiencies[r] = approxs[r]/exacts[r];

        if (adaptivity==0)
        {
            mp.uniformRefine();
        }
        else if (adaptivity > 0)
        {
            gsFileData<> fd_mesher(mesherOptionsFile);
            gsOptionList mesherOpts;
            fd_mesher.getFirst<gsOptionList>(mesherOpts);

            // gsDebugVar(basisL.basis(0).minDegree());
            // gsDebugVar(mesherOpts.getInt("Convergence_alpha"));
            // gsDebugVar(mesherOpts.getInt("Convergence_beta"));
            // gsDebugVar(math::pow(1/2.,mesherOpts.getInt("Convergence_alpha") * basisL.basis(0).minDegree() + mesherOpts.getInt("Convergence_beta")));
            // gsDebugVar(exacts[r]);

            // gsDebugVar(exacts[r] * math::pow(1/2.,mesherOpts.getInt("Convergence_alpha") * basisL.basis(0).minDegree() + mesherOpts.getInt("Convergence_beta")));


            // gsDebugVar(approxs[r]);
            // gsDebugVar(exacts[r]);
            // gsDebugVar(numGoal[r]);
            elErrors = DWR->computeErrorElements(dualL, dualH,true);
            real_t error = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
            gsDebugVar(error);
            // gsDebugVar(*std::min_element(elErrors.begin(),elErrors.end()));
            // for (std::vector<real_t>::iterator it = elErrors.begin(); it != elErrors.end(); it++)
            // {
            //     *it = std::abs(*it);
            // }

            real_t abserror = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
            gsDebugVar(abserror);


            // Make container of the boxes
            gsHBoxContainer<2,real_t> markRef, markCrs;

            gsHBoxContainer<2,real_t> elts;
            index_t c = 0;
            real_t newError = 0;
            real_t oldError = 0;
            for (index_t p=0; p < mp.nPieces(); ++p)
            {
            // index_t p=0;
                typename gsBasis<real_t>::domainIter domIt = mp.basis(p).makeDomainIterator();
                gsHDomainIterator<real_t,2> * domHIt = nullptr;
                domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
                GISMO_ENSURE(domHIt!=nullptr,"Domain not loaded");
                for (; domHIt->good(); domHIt->next())
                {
                    gsHBox<2,real_t> box(domHIt,p);
                    box.setAndProjectError(elErrors[c],mesherOpts.getInt("Convergence_alpha"),mesherOpts.getInt("Convergence_beta"));
                    elts.add(box);
                    markRef.add(box);
                    oldError += box.error();
                    newError += box.projectedErrorRef();
                    c++;
                }
            }

            gsDebugVar(oldError);
            gsDebugVar(newError);


            gsElementErrorPlotter<real_t> err_eh(mp.basis(0),elErrors);
            const gsField<> elemError_eh( mp.patch(0), err_eh, true );
            gsWriteParaview<>( elemError_eh, "error_elem_ref" + util::to_string(r), 1000, true);
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,".vts");
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,"_mesh.vtp");


            mesher.markRef_into(elErrors,markRef);
            mesher.refine(markRef);
            // gsDebugVar(markRef);

            if (adaptivity>1)
            {
                mesher.markCrs_into(elErrors,markRef,markCrs);
                mesher.refine(markRef);
                mesher.unrefine(markCrs);
                // gsDebugVar(markCrs);
            }
            mesher.rebuild();
            mp_def = mp;
        }
        // else if (adaptivity > 0)
        // {
        //     gsFileData<> fd_mesher(mesherOptionsFile);
        //     gsOptionList mesherOpts;
        //     fd_mesher.getFirst<gsOptionList>(mesherOpts);

        //     // gsDebugVar(basisL.basis(0).minDegree());
        //     // gsDebugVar(mesherOpts.getInt("Convergence_alpha"));
        //     // gsDebugVar(mesherOpts.getInt("Convergence_beta"));
        //     // gsDebugVar(math::pow(1/2.,mesherOpts.getInt("Convergence_alpha") * basisL.basis(0).minDegree() + mesherOpts.getInt("Convergence_beta")));
        //     // gsDebugVar(exacts[r]);

        //     // gsDebugVar(exacts[r] * math::pow(1/2.,mesherOpts.getInt("Convergence_alpha") * basisL.basis(0).minDegree() + mesherOpts.getInt("Convergence_beta")));


        //     // gsDebugVar(approxs[r]);
        //     // gsDebugVar(exacts[r]);
        //     // gsDebugVar(numGoal[r]);
        //     elErrors = DWR->computeErrorElements(dualL, dualH);
        //     real_t error = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
        //     gsDebugVar(error);
        //     // gsDebugVar(*std::min_element(elErrors.begin(),elErrors.end()));
        //     // for (std::vector<real_t>::iterator it = elErrors.begin(); it != elErrors.end(); it++)
        //     // {
        //     //     *it = std::abs(*it);
        //     // }

        //     real_t abserror = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
        //     gsDebugVar(abserror);


        //     // Make container of the boxes
        //     gsHBoxContainer<2,real_t> markRef, markCrs;
        //     gsHBoxContainer<2,real_t> elts;
        //     index_t c = 0;
        //     real_t newError = 0;
        //     real_t oldError = 0;
        //     for (index_t p=0; p < mp.nPieces(); ++p)
        //     {
        //     // index_t p=0;
        //         typename gsBasis<real_t>::domainIter domIt = mp.basis(p).makeDomainIterator();
        //         gsHDomainIterator<real_t,2> * domHIt = nullptr;
        //         domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
        //         GISMO_ENSURE(domHIt!=nullptr,"Domain not loaded");
        //         for (; domHIt->good(); domHIt->next())
        //         {
        //             gsHBox<2,real_t> box(domHIt,p);
        //             box.setAndProjectError(elErrors[c],mesherOpts.getInt("Convergence_alpha"),mesherOpts.getInt("Convergence_beta"));
        //             elts.add(box);
        //             markCrs.add(box);
        //             oldError += box.error();
        //             newError += box.projectedErrorCrs();
        //             c++;
        //         }
        //     }

        //     gsDebugVar(oldError);
        //     gsDebugVar(newError);

        //     mesher.markCrs_into(elErrors,markCrs);
        //     mesher.unrefine(markCrs);
        //     mesher.rebuild();
        //     mp_def = mp;
        // }
        else if (adaptivity < 0)
        {
            gsFileData<> fd_mesher(mesherOptionsFile);
            gsOptionList mesherOpts;
            fd_mesher.getFirst<gsOptionList>(mesherOpts);

            elErrors = DWR->computeErrorElements(dualL, dualH);
            real_t error = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
            gsDebugVar(error);
            gsDebugVar(*std::min_element(elErrors.begin(),elErrors.end()));
            for (std::vector<real_t>::iterator it = elErrors.begin(); it != elErrors.end(); it++)
            {
                *it = std::abs(*it);
            }

            real_t abserror = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
            gsDebugVar(abserror);

            gsElementErrorPlotter<real_t> err_eh(mp.basis(0),elErrors);
            const gsField<> elemError_eh( mp.patch(0), err_eh, true );
            gsWriteParaview<>( elemError_eh, "error_elem_ref" + util::to_string(r), 1000, true);
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,".vts");
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,"_mesh.vtp");

            std::vector<bool> elMarked( elErrors.size() );
            gsHBoxContainer<2,real_t> markRef, markCrs;

            // Make container of the boxes
            gsHBoxContainer<2,real_t> elts;
            index_t c = 0;
            for (index_t p=0; p < mp.nPieces(); ++p)
            {
            // index_t p=0;
                typename gsBasis<real_t>::domainIter domIt = mp.basis(p).makeDomainIterator();
                gsHDomainIterator<real_t,2> * domHIt = nullptr;
                domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
                GISMO_ENSURE(domHIt!=nullptr,"Domain not loaded");
                for (; domHIt->good(); domHIt->next())
                {
                    gsHBox<2,real_t> box(domHIt,p);
                    box.setAndProjectError(elErrors[c],mesherOpts.getInt("Convergence_alpha"),mesherOpts.getInt("Convergence_beta"));
                    elts.add(box);
                    c++;
                }
            }
            typename gsHBoxContainer<2,real_t>::Container container = elts.toContainer();
            typename gsHBoxContainer<2,real_t>::SortedContainer scontainer = gsHBoxUtils<2,real_t>::Sort(container);

            gsWriteParaview(elts,"elts");

            // --------------- adaptive refinement ---------------
            cummImprovement = 0;
            index_t index = std::rand()%container.size();
            // gsDebugVar(index);
            gsHBox<2,real_t> refbox = *std::next(container.begin(),index);
            typename gsHBoxContainer<2,real_t>::HContainer hcontainer = gsHBoxUtils<2,real_t>::markAdmissible(refbox,2);
            gsHBoxContainer<2,real_t>refined(hcontainer);

            gsDebugVar(refined);

            for (typename gsHBoxContainer<2,real_t>::HIterator hit = refined.begin(); hit!=refined.end(); hit++)
                for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
                {
                    auto pred = [it](auto b){return it->isSame(b);};
                    typename gsHBox<2,real_t>::SortedContainer::iterator found = std::find_if(scontainer.begin(),scontainer.end(),pred);
                    if (found!=scontainer.end())
                        cummImprovement += found->error()-found->projectedErrorRef();
                    else
                        GISMO_ERROR("Hi");
                }

            mesher.refine(refined);
            gsWriteParaview(refined,"refined");

            mesher.rebuild();
            mp_def = mp;

            gsDebugVar(cummImprovement);
            if (r>0) gsDebugVar(approxs[r]-approxs[r-1]);
        }

        delete DWR;

    }

    if (plot)
    {
        collection.save();
        indicator.save();
        errors.save();
    }

    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    gsInfo<<"Ref.\tApprox    \tExact     \tEfficiency\tNumGoal   \tEstGoal   \texGoal    \t#elements \thmin \tdApprox      \n";
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    for(index_t r=0; r!=numRefine+1; r++)
    {
        gsInfo  <<std::setw(4 )<<std::left<<r<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<approxs[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<exacts[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<efficiencies[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<numGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<estGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<exGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<DoFs[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<hmins[r];
        if (r>0)
            gsInfo  <<std::setw(10)<<std::left<<approxs[r-1]-approxs[r];
        gsInfo  <<"\n";
    }
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";

    if (write)
    {
        std::string filename;
        filename = "example_shell3D_DWR_r" + std::to_string(numRefine) + "_e" + std::to_string(numElevate) + "_g" + std::to_string(goal) + "_C" + std::to_string(component);
        filename = filename + ".csv";
        std::ofstream file_out;
        file_out.open (filename);

        file_out<<"Ref,Approx,Exact,Efficiency,NumGoal,EstGoal,exGoal,DoFs,hmin\n";
        for(index_t r=0; r!=numRefine+1; r++)
        {
            file_out<<r<<","<<approxs[r]<<","<<exacts[r]<<","<<efficiencies[r]<<","<<numGoal[r]<<","<<estGoal[r]<<","<<exGoal[r]<<","<<DoFs[r]<<","<<hmins[r]<<"\n";
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


    delete materialMatrix;
    return EXIT_SUCCESS;

}// end main
