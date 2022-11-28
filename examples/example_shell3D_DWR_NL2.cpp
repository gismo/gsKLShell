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
#include <gsAssembler/gsAdaptiveMeshing.h>
#include <gsAssembler/gsAdaptiveMeshingUtils.h>

#include <gsStructuralAnalysis/gsALMLoadControl.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool write = false;
    index_t numRefine  = 1;
    index_t numRefineIni = 0;
    index_t numElevate = 1;
    index_t goal = 1;
    index_t component = 9;
    bool loop = false;
    std::string fn;

    // real_t E_modulus = 1;
    // real_t PoissonRatio = 0.3;
    // real_t thickness = 0.01;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t thickness = 1e-3;

    index_t steps = 10;

    int adaptivity = 0;

    std::string dirname  = "Static_Pointload";

    std::string mesherOptionsFile("options/shell_mesher_options.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt("R", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", numRefineIni);
    cmd.addInt("r", "refine", "Maximum number of adaptive refinement steps to perform",
               numRefine);

    cmd.addInt("N", "steps", "Number of ALM steps", steps);

    cmd.addInt("A", "adaptivity", "Adaptivity scheme: 0) uniform refinement, 1) adaptive refinement, 2) adaptive refinement and coarsening", adaptivity);

    cmd.addInt( "g", "goal", "Goal function to use", goal );
    cmd.addInt( "C", "comp", "Component", component );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write convergence to file", write);
    cmd.addSwitch("loop", "Uniform Refinement loop", loop);
    cmd.addString( "O", "mesherOpt", "Input XML file for mesher options", mesherOptionsFile );
    cmd.addString( "o", "output", "output directory", dirname );


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    // Unit square
    // real_t L = 2;
    real_t L = 0.5;
    real_t B = 0.5;
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
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    real_t load = 1e-7;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    bc.addCondition(boundary::north,condition_type::dirichlet,0,0,false,1);
    bc.addCondition(boundary::north,condition_type::clamped,0,0,false,0);
    bc.addCondition(boundary::north,condition_type::clamped,0,0,false,2);

    bc.addCondition(boundary::east,condition_type::dirichlet,0,0,false,0);
    bc.addCondition(boundary::east,condition_type::clamped,0,0,false,1);
    bc.addCondition(boundary::east,condition_type::clamped,0,0,false,2);

    bc.addCornerValue(boundary::southwest, 0.0, 0, 0, -1); // (corner,value, patch, unknown)

    gsVector<> pointvec(2);
    pointvec<< 1.0, 1.0 ;
    gsVector<> loadvec (3);
    loadvec << 0.0, 0.0, load ;
    pLoads.addLoad(pointvec, loadvec, 0 );

    // gsConstantFunction<> force(tmp,3);
    // gsFunctionExpr<> force("0","0","if (sqrt((x-0.5)^2+(y-0.5)^2)<0.1){-1e-4} else{0}",3);
    gsFunctionExpr<> force("0","0","0",3);
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

    gsSparseSolver<real_t>::uPtr solver;
#ifdef GISMO_WITH_PARDISO
    solver = gsSparseSolver<real_t>::get( "PardisoLU");
#else
    solver = gsSparseSolver<real_t>::get( "SimplicialLDLT");
#endif

    gsMatrix<> points(2,1);
    points.col(0).setConstant(0.5);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<real_t> exacts(numRefine+1);
    std::vector<real_t> approxs(numRefine+1);
    std::vector<real_t> efficiencies(numRefine+1);
    std::vector<real_t> numGoal(numRefine+1);
    std::vector<real_t> estGoal(numRefine+1);
    std::vector<real_t> exGoal(numRefine+1);
    std::vector<real_t> Uz(numRefine+1);
    std::vector<real_t> DoFs(numRefine+1);

    gsVector<> solVector;
    gsMultiPatch<> primalL,dualL,dualH;

    gsThinShellAssemblerDWRBase<real_t> * DWR;

    std::string commands = "mkdir -p " + dirname;
    const char *command  = commands.c_str();
    system(command);

    gsParaviewCollection collection(dirname + "/" + "solution");
    gsParaviewCollection errors(dirname + "/" + "error_elem_ref");

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

    for (index_t r=0; r!=numRefine+1; r++)
    {
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

        DWR->assemblePrimalL();
        gsVector<> Force = DWR->primalL();

        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
        // Function for the Jacobian
        Jacobian_t Jacobian = [&DWR,&mp_def](gsVector<real_t> const &x)
        {
          DWR->constructSolutionL(x,mp_def);
          DWR->assembleMatrixL(mp_def);
          gsSparseMatrix<real_t> m = DWR->matrixL();
          return m;
        };
        // Function for the Residual
        ALResidual_t ALResidual = [&DWR,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
        {
          DWR->constructSolutionL(x,mp_def);
          DWR->assemblePrimalL(mp_def);
          gsVector<real_t> Fint = -(DWR->primalL() - force);
          gsVector<real_t> result = Fint - lam * force;
          return result; // - lam * force;
        };

        gsInfo << "Solving primal, size ="<<DWR->matrixL().rows()<<","<<DWR->matrixL().cols()<<"... "<< "\n";
        real_t dL = 1.0/steps;
        gsALMLoadControl<real_t> arcLength(Jacobian, ALResidual, Force);
#ifdef GISMO_WITH_PARDISO
        arcLength.options().setString("Solver","PardisoLU"); // LDLT solver
#else
        arcLength.options().setString("Solver","SimplicialLDLT"); // LDLT solver
#endif
        arcLength.options().setReal("Length",dL);
        // arcLength.options().setReal("Tol",tol);
        // arcLength.options().setReal("TolU",tolU);
        // arcLength.options().setReal("TolF",tolF);
        arcLength.options().setInt("MaxIter",50);
        arcLength.options().setSwitch("Verbose",true);
        arcLength.options().setInt("BifurcationMethod",gsALMBase<real_t>::bifmethod::Nothing);
        gsInfo<<arcLength.options();
        arcLength.applyOptions();
        arcLength.initialize();

        real_t dL0 = dL;
        real_t Lold = 0;
        real_t L = 0;
        index_t k = 0;
        gsMatrix<> Uold(Force.rows(),1);
        Uold.setZero();
        while (L < 1 && std::abs(L-1)>1e-14)
        {
            gsInfo<<"Load step "<< k<<"\n";
            arcLength.setLength(dL);
            arcLength.step();

            if (!(arcLength.converged()))
            {
              gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
              dL /= 2.;
              arcLength.setLength(dL);
              arcLength.setSolution(Uold,Lold);
              continue;
            }
            dL = dL0;
            solVector = arcLength.solutionU();
            Uold = solVector;
            L = Lold = arcLength.solutionL();
            k++;
            dL0 = dL = std::min(1-L,dL0);
            gsInfo<<"dL = "<<dL<<"; 1-L = "<<1-L<<"\n";
            if (dL > 1-L) dL0 = dL = 1-L;
        }
        DWR->constructMultiPatchL(solVector,primalL);
        DWR->constructSolutionL(solVector,mp_def);
        gsInfo << "done.\n";

        gsInfo << "Assembling dual matrix (L)... "<< std::flush;
        DWR->assembleMatrixL(mp_def);
        solver->compute(DWR->matrixL());
        gsInfo << "Assembling dual vector (L)... "<< std::flush;
        gsVector<> rhsL(DWR->numDofsL());
        rhsL.setZero();
        // DWR->assembleDualL(primalL);
        // rhsL += DWR->dualL();
        DWR->assembleDualL(points,primalL);
        rhsL += DWR->dualL();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (L), size = "<<DWR->matrixL().rows()<<","<<DWR->matrixL().cols()<<"... "<< std::flush;
        solVector = solver->solve(rhsL);

        DWR->constructMultiPatchL(solVector,dualL);
        gsInfo << "done.\n";
        // gsInfo << "done." << " --> ";
        // gsInfo <<"Dual L error: \t"<<evL.integral(((dual_exL - zL_sol).norm()*meas(mapL)))<<"\n";

        gsInfo << "Assembling dual matrix (H)... "<< std::flush;
        DWR->assembleMatrixH(mp_def);
        gsInfo << "done.\n";

        gsInfo << "Assembling dual vector (H)... "<< std::flush;
        gsVector<> rhsH(DWR->numDofsH());
        rhsH.setZero();
        // DWR->assembleDualH(primalL);
        // rhsH += DWR->dualH();
        DWR->assembleDualH(points,primalL);
        rhsH += DWR->dualH();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (H), size = "<<DWR->matrixH().rows()<<","<<DWR->matrixH().cols()<<"... "<< std::flush;
        solver->compute(DWR->matrixH());
        solVector = solver->solve(rhsH);
        DWR->constructMultiPatchH(solVector,dualH);
        gsInfo << "done.\n";

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (plot)
        {
            gsField<> Def(mp_def,primalL, true);
            std::string fileName = dirname + "/" + "solution" + util::to_string(r);
            gsWriteParaview<>(Def, fileName, 5000, true);
            fileName = "solution" + util::to_string(r) + "0";
            collection.addTimestep(fileName,r,".vts");
            collection.addTimestep(fileName,r,"_mesh.vtp");
        }


        exacts[r] = 0;
        numGoal[r] = 0;
        // numGoal[r] += DWR->computeGoal(mp_def);
        numGoal[r] += DWR->computeGoal(points,mp_def);
        DoFs[r] = basisL.basis(0).numElements();

        approxs[r] = DWR->computeError(dualL,dualH,mp_def);
        gsInfo<<"Error = "<<approxs[r]<<"\n";
        estGoal[r] = numGoal[r]+approxs[r];

        gsVector<> Uz_pt(2);
        Uz_pt.setConstant(0.25);
        gsMatrix<> tmp;
        mp_def.patch(0).eval_into(Uz_pt,tmp);
        Uz[r] = tmp(2,0);

        if (adaptivity==0)
        {
            mp.uniformRefine();
        }
        else if (adaptivity > 0)
        {
            gsFileData<> fd_mesher(mesherOptionsFile);
            gsOptionList mesherOpts;
            fd_mesher.getFirst<gsOptionList>(mesherOpts);

            elErrors = DWR->computeErrorElements(dualL, dualH,mp_def);
            real_t error = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
            gsInfo<<"Accumulated error = "<<error<<"\n";

            gsDebugVar(*std::min_element(elErrors.begin(),elErrors.end()));
            for (std::vector<real_t>::iterator it = elErrors.begin(); it != elErrors.end(); it++)
            {
                *it = std::abs(*it);
            }

            real_t abserror = std::accumulate(elErrors.begin(),elErrors.end(),0.0);
            gsInfo<<"Absolute accumulated error = "<<abserror<<"\n";


            // Make container of the boxes
            gsHBoxContainer<2,real_t> markRef, markCrs;

            // gsHBoxContainer<2,real_t> elts;
            // index_t c = 0;
            // real_t newError = 0;
            // real_t oldError = 0;
            // for (index_t p=0; p < mp.nPieces(); ++p)
            // {
            // // index_t p=0;
            //     typename gsBasis<real_t>::domainIter domIt = mp.basis(p).makeDomainIterator();
            //     gsHDomainIterator<real_t,2> * domHIt = nullptr;
            //     domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
            //     GISMO_ENSURE(domHIt!=nullptr,"Domain not loaded");
            //     for (; domHIt->good(); domHIt->next())
            //     {
            //         gsHBox<2,real_t> box(domHIt,p);
            //         box.setAndProjectError(elErrors[c],mesherOpts.getInt("Convergence_alpha"),mesherOpts.getInt("Convergence_beta"));
            //         elts.add(box);
            //         markRef.add(box);
            //         oldError += box.error();
            //         newError += box.projectedErrorRef();
            //         c++;
            //     }
            // }

            // gsDebugVar(oldError);
            // gsDebugVar(newError);


            if (plot)
            {
                gsElementErrorPlotter<real_t> err_eh(mp.basis(0),elErrors);
                const gsField<> elemError_eh( mp.patch(0), err_eh, true );
                std::string fileName = dirname + "/" + "error_elem_ref" + util::to_string(r);
                gsWriteParaview<>( elemError_eh, fileName, 5000, true);
                fileName = "error_elem_ref" + util::to_string(r) + "0";
                errors.addTimestep(fileName,r,".vts");
                errors.addTimestep(fileName,r,"_mesh.vtp");
            }

            mesher.markRef_into(elErrors,markRef);
            mesher.refine(markRef);
            gsInfo<<"------------------------------------------------------------\n";
            gsInfo<<"-------------------Marked elements for refinement-----------\n";
            gsInfo<<"------------------------------------------------------------\n";
            gsInfo<<markRef<<"\n";
            gsInfo<<"------------------------------------------------------------\n";
            gsInfo<<"------------------------------------------------------------\n";
            gsInfo<<"------------------------------------------------------------\n";

            if (adaptivity>1)
            {
                mesher.markCrs_into(elErrors,markRef,markCrs);
                mesher.refine(markRef);
                mesher.unrefine(markCrs);
                gsInfo<<"------------------------------------------------------------\n";
                gsInfo<<"-------------------Marked elements for coarsening-----------\n";
                gsInfo<<"------------------------------------------------------------\n";
                gsInfo<<markRef<<"\n";
                gsInfo<<"------------------------------------------------------------\n";
                gsInfo<<"------------------------------------------------------------\n";
                gsInfo<<"------------------------------------------------------------\n";


                // gsDebugVar(markCrs);
            }
            mesher.rebuild();
        }
    }

    if (plot)
    {
        collection.save();
        errors.save();
    }

    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    gsInfo<<"Ref.\tApprox     \tEfficiency\tNumGoal   \texGoal    \tUz        \t#elements \n";
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    for(index_t r=0; r!=numRefine+1; r++)
    {
        gsInfo  <<std::setw(4 )<<std::left<<r<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<approxs[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<efficiencies[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<numGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<estGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<Uz[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<DoFs[r]<<"\n";
    }
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";


    if (write)
    {
        std::string filename;
        filename = dirname + "/" + "example_shell3D_DWR_NL_r" + std::to_string(numRefine) + "_e" + std::to_string(numElevate) + "_g" + std::to_string(goal) + "_C" + std::to_string(component);
        filename = filename + ".csv";
        std::ofstream file_out;
        file_out.open (filename);

        file_out<<"Ref,Approx,Efficiency,NumGoal,exGoal,Uz,DoFs\n";
        for(index_t r=0; r!=numRefine+1; r++)
        {
            file_out<<r<<","<<approxs[r]<<","<<efficiencies[r]<<","<<numGoal[r]<<","<<estGoal[r]<<","<<Uz[r]<<","<<DoFs[r]<<"\n";
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
    delete materialMatrix;
    return EXIT_SUCCESS;

}// end main
