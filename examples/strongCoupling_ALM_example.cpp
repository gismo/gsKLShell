/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsUnstructuredSplines/gsMPBESBasis.h>
#include <gsUnstructuredSplines/gsMPBESSpline.h>
#include <gsUnstructuredSplines/gsDPatch.h>
#include <gsUnstructuredSplines/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/gsC1SurfSpline.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsUtils/gsQuasiInterpolate.h>


#include <gsAssembler/gsExprAssembler.h>

// #ifdef GISMO_STRUCTURALANALYSIS
#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>
#include <gsStructuralAnalysis/gsALMConsistentCrisfield.h>
// #endif

using namespace gismo;

int main(int argc, char *argv[])
{
// #ifndef GISMO_STRUCTURALANALYSIS
//     GISMO_ERROR("This code should be compiled with cmake flag -DGISMO_STRUCTURALANALYSIS=ON");
// #else
    bool plot       = false;
    bool mesh       = false;
    bool last       = false;
    bool info       = false;
    bool writeMatrix= false;
    bool nonlinear  = false;
    bool SingularPoint = false;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    index_t numRefine  = 2;
    index_t numElevate = 1;
    index_t geometry = 1;
    index_t smoothing = 0;
    std::string input;

    int step          = 10;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method)
    real_t dL         = 0.5; // Arc length
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;
    real_t tau = 1e4;

    std::string fn1,fn2,fn3;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "s", "smooth", "Smoothing method to use",  smoothing );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("last", "last case only",last);
    cmd.addSwitch("writeMat", "Write projection matrix",writeMatrix);
    cmd.addSwitch( "info", "Print information", info );
    cmd.addSwitch( "nl", "Print information", nonlinear );

    cmd.addReal("F","factor", "factor for bifurcation perturbation", tau);
    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addSwitch("bifurcation", "Compute singular points and bifurcation paths", SingularPoint);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);

    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    gsBoundaryConditions<> bc;

    /*
        to do:
        - remove hard-coded IDs from XML reader
     */

    gsFileData<> fd;
    gsInfo<<"Reading geometry from "<<fn1<<"...";
    gsReadFile<>(fn1, mp);
    if (mp.nInterfaces()==0 && mp.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it...";
        mp.computeTopology();
    }
    gsInfo<<"Finished\n";

    fd.read(fn2);
    index_t num = 0;
    gsInfo<<"Reading BCs from "<<fn2<<"...";
    num = fd.template count<gsBoundaryConditions<>>();
    GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
    fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
    gsInfo<<"Finished\n";

    bc.setGeoMap(mp);

    // Loads
    gsFunctionExpr<> force, pressure;
    gsInfo<<"Reading force function from "<<fn2<<" (ID=21) ...";
    fd.getId(21, force); // id=1: source function
    gsInfo<<"Finished\n";
    // fd.getId(22, pressure); // id=1: source function ------- TO DO!
    // gsInfo<<"Pressure function "<< force << "\n";

    // Loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsInfo<<"Reading point load locations from "<<fn2<<" (ID=30) ...";
    fd.getId(30,points);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point loads from "<<fn2<<" (ID=31) ...";
    fd.getId(31,loads);
    gsInfo<<"Finished\n";
    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), 0 ); // in parametric domain!

    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    gsInfo<<"Reading thickness from "<<fn2<<" (ID=10) ...";
    fd.getId(10,t);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Young's Modulus from "<<fn2<<" (ID=11) ...";
    fd.getId(11,E);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Poisson ratio from "<<fn2<<" (ID=12) ...";
    fd.getId(12,nu);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading density from "<<fn2<<" (ID=13) ...";
    fd.getId(13,rho);
    gsInfo<<"Finished\n";

    if (mp.domainDim()==2)
        mp.embed(3);

    gsMultiPatch<> geom = mp;

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    if(smoothing!=3)
    {
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine();
    }
    else
    {
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine(1, mp.patch(0).degree(0) -1);
    }

    if (plot) gsWriteParaview(mp,"mp",1000,true,false);
    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    //     gsDebugVar(mp.patch(p));

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters,rho);

    gsThinShellAssembler<3, real_t, true> assembler;

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> solVector;

    gsMappedBasis<2,real_t> bb2;

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;

    gsMultiBasis<> dbasis(mp);
    gsWriteParaview(mp.basis(0),"basis",1000);

    if (smoothing==0)
    {
        gsMPBESSpline<2,real_t> cgeom(mp,3);
        gsMappedBasis<2,real_t> basis = cgeom.getMappedBasis();

        global2local = basis.getMapper().asMatrix();
        geom = cgeom.exportToPatches();
        auto container = basis.getBasesCopy();
        dbasis = gsMultiBasis<>(container,mp.topology());
    }
    else if (smoothing==1)
    {
        gsDPatch<2,real_t> dpatch(mp);
        dpatch.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = dpatch.exportToPatches();
        dbasis = dpatch.localBasis();
    }
    else if (smoothing==2) // Pascal
    {
        mp.embed(2);
        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        approxC1.options().setSwitch("info",info);
        // approxC1.options().setSwitch("plot",plot);
        // approxC1.options().setInt("gluingDataDegree",)
        // approxC1.options().setInt("gluingDataRegularity",)

        gsDebugVar(approxC1.options());

        approxC1.init();
        approxC1.compute();
        mp.embed(3);

        global2local = approxC1.getSystem();
        global2local = global2local.transpose();
        global2local.pruned(1,1e-10);
        geom = mp;
        approxC1.getMultiBasis(dbasis);
    }
    else if (smoothing==3) // Andrea
    {
        gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
        smoothC1.init();
        smoothC1.compute();

        global2local = smoothC1.getSystem();
        global2local = global2local.transpose();
        smoothC1.getMultiBasis(dbasis);
    }
    else
        GISMO_ERROR("Option "<<smoothing<<" for smoothing does not exist");

    if (writeMatrix)
    {
        gsWrite(global2local,"mat");
        //gsWrite(geom,"geom");
        //gsWrite(dbasis,"dbasis");
    }

    bb2.init(dbasis,global2local);
    // gsMappedSpline<2,real_t> mspline(bb2,coefs);
    // geom = mspline.exportToPatches();

    assembler = gsThinShellAssembler<3, real_t, true>(geom,dbasis,bc,force,&materialMatrix);
    if (smoothing==1)
        assembler.options().setInt("Continuity",-1);
    else if (smoothing==2)
        assembler.options().setInt("Continuity",-1);
    assembler.setSpaceBasis(bb2);
    assembler.setPointLoads(pLoads);
    // gsOptionList options = assembler.options();
    // options.setInt("Continuity",1);
    // assembler.setOptions(options);

    // Initialize the system
    // Linear
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();

    // gsDebugVar(matrix.toDense());
    // gsDebugVar(vector.transpose());

    assembler.assemble();
    gsVector<> Force = assembler.rhs();
    // Nonlinear
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >                                     Residual_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t) >                             ALResidual_t;
    Jacobian_t Jacobian = [&mp,&bb2,&assembler](gsVector<real_t> const &x)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&mp,&mspline);

        assembler.assembleMatrix(def);
        // gsSparseMatrix<real_t> m =
        return assembler.matrix();
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&Force,&geom,&mp,&bb2,&assembler](gsVector<real_t> const &x, real_t lam)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&mp,&mspline);

        assembler.assembleVector(def);
        gsVector<real_t> Fint = -(assembler.rhs() - Force);
        gsVector<real_t> result = Fint - lam * Force;
        return result; // - lam * force;
    };

    gsALMBase<real_t> * arcLength;
    if (method==0)
      arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (method==1)
      arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (method==2)
      arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else if (method==3)
      arcLength = new gsALMConsistentCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
      GISMO_ERROR("Method unknown");

    arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
    arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    arcLength->options().setReal("Length",dL);
    // arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
    arcLength->options().setInt("AdaptiveIterations",5);
    arcLength->options().setReal("Perturbation",tau);
    // arcLength->options().setReal("Scaling",0.0);
    arcLength->options().setReal("Tol",tol);
    arcLength->options().setReal("TolU",tolU);
    arcLength->options().setReal("TolF",tolF);
    arcLength->options().setSwitch("Verbose",true);
    if (quasiNewtonInt>0)
    {
      quasiNewton = true;
      arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength->options().setSwitch("Quasi",quasiNewton);

    gsInfo<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR
    bool bisected = false;
    real_t dL0 = dL;
    real_t dLb = dL;

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";
    gsParaviewCollection collection(dirname + "/" + output);

    for (index_t k=0; k<step; k++)
    {
        gsInfo<<"Load step "<< k<<"\n";
        // assembler->constructSolution(solVector,solution);
        arcLength->step();

        // gsInfo<<"m_U = "<<arcLength->solutionU()<<"\n";
        if (!(arcLength->converged()))
        {
          gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
          dL = dL / 2.;
          arcLength->setLength(dL);
          arcLength->setSolution(Uold,Lold);
          bisected = true;
          k -= 1;
          continue;
        }
        arcLength->computeStability(arcLength->solutionU(),quasiNewton);

        if (arcLength->stabilityChange())
        {
            gsInfo<<"Bifurcation spotted!"<<"\n";
            arcLength->computeSingularPoint(1e-4, 5, Uold, Lold, 1e-7, 0, false);
            arcLength->switchBranch();
            dL0 = dLb = dL;
            arcLength->setLength(dLb);
        }
        indicator = arcLength->indicator();

        solVector = arcLength->solutionU();
        Uold = solVector;
        Lold = arcLength->solutionL();

        if (plot)
        {
            /// Make a gsMappedSpline to represent the solution
            // 1. Get all the coefficients (including the ones from the eliminated BCs.)
            gsMatrix<real_t> solFull = assembler.fullSolutionVector(solVector);

            // 2. Reshape all the coefficients to a Nx3 matrix
            GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
            solFull.resize(solFull.rows()/3,3);

            // 3. Make the mapped spline
            gsMappedSpline<2,real_t> mspline(bb2,solFull);

            gsFunctionSum<real_t> def(&mp,&mspline);

            // 4. Plot the mapped spline on the original geometry
            gsField<> solField(geom, mspline,true);

            std::string fileName = dirname + "/" + output + util::to_string(k);
            gsWriteParaview<>(solField, fileName, 1000,mesh);
            for (index_t p = 0; p!=mp.nPatches(); p++)
            {
                fileName = output + util::to_string(k);
                collection.addTimestep(fileName,p,k,".vts");
                if (mesh)
                    collection.addTimestep(fileName,p,k,"_mesh.vtp");
            }
        }

        if (!bisected)
        {
          dL = dL0;
          arcLength->setLength(dL);
        }
        bisected = false;
    }

    if (plot)
    {
        collection.save();
    }

    delete arcLength;

    return EXIT_SUCCESS;
// #endif


}