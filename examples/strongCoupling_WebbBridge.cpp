/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsUnstructuredSplines/src/gsDPatch.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsAssembler/gsExprAssembler.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

#include <gsKLShell/gsThinShellUtils.h>
#include <gsKLShell/getMaterialMatrix.h>

#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsStaticNewton.h>
#include <gsStructuralAnalysis/gsControlDisplacement.h>
#include <gsUtils/gsL2Projection.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot       = false;
    bool mesh       = false;
    bool stress     = false;
    bool write      = false;
    index_t verbose = 0;

    index_t numRefine  = 1;
    index_t numElevate = 1;

    bool Compressibility = false;
    index_t material = 0;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t E_modulus = 1.0e9;
    real_t PoissonRatio = 1.0/3.0;
    real_t Density = 1.0;
    real_t thickness = 1.0e-3;

    index_t maxIt     = 1e3;
    // Arc length method options
    real_t tol        = 1e-6;

    real_t alpha = 1.0;
    real_t damping = 0.1;

    real_t dt = 0.1;


    gsCmdLine cmd("Composite basis tests.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("comp", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);

    cmd.addInt( "N", "maxit", "maxit",  maxIt );

    cmd.addReal( "d", "dt", "dt",  dt );
    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);



    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t H = 25e-3;
    real_t W  = 75e-3;
    real_t W1 = 2./5.*W;
    real_t W2 = 1./5.*W;

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,W1,H));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(W1,0,W1+W2,H));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(W1+W2,0,W,H));

    mp.degreeReduce(1);
    mp.patch(0).uniformRefine(1,1,0);
    mp.patch(2).uniformRefine(1,1,0);

    mp.embed(3);
    mp.computeTopology();
    gsDebugVar(mp.patch(0));

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    index_t ampl = -6;
    gsFunctionExpr<> perturbation("10^(" + std::to_string(ampl) + ")" + "*sin(10*x)*sin(10*y)",2);
    gsInfo<<perturbation<<"\n";
    gsMatrix<> coefs;
    gsL2Projection<real_t>::projectFunction(dbasis,perturbation,mp,coefs);
    // coefs = mp.coefs().col(2) + coefs;
    index_t offset = 0;
    index_t patchsize = 0;
    for (index_t p=0; p!=mp.nPatches(); p++)
    {
        patchsize = mp.patch(p).coefs().rows();
        mp.patch(p).coefs().col(2) = coefs.block(offset,0,patchsize,1);
        offset += patchsize;
    }

    if (plot) gsWriteParaview<>(mp,"mp");

    // Set the deformed configuration
    mp_def = mp;


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";
    //! [Make geometry and refine/elevate]

    // Reference points
    gsMatrix<index_t> refPatches(1,1);
    refPatches<<1;
    gsMatrix<> refPars(2,1);
    refPars<<0.5,0;

    //! [Make material functions]
    // Linear isotropic material model and Neo-Hookean material
    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    gsConstantFunction<> mu1(749.18,3);
    gsConstantFunction<> alpha1(17.14,3);

    // Define parameters vector depending on material law
    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK & Composites
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==1 || material==2) // NH & NH_ext
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==4) // OG
    {
      parameters.resize(4);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
    }
    else
        GISMO_ERROR("No material law known");
    //! [Make material functions]

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    // Make gsMaterialMatrix depending on the user-defined choices

    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    //! [Define jacobian and residual]
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;

    //! [Set boundary conditions]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsConstantFunction<> DY(-6e-3,3);
    bc.addCondition(0,boundary::west , condition_type::dirichlet, 0, 0 ,false, -1);
    // bc.addCondition(1,boundary::south, condition_type::collapsed, 0, 0 ,false, 1);
    bc.addCondition(1,boundary::south, condition_type::dirichlet, &DY, 0 ,false, 1);
    bc.addCondition(2,boundary::east , condition_type::dirichlet, 0, 0 ,false, -1);

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsVector<> points(2);
    points<<0.5,0;
    gsVector<> loads(3);
    loads<<0,-1e3,0;
    pLoads.addLoad(points, loads, 1 ); // in parametric domain!


    gsMultiPatch<> geom;
    gsMappedBasis<2,real_t> bb2;
    gsSparseMatrix<> global2local;

    geom = mp;
    gsDPatch<2,real_t> dpatch(geom);
    dpatch.compute();
    dpatch.matrix_into(global2local);

    global2local = global2local.transpose();
    geom = dpatch.exportToPatches();
    dbasis = dpatch.localBasis();
    bb2.init(dbasis,global2local);

    gsThinShellAssembler<3, real_t, true> assembler(geom,dbasis,bc,force,materialMatrix);
    assembler.setSpaceBasis(bb2);
    assembler.setPointLoads(pLoads);

    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();
    gsVector<> F = assembler.rhs();
    assembler.assembleMass(true);
    gsVector<> M = assembler.rhs();

    // Nonlinear
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t) > ALResidual_t;
    Jacobian_t Jacobian = [&geom,&bb2,&assembler](gsVector<real_t> const &x)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleMatrix(def);
        // gsSparseMatrix<real_t> m =
        return assembler.matrix();
    };
    // Function for the Residual
    Residual_t Residual = [&geom,&bb2,&assembler](gsVector<real_t> const &x)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleVector(def,false);
        return assembler.rhs();
    };

    // Function for the ALResidual
    ALResidual_t ALResidual = [&DY,&geom,&bb2,&assembler](gsVector<real_t> const &x, real_t lam)
    {
        DY.setValue(lam,3);
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleVector(def,false);
        return assembler.rhs();
    };


    gsStaticDR<real_t> DRM(M,F,ALResidual);
    gsOptionList DROptions = DRM.options();
    DROptions.setReal("damping",damping);
    DROptions.setReal("alpha",alpha);
    DROptions.setInt("maxIt",maxIt);
    DROptions.setReal("tol",1e-2);
    DROptions.setReal("tolE",1e-4);
    DROptions.setInt("verbose",verbose);
    DRM.setOptions(DROptions);
    DRM.initialize();
    gsControlDisplacement<real_t> controlDR(&DRM);

    // gsStaticNewton<real_t> NRM(K,F,Jacobian,ALResidual);
    // gsOptionList NROptions = NRM.options();
    // NROptions.setInt("verbose",verbose);
    // NROptions.setInt("maxIt",maxIt);
    // NROptions.setReal("tol",1e-2);
    // NRM.setOptions(NROptions);
    // NRM.initialize();
    // gsControlDisplacement<real_t> controlDR(&NRM);


    index_t maxSteps = 20;
    real_t dL = -6e-3/maxSteps;
    gsParaviewCollection collection("solution");
    for (index_t step = 0; step != maxSteps; step++)
    {
        gsInfo<<"Step "<<step<<"\n";
        controlDR.step(dL);

        solVector = DRM.solution();
        // ! [Solver loop]

        /// Make a gsMappedSpline to represent the solution
        // 1. Get all the coefficients (including the ones from the eliminated BCs.)
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(solVector);
        gsMatrix<real_t> solZero = solFull;
        solZero.setZero();

        // 2. Reshape all the coefficients to a Nx3 matrix
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solZero.resize(solZero.rows()/3,3);
        solFull.resize(solFull.rows()/3,3);

        // 3. Make the mapped spline
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsInfo <<"Maximum deformation coef: "
               << solFull.colwise().maxCoeff() <<".\n";
        gsInfo <<"Minimum deformation coef: "
               << solFull.colwise().minCoeff() <<".\n";
        gsDebugVar(solFull.colwise().maxCoeff());

        gsField<> solField(geom, mspline,true);

        gsInfo << "Number of Dofs: " << assembler.numDofs() << "\n";

        //! [Export visualization in ParaView]
        if (plot)
        {
            // 4. Plot the mapped spline on the original geometry
            gsField<> solField(geom, mspline,true);
            gsInfo<<"Plotting in Paraview...\n";
            std::string fileName = "solution_" + util::to_string(step) + "_";
            gsWriteParaview<>( solField,fileName, 1000, mesh);
            for (index_t p=0; p!=mp.nPatches(); p++)
            {
                collection.addPart(fileName + "vts",step,"solution",p);
            }

            // // 4. Plot the mapped spline on the original geometry
            // gsField<> solField2(mp_def, def,true);
            // gsInfo<<"Plotting in Paraview...\n";
            // gsWriteParaview<>( solField2, "Deformation_", 1000, true);


            // // 5. Plot the mapped spline on the deformed geometry
            // gsField<> defField(def, def,true);
            // gsInfo<<"Plotting in Paraview...\n";
            // gsWriteParaview<>( defField, "mp_def", 1000, true);

            // gsMultiPatch<> mpatches = mbasis.exportToPatches(tmp);

            // gsField<> solfield(geom,def,true);
            // gsWriteParaview(solfield,"solfield",1000,true);
        }

    }

    collection.save();

    // gsMappedSpline<2,real_t> mspline_ori(bb2,solZero);
    // gsMappedSpline<2,real_t> mspline_def(bb2,solFull);

    // gsFunctionSum<real_t> ori(&geom,&mspline_ori);
    // gsFunctionSum<real_t> def(&geom,&mspline_def);

    // if (stress)
    // {
    //     gsPiecewiseFunction<> membraneStresses;
    //     gsDebugVar("MembraneStress construction");
    //     assembler.constructStress(ori,def,membraneStresses,stress_type::membrane);
    //     gsWriteParaview(ori,membraneStresses,"MembraneStress",5000);

    //     gsPiecewiseFunction<> membraneStressesVM;
    //     gsDebugVar("MembraneStress (VM) construction");
    //     assembler.constructStress(ori,def,membraneStressesVM,stress_type::von_mises_membrane);
    //     gsWriteParaview(ori,membraneStressesVM,"MembraneStressVM",5000);

    //     gsPiecewiseFunction<> flexuralStresses;
    //     gsDebugVar("FlexuralStress construction");
    //     assembler.constructStress(ori,def,flexuralStresses,stress_type::flexural);
    //     gsWriteParaview(geom,flexuralStresses,"FlexuralStress",5000);
    // }
    // if (write)
    // {
    //     std::vector<std::string> headers = {"u_x","u_y","u_z"};
    //     gsStaticOutput<real_t> ptsWriter(dirname + "pointcoordinates.csv",refPoints);
    //     ptsWriter.init(headers);
    //     gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
    //     for (index_t p=0; p!=refPoints.cols(); p++)
    //         pointResults.col(p) = mp.piece(refPatches(0,p)).eval(refPoints.col(p));
    //     ptsWriter.add(pointResults);

    //     gsStaticOutput<real_t> numWriter(dirname + "numerical.csv",refPoints);
    //     numWriter.init(headers);
    //     for (index_t p=0; p!=refPoints.cols(); p++)
    //         pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
    //     numWriter.add(pointResults);

    //     gsStaticOutput<real_t> refWriter(dirname + "reference.csv",refPoints);
    //     refWriter.init(headers);
    //     refWriter.add(refValue);
    // }
    // //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}
