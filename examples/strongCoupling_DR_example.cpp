/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESSpline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsControlDisplacement.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{
// #ifndef GISMO_STRUCTURALANALYSIS
//     GISMO_ERROR("This code should be compiled with cmake flag -DGISMO_STRUCTURALANALYSIS=ON");
// #else

    bool plot       = false;
    bool stress     = false;
    bool write      = false;
    bool mesh       = false;
    bool last       = false;
    bool info       = false;
    bool writeMatrix= false;
    bool nonlinear  = false;
    index_t numRefine  = 2;
    index_t degree = 3;
    index_t smoothness = 2;
    index_t geometry = 1;
    index_t method = 0;
    std::string input;
    real_t perturbation = 0;

    index_t step      = 1e1;
    index_t maxIt     = 1e3;
    real_t tolE = 1e-3;
    real_t tolF = 1e-5;
    real_t alpha = 1.0;
    real_t damping = 0.1;
    real_t dt = 0.1;

    real_t Lmax       = std::numeric_limits<real_t>::max();

    index_t verbose = 0;

    std::string fn1,fn2,fn3;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";
    std::string dirname = "DynamicRelaxationResults";

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    std::vector<index_t> modevec;

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );

    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addString( "o", "out", "Dir name of the output",  dirname );
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addInt( "m", "method", "Smoothing method to use", method );

    cmd.addInt( "n", "maxIt", "maxIt",  maxIt );
    cmd.addInt( "N", "step", "number of load steps",  step );

    cmd.addReal( "d", "dt", "dt",  dt );
    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );
    cmd.addReal("P","perturbation", "perturbation factor", perturbation);
    cmd.addReal("b","Lmax", "Maximum L", Lmax);
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);

    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("last", "last case only",last);
    cmd.addSwitch("writeMat", "Write projection matrix",writeMatrix);
    cmd.addSwitch( "info", "Print information", info );


    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
    GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");
    if (method==3)
        GISMO_ENSURE(smoothness>=1 || smoothness <= degree-2,"Exact C1 method only works for smoothness <= p-2, but smoothness="<<smoothness<<" and p-2="<<degree-2);
    if (method==2 || method==3)
        GISMO_ENSURE(degree > 2,"Degree must be larger than 2 for the approx and exact C1 methods, but it is "<<degree);


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
    gsMatrix<index_t> pid_ploads;
    if ( fd.hasId(30) )
        fd.getId(30,points);
    if ( fd.hasId(31) )
        fd.getId(31,loads);

    if ( fd.hasId(32) )
        fd.getId(32,pid_ploads);
    else
        pid_ploads = gsMatrix<index_t>::Zero(1,points.cols());

    if (points.cols()>0)
    {
        // gsMatrix<> pointLoadPoints(3,points.cols());
        // for (index_t k =0; k!=points.cols(); k++)
        //     pointLoadPoints.col(k) = mp.patch(pid_ploads.at(k)).eval(points.col(k));
        // gsWriteParaviewPoints(pointLoadPoints,"pointLoadPoints");

        for (index_t k =0; k!=points.cols(); k++)
            pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!

    }

    gsInfo<<pLoads;

    // // Loads
    // gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    // gsMatrix<> points,loads;
    // gsInfo<<"Reading point load locations from "<<fn2<<" (ID=30) ...";
    // fd.getId(30,points);
    // gsInfo<<"Finished\n";
    // gsInfo<<"Reading point loads from "<<fn2<<" (ID=31) ...";
    // fd.getId(31,loads);
    // gsInfo<<"Finished\n";
    // for (index_t k =0; k!=points.cols(); k++)
    //     pLoads.addLoad(points.col(k), loads.col(k), 0 ); // in parametric domain!

    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations from "<<fn2<<" (ID=50) ...";
    if ( fd.hasId(50) )
        fd.getId(50,refPoints);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches from "<<fn2<<" (ID=51) ...";
    if ( fd.hasId(51) )
        fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    GISMO_ENSURE(refPatches.cols()==refPoints.cols(),"Number of reference points and patches do not match");

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

    if (mp.geoDim()==2)
        mp.embed(3);

    gsOptionList solverOptions;
    fd.read(fn3);
    fd.template getFirst<gsOptionList>(solverOptions);

    gsMultiPatch<> geom = mp;

    GISMO_ENSURE(degree>=mp.patch(0).degree(0),"Degree must be larger than or equal to the degree of the initial geometry, but degree = "<<degree<<" and the original degree = "<<mp.patch(0).degree(0));
    mp.degreeElevate(degree-mp.patch(0).degree(0));

    // h-refine
    for (index_t r =0; r < numRefine; ++r)
        mp.uniformRefine(1,degree-smoothness);


    if (plot) gsWriteParaview(mp,"mp",1000,true,false);

    if (perturbation!=0)
        for (size_t p=0; p!=mp.nPatches(); p++)
        {
            mp.patch(p).coefs().col(2)  = gsMatrix<>::Random(mp.patch(p).coefs().rows(),1);
            mp.patch(p).coefs().col(2) *= perturbation;
        }

    if (plot) gsWriteParaview(mp,"mp_perturbed",1000,true,false);

    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    //     gsDebugVar(mp.patch(p));

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsDebugVar(rho);
    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters,rho);

    gsThinShellAssembler<3, real_t, true> assembler;

    //! [Solver loop]
    gsVector<> solVector;

    gsMappedBasis<2,real_t> bb2;

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;

    gsMultiBasis<> dbasis(mp);

    if (method==-1)
    {
        // identity map
        global2local.resize(dbasis.totalSize(),dbasis.totalSize());
        for (size_t k=0; k!=dbasis.totalSize(); ++k)
            global2local.coeffRef(k,k) = 1;
        geom = mp;
        bb2.init(dbasis,global2local);
    }
    else if (method==0)
    {
        gsMPBESSpline<2,real_t> cgeom(mp,3);
        gsMappedBasis<2,real_t> basis = cgeom.getMappedBasis();

        global2local = basis.getMapper().asMatrix();
        geom = cgeom.exportToPatches();
        auto container = basis.getBasesCopy();
        dbasis = gsMultiBasis<>(container,mp.topology());
        bb2.init(dbasis,global2local);
    }
    else if (method==1)
    {
        geom = mp;
        gsDPatch<2,real_t> dpatch(geom);
        dpatch.options().setSwitch("SharpCorners",false);
        dpatch.compute();
        dpatch.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = dpatch.exportToPatches();
        dbasis = dpatch.localBasis();
        bb2.init(dbasis,global2local);
    }
    else if (method==2) // Pascal
    {
        // The approx. C1 space
        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        // approxC1.options().setSwitch("info",info);
        // approxC1.options().setSwitch("plot",plot);
        approxC1.options().setSwitch("interpolation",true);
        approxC1.options().setInt("gluingDataDegree",-1);
        approxC1.options().setInt("gluingDataSmoothness",-1);
        approxC1.update(bb2);
        geom = mp;
    }
    else if (method==3) // Andrea
    {
        dbasis = gsMultiBasis<>(mp);
        gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
        smoothC1.init();
        smoothC1.compute();

        global2local = smoothC1.getSystem();
        global2local = global2local.transpose();
        smoothC1.getMultiBasis(dbasis);
        bb2.init(dbasis,global2local);
        geom = mp;
    }
    else if (method==4)
    {
        gsAlmostC1<2,real_t> almostC1(mp);
        almostC1.compute();
        almostC1.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = almostC1.exportToPatches();
        dbasis = almostC1.localBasis();
        bb2.init(dbasis,global2local);
    }
    else
        GISMO_ERROR("Option "<<method<<" for method does not exist");

    if (writeMatrix)
    {
        gsWrite(global2local,"mat");
        //gsWrite(geom,"geom");
        //gsWrite(dbasis,"dbasis");
    }

    if (plot) gsWriteParaview(geom,"geom",1000,true,false);

    assembler = gsThinShellAssembler<3, real_t, true>(geom,dbasis,bc,force,&materialMatrix);
    assembler.setOptions(solverOptions); //sets solverOptions
    assembler.options().setInt("Continuity",-1);
    assembler.options().setReal("WeakDirichlet",bcDirichlet);
    assembler.options().setReal("WeakClamped",bcClamped);
    assembler.setSpaceBasis(bb2);
    assembler.setPointLoads(pLoads);

    // Initialize the system
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();
    gsVector<> F = assembler.rhs();
    assembler.assembleMass(true);
    gsVector<> M = assembler.rhs();

    // Nonlinear
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t) > ALResidual_t;
    // Function for the Residual
    ALResidual_t ALResidual = [&geom,&bb2,&assembler,&F /*,&Force_const*/](gsVector<real_t> const &x, real_t lam)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleVector(def);
        gsVector<real_t> Fint = -(assembler.rhs() - F);
        gsVector<real_t> result = Fint - lam * F/*- Force_const*/;
        return result; // - lam * force;
    };

    gsStaticDR<real_t> DRM(M,F,ALResidual);
    gsOptionList DROptions = DRM.options();
    DROptions.setReal("damping",damping);
    DROptions.setReal("alpha",alpha);
    DROptions.setInt("maxIt",maxIt);
    DROptions.setReal("tol",tolF);
    DROptions.setReal("tolE",tolE);
    DROptions.setInt("verbose",verbose);
    DRM.setOptions(DROptions);

    std::string output = "solution";
    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection membraneStressCollection(dirname + "/MembraneStress");
    gsParaviewCollection flexuralStressCollection(dirname + "/FlexuralStress");

    gsALMOutput<real_t> numWriter(dirname + "/pointdata.csv",refPoints);
    std::vector<std::string> headers = {"u_x","u_y","u_z"};
    if (write)
        numWriter.init(headers);

    DRM.initialize();
    gsControlDisplacement<real_t> control(&DRM);

    index_t k=0;
    index_t L=0;
    while (k < step && L < Lmax)
    {
        gsInfo<<"Load step "<< k<<", Load = "<<L<<"\n";
        control.step(1.0);
        gsInfo<<"Step finished in "<<DRM.iterations()<<" iterations";

        solVector = control.solutionU();

        if (plot || write || stress)
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

            if (plot)
            {
                std::string fileName = dirname + "/" + output + util::to_string(k) + "_";
                gsWriteParaview<>(solField, fileName, 1000,mesh);
                for (index_t p = 0; p!=mp.nPatches(); p++)
                {
                    fileName = output + util::to_string(k) + "_";
                    collection.addTimestep(fileName,p,k,".vts");
                    if (mesh)
                        collection.addTimestep(fileName,p,k,"_mesh.vtp");
                }
            }

            if (write)
            {
                if (refPoints.cols()!=0)
                {
                    gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
                    for (index_t p=0; p!=refPoints.cols(); p++)
                        pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
                    numWriter.add(pointResults,control.solutionL());
                }
            }

            if (stress)
            {
                std::string fileName;
                gsPiecewiseFunction<> membraneStresses;
                assembler.constructStress(def,membraneStresses,stress_type::membrane);
                fileName = dirname + "/MembraneStress" + util::to_string(k) + "_";
                gsWriteParaview(def,membraneStresses,fileName,1000);
                for (index_t p = 0; p!=mp.nPatches(); p++)
                {
                    fileName = "MembraneStress" + util::to_string(k) + "_";
                    membraneStressCollection.addTimestep(fileName,p,k,".vts");
                    if (mesh)
                        membraneStressCollection.addTimestep(fileName,p,k,"_mesh.vtp");
                }

                gsPiecewiseFunction<> flexuralStresses;
                assembler.constructStress(def,flexuralStresses,stress_type::flexural);
                fileName = dirname + "/FlexuralStresses" + util::to_string(k) + "_";
                gsWriteParaview(def,flexuralStresses,fileName,1000);
                for (index_t p = 0; p!=mp.nPatches(); p++)
                {
                    fileName = "FlexuralStress" + util::to_string(k) + "_";
                    flexuralStressCollection.addTimestep(fileName,p,k,".vts");
                    if (mesh)
                        flexuralStressCollection.addTimestep(fileName,p,k,"_mesh.vtp");
                }

            }
        }

        L = control.solutionL();
        k++;
    }

    if (plot)
    {
        collection.save();
    }
    if (stress)
    {
        membraneStressCollection.save();
        flexuralStressCollection.save();
    }

    return EXIT_SUCCESS;
// #endif
}
