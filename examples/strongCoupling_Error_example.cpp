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

#include <gsUtils/gsQuasiInterpolate.h>


#include <gsAssembler/gsExprAssembler.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot       = false;
    bool mesh       = false;
    bool stress     = false;
    bool write      = false;
    bool last       = false;
    bool info       = false;
    bool writeMatrix= false;
    bool nonlinear  = false;
    index_t numRefine  = 2;
    index_t numRefine0 = 1;
    index_t degree = 3;
    index_t smoothness = 2;
    index_t geometry = 1;
    index_t method = 0;
    std::string input;

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    std::string fn1,fn2,fn3;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );
    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addInt( "R", "preRefine", "Refinement before the loop.",  numRefine0);
    cmd.addInt( "m", "method", "Smoothing method to use", method );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("last", "last case only",last);
    cmd.addSwitch("writeMat", "Write projection matrix",writeMatrix);
    cmd.addSwitch( "info", "Print information", info );
    cmd.addSwitch( "nl", "Print information", nonlinear );

    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    gsBoundaryConditions<> bc;

    if (method==3 && degree-smoothness > 1)
        gsWarn<<"Smoothing method 3 with a degree "<<degree<<" and smoothness "<<smoothness<<" does not work, degree-smoothness=regularity=1 is used.\n";

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

    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!

    gsInfo<<pLoads;

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
    gsInfo<<"Reading reference values from "<<fn2<<" (ID=52) ...";
    if ( fd.hasId(52) )
        fd.getId(52,refValue);
    else
        refValue = gsMatrix<>::Zero(mp.geoDim(),refPoints.cols());
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

    gsMultiPatch<> geom = mp;

    // p-refine
    GISMO_ENSURE(degree>=mp.patch(0).degree(0),"Degree must be larger than or equal to the degree of the initial geometry, but degree = "<<degree<<" and the original degree = "<<mp.patch(0).degree(0));
    mp.degreeElevate(degree-mp.patch(0).degree(0));

    // h-refine each basis
    if(method!=3)
    {
        // h-refine each basis
        for (int r =0; r < numRefine0; ++r)
        {
            mp.uniformRefine(1,degree-smoothness);
        }
    }
    else
    {
        // Always regularity 1
        for (int r =0; r < numRefine0; ++r)
            mp.uniformRefine(1, degree-1);
    }
    numRefine -= numRefine0;

    if (last)
    {
        // h-refine
        if(method!=3)
        {
            // h-refine each basis
            for (int r =0; r < numRefine; ++r)
            {
                mp.uniformRefine(1,degree-smoothness);
            }
        }
        else
        {
            // Always regularity 1
            for (int r =0; r < numRefine; ++r)
                mp.uniformRefine(1, degree-1);
        }
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
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), linferr(numRefine+1),
        b2err(numRefine+1), b1err(numRefine+1), binferr(numRefine+1);

    gsVector<> numDofs(numRefine+1);
    gsMatrix<> refs(numRefine+1,3*refPoints.cols());

    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> solVector;

    gsMappedBasis<2,real_t> bb2;

    gsSparseMatrix<> global2local;
    gsStopwatch time;

    gsMultiBasis<> dbasis(mp);

    for( index_t r = 0; r<=numRefine; ++r)
    {

        gsInfo<<"--------------------------------------------------------------\n";
        time.restart();
        if (method==-1)
        {
            // identity map
            global2local.resize(dbasis.totalSize(),dbasis.totalSize());
            for (size_t k=0; k!=dbasis.totalSize(); ++k)
                global2local.coeffRef(k,k) = 1;
            geom = mp;
            gsInfo << "Basis Patch: " << dbasis.basis(0).component(0) << "\n";
        }
        else if (method==0)
        {
            gsMPBESSpline<2,real_t> cgeom(mp,3);
            gsMappedBasis<2,real_t> basis = cgeom.getMappedBasis();

            global2local = basis.getMapper().asMatrix();
            geom = cgeom.exportToPatches();
            auto container = basis.getBasesCopy();
            dbasis = gsMultiBasis<>(container,mp.topology());
        }
        else if (method==1)
        {
            geom = mp;
            gsDPatch<2,real_t> dpatch(geom);
            dpatch.matrix_into(global2local);

            global2local = global2local.transpose();
            geom = dpatch.exportToPatches();
            dbasis = dpatch.localBasis();
        }
        else if (method==2) // Pascal
        {
            mp.embed(2);
            gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
            approxC1.options().setSwitch("info",info);
            approxC1.options().setSwitch("plot",plot);
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
        else if (method==3) // Andrea
        {
            gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
            smoothC1.init();
            smoothC1.compute();

            global2local = smoothC1.getSystem();
            global2local = global2local.transpose();
            smoothC1.getMultiBasis(dbasis);
        }
        else if (method==4)
        {
            geom = mp;
            gsAlmostC1<2,real_t> almostC1(geom);
            almostC1.matrix_into(global2local);

            global2local = global2local.transpose();
            geom = almostC1.exportToPatches();
            dbasis = almostC1.localBasis();
        }
        else
            GISMO_ERROR("Option "<<method<<" for method does not exist");

        gsInfo<<"\tAssembly of mapping:\t"<<time.stop()<<"\t[s]\n";

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
        if (method==1)
            assembler.options().setInt("Continuity",-1);
        else if (method==2)
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


        // Nonlinear
        // Function for the Jacobian
        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
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
        Residual_t Residual = [&geom,&mp,&bb2,&assembler](gsVector<real_t> const &x)
        {
            gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
            GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
            solFull.resize(solFull.rows()/3,3);

            gsMappedSpline<2,real_t> mspline(bb2,solFull);
            gsFunctionSum<real_t> def(&mp,&mspline);

            assembler.assembleVector(def);
            return assembler.rhs();
        };

        // Linear solve
        solver.compute( matrix );
        solVector = solver.solve(vector);

        if (nonlinear)
        {
            real_t residual = vector.norm();
            real_t residual0 = residual;
            real_t residualOld = residual;
            gsVector<real_t> updateVector = solVector;
            gsVector<real_t> resVec = Residual(solVector);
            gsSparseMatrix<real_t> jacMat;
            for (index_t it = 0; it != 100; ++it)
            {
                jacMat = Jacobian(solVector);
                solver.compute(jacMat);
                updateVector = solver.solve(resVec); // this is the UPDATE
                solVector += updateVector;

                resVec = Residual(solVector);
                residual = resVec.norm();

                gsInfo<<"Iteration: "<< it
                   <<", residue: "<< residual
                   <<", update norm: "<<updateVector.norm()
                   <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
                   <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
                   <<"\n";

                residualOld = residual;

                if (updateVector.norm() < 1e-6)
                    break;
                else if (it+1 == it)
                    gsWarn<<"Maximum iterations reached!\n";
            }
        }


        /// Make a gsMappedSpline to represent the solution
        // 1. Get all the coefficients (including the ones from the eliminated BCs.)
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(solVector);

        // 2. Reshape all the coefficients to a Nx3 matrix
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        // 3. Make the mapped spline
        gsMappedSpline<2,real_t> mspline(bb2,solFull);

        gsFunctionSum<real_t> def(&mp,&mspline);

        gsField<> solField(geom, mspline,true);

        if (refPoints.cols()!=0)
        {

            gsMatrix<> ppoints(3,1), result;
            ppoints<<0.5,0,0.25;
            mp.patch(refPatches(0,0)).invertPoints(ppoints,result);
            for (index_t p=0; p!=refPoints.cols(); p++)
                // refs.block(r,p*mp.geoDim(),1,mp.geoDim()) = def.piece(refPatches(0,p)).eval(refPoints.col(p)).transpose();
                refs.block(r,p*mp.geoDim(),1,mp.geoDim()) = solField.value(refPoints.col(p),refPatches(0,p)).transpose();
        }

        numDofs[r] = assembler.numDofs();

        // h-refine
        if(method!=3)
        {
            mp.uniformRefine(1,degree-smoothness);
        }
        else
        {
            mp.uniformRefine(1, degree-1);
        }



        dbasis = gsMultiBasis<>(mp);
    }
    //! [Solver loop]

    gsInfo<<"numDoFs";
    for (index_t p=0; p!=refPoints.cols(); ++p)
        gsInfo<<",x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p);
    gsInfo<<"\n";

    for (index_t k=0; k<=numRefine; ++k)
    {
        gsInfo<<numDofs(k);
        for (index_t p=0; p!=refPoints.cols(); ++p)
        {
            gsInfo<<","<<refs(k,3*p)<<","<<refs(k,3*p+1)<<","<<refs(k,3*p+2)<<",";
        }
        gsInfo<<"\n";
    }

    //! [Export visualization in ParaView]
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
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);

        // 5. Plot the mapped spline on the deformed geometry
        gsField<> defField(geom, def,true);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( defField, "mp_def", 1000, true);

        // QUASI INTERPOLATION
        // /*

        // gsMultiPatch<> mpatches = mbasis.exportToPatches(tmp);
        gsMultiPatch<> mp2;
        for (size_t p = 0; p!=mp.nPatches(); p++)
        {
            gsMatrix<> coefs;
            gsQuasiInterpolate<real_t>::localIntpl(mp.basis(p), mspline.piece(p), coefs);
            mp2.addPatch(mp.basis(p).makeGeometry( give(coefs) ));
        }

        gsField<> solfield(mp,mp2,true);
        gsWriteParaview(solfield,"solfield");
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}