/** @file strongCoupling_example.cpp

    @brief Multi-patch shell analysis via unstructured splines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#ifdef gsUnstructuredSplines_ENABLED
#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESSpline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>
#endif

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsUtils/gsQuasiInterpolate.h>

#include <gsAssembler/gsExprAssembler.h>

using namespace gismo;

void writeToFile(const std::string & bufferstring, std::ofstream & file, const std::string & name)
{
    file.open(name);
    file<<bufferstring;
    file.close();
    gsInfo<<"Data written to "<<name<<"\n";
}

#ifdef gsUnstructuredSplines_ENABLED
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool mesh       = false;
    bool stress     = false;
    bool write      = false;
    bool nonlinear  = false;
    bool homogeneous= false;

    index_t numRefine  = 2;
    index_t degree = 3;
    index_t smoothness = 2;

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    std::string fn1,fn2,fn3,fn4;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";
    std::string dirname = ".";

    index_t method = 0;

    gsCmdLine cmd("Multi-patch shell analysis via unstructured splines.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );

    cmd.addInt( "m", "method", "Smoothing method to use", method );

    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addString( "o", "out", "Dir name of the output",  dirname );

    cmd.addString( "","fullbasis", "File containing all basis details: 1) a multi-patch containing the geometry, 2) a multi-basis containing the local basis, and 3) a sparse matrix being a basis mapper",  fn4 );

    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );

    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch( "nl", "Nonlinear analysis", nonlinear );
    cmd.addSwitch("homogeneous", "homogeneous dirichlet BCs",homogeneous);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
    GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");
    if (method==3)
        GISMO_ENSURE(smoothness>=1 || smoothness <= degree-2,"Exact C1 method only works for smoothness <= p-2, but smoothness="<<smoothness<<" and p-2="<<degree-2);
    if (method==2 || method==3)
        GISMO_ENSURE(degree > 2,"Degree must be larger than 2 for the approx and exact C1 methods, but it is "<<degree);
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> fd;

    // Geometry
    gsMultiPatch<> mp;

    gsInfo<<"Reading geometry from "<<fn1<<"...";
    gsReadFile<>(fn1, mp);
    if (mp.nInterfaces()==0 && mp.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it...";
        mp.computeTopology();
    }
    gsInfo<<"Finished\n";

    if (mp.geoDim()==2)
        mp.embed(3);


    // Boundary conditions
    gsBoundaryConditions<> bc;
    if (homogeneous)
    {
        for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        {
            bc.addCondition(*bit, condition_type::dirichlet, 0, 0, false, 0);
            bc.addCondition(*bit, condition_type::dirichlet, 0, 0, false, 1);
            bc.addCondition(*bit, condition_type::dirichlet, 0, 0, false, 2);
        }
    }
    else
    {
        fd.read(fn2);
        index_t num = 0;
        gsInfo<<"Reading BCs from "<<fn2<<"...";
        num = fd.template count<gsBoundaryConditions<>>();
        GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
        fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
        gsInfo<<"Finished\n";

    }
    bc.setGeoMap(mp);

    // Distributed load
    gsFunctionExpr<> force;
    gsInfo<<"Reading force function from "<<fn2<<" (ID=21) ...";
    fd.getId(21, force);
    gsInfo<<"Finished\n";

    // Point loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsMatrix<index_t> pid_ploads;
    gsInfo<<"Reading point load point locations from "<<fn2<<" (ID=30) ...";
    if ( fd.hasId(30) ) fd.getId(30,points);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load point vectors from "<<fn2<<" (ID=31) ...";
    if ( fd.hasId(31) ) fd.getId(31,loads);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load point patch indices from "<<fn2<<" (ID=32) ...";
    if ( fd.hasId(32) ) fd.getId(32,pid_ploads);
    gsInfo<<"Finished\n";

    if ( !fd.hasId(30) || !fd.hasId(31) || !fd.hasId(32) )
        pid_ploads = gsMatrix<index_t>::Zero(1,points.cols());

    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!

    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refPars, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations from "<<fn2<<" (ID=50) ...";
    if ( fd.hasId(50) ) fd.getId(50,refPoints);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches from "<<fn2<<" (ID=51) ...";
    if ( fd.hasId(51) ) fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference values from "<<fn2<<" (ID=52) ...";
    if ( fd.hasId(52) ) fd.getId(52,refValue);

    if ( !fd.hasId(50) || !fd.hasId(51) || !fd.hasId(52) )
        refValue = gsMatrix<>::Zero(mp.geoDim(),refPoints.cols());
    GISMO_ENSURE(refPatches.cols()==refPoints.cols(),"Number of reference points and patches do not match");

    if (refPoints.rows()==2)
    {
        refPars = refPoints;
        gsInfo<<"Reference points are provided in parametric coordinates.\n";
    }
    else if (refPoints.rows()==3)
        gsInfo<<"Reference points are provided in physical coordinates.\n";
    else
        gsInfo<<"No reference points are provided.\n";

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

    // Solver options
    gsInfo<<"Reading solver options from "<<fn1<<"...";
    gsOptionList solverOptions;
    fd.read(fn3);
    fd.template getFirst<gsOptionList>(solverOptions);
    gsInfo<<"Finished\n";
    //! [Read input file]

    //! [Create output directory]
    if ((plot || write) && !dirname.empty())
        gsFileManager::mkdir(dirname);
    //! [Create output directory]

    //! [Refine and elevate]
    gsMultiPatch<> geom;
    gsMultiBasis<> dbasis(mp);

    GISMO_ENSURE(degree>=mp.patch(0).degree(0),"Degree must be larger than or equal to the degree of the initial geometry, but degree = "<<degree<<" and the original degree = "<<mp.patch(0).degree(0));
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    if (method != -1 && method != 2)
        mp.degreeElevate(degree-mp.patch(0).degree(0));
    else
        dbasis.setDegree( degree); // preserve smoothness

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
    {
        if (method != -1 && method != 2)// && method != 3)
            mp.uniformRefine(1,degree-smoothness);
        else
            dbasis.uniformRefine(1,degree-smoothness);
    }
    if (plot) gsWriteParaview<>( mp    , dirname + "/mp", 1000, mesh);

    //! [Construct unstructured spline]
    // Mapped basis
    gsMappedBasis<2,real_t> bb2;
    // Transfer matrix
    gsSparseMatrix<> global2local;

    if (!fn4.empty())
    {
        gsInfo<<"Loading the basis data from file: "<<fn4<<"\n";
        fd.read(fn4); //filename: "square_knt"
        fd.getFirst(geom);
        fd.getFirst(dbasis);
        fd.getFirst(global2local);//mb.setTopology(mp);
        bb2.init(dbasis,global2local);
    }
    else
    {
        if (method==-1)
        {
            // identity map
            global2local.resize(dbasis.totalSize(),dbasis.totalSize());
            for (size_t k=0; k!=dbasis.totalSize(); ++k)
                global2local.coeffRef(k,k) = 1;
            geom = mp;
            gsInfo << "Basis Patch: " << dbasis.basis(0).component(0) << "\n";
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
            approxC1.options().setSwitch("interpolation",true);
            approxC1.options().setInt("gluingDataDegree",-1);
            approxC1.options().setInt("gluingDataSmoothness",-1);
            approxC1.update(bb2);
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
        }
        else if (method==4)
        {
            geom = mp;
            gsAlmostC1<2,real_t> almostC1(geom);
            almostC1.compute();
            almostC1.matrix_into(global2local);

            global2local = global2local.transpose();
            geom = almostC1.exportToPatches();
            dbasis = almostC1.localBasis();
            bb2.init(dbasis,global2local);
        }
        else
            GISMO_ERROR("Option "<<method<<" for method does not exist");

    }

    if (plot) gsWriteParaview(geom,dirname + "/geom",1000,mesh);

    //! [Construct unstructured spline]

    //! [Make assembler]
    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(geom,t,parameters,rho);

    // Construct the gsThinShellAssembler
    gsThinShellAssembler<3, real_t, true> assembler(geom,dbasis,bc,force,&materialMatrix);
    // if (method==1)
    assembler.options().setInt("Continuity",-1);
    assembler.options().setReal("WeakDirichlet",bcDirichlet);
    assembler.options().setReal("WeakClamped",bcClamped);
    assembler.setSpaceBasis(bb2);

    assembler.setPointLoads(pLoads);
    //! [Make assembler]

    //! [Define jacobian and residual]
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
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

        assembler.assembleVector(def);
        return assembler.rhs();
    };
    //! [Define jacobian and residual]

    //! [Assemble linear part]
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]

    //! [Solve non-linear problem]
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
    //! [Solve non-linear problem]

    //! [Construct and evaluate solution]
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
    //! [Construct and evaluate solution]

    //! [Export reference point data]
    gsField<> solField(geom, mspline,true);
    if (refPoints.cols()!=0)
    {
        std::stringstream buffer;
        std::ofstream file;
        gsMatrix<> result;
        if (refPoints.rows()==3) // then they are provided in the physical domain and should be mapped to the parametric domain
        {
            refPars.resize(2,refPoints.cols());
            for (index_t p = 0; p!=refPoints.cols(); p++)
            {
                geom.patch(refPatches(0,p)).invertPoints(refPoints.col(p),result,1e-10);
                if (result.at(0)==std::numeric_limits<real_t>::infinity()) // if failed
                    gsWarn<<"Point inversion failed\n";
                refPars.col(p) = result;
            }
        }

        gsInfo<<"Physical coordinates of points\n";
        for (index_t p=0; p!=refPars.cols(); p++)
            buffer<<"x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p)<<",";
        buffer<<"\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
        {
            geom.patch(refPatches(0,p)).eval_into(refPars.col(p),result);
            buffer<<result.row(0)<<","<<result.row(1)<<","<<result.row(2)<<",";
        }
        buffer<<"\n";
        if (write) writeToFile(buffer.str(),file,dirname + "/pointcoordinates.csv");
        buffer.str(std::string());

        gsMatrix<> refs(1,mp.geoDim()*refPoints.cols());
        for (index_t p=0; p!=refPars.cols(); p++)
            refs.block(0,p*mp.geoDim(),1,mp.geoDim()) = mp.piece(refPatches(0,p)).eval(refPars.col(p)).transpose();
        gsInfo<<"Reference point coordinates\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            buffer<<"x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p)<<",";
        buffer<<"\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            buffer<<refs(0,mp.geoDim()*p)<<","<<refs(0,mp.geoDim()*p+1)<<","<<refs(0,mp.geoDim()*p+2)<<",";
        buffer<<"\n";
        gsInfo<<buffer.str();
        if (write) writeToFile(buffer.str(),file,dirname + "/refpointcoordinates.csv");
        buffer.str(std::string());

        for (index_t p=0; p!=refPars.cols(); p++)
            refs.block(0,p*mp.geoDim(),1,mp.geoDim()) = solField.value(refPars.col(p),refPatches(0,p)).transpose();
        gsInfo<<"Computed values\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            buffer<<"x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p)<<",";
        buffer<<"\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            buffer<<refs(0,mp.geoDim()*p)<<","<<refs(0,mp.geoDim()*p+1)<<","<<refs(0,mp.geoDim()*p+2)<<",";
        buffer<<"\n";
        gsInfo<<buffer.str();
        if (write) writeToFile(buffer.str(),file,dirname + "/solution.csv");
        buffer.str(std::string());

        buffer<<"Reference values\n"; // provided as mp.geoDim() x points.cols() matrix
        for (index_t p=0; p!=refValue.cols(); ++p)
            buffer<<"x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p)<<",";
        buffer<<"\n";
        for (index_t p=0; p!=refValue.cols(); ++p)
            for (index_t d=0; d!=mp.geoDim(); d++)
                buffer<<refValue(d,p)<<",";
        buffer<<"\n";
        gsInfo<<buffer.str();
        if (write) writeToFile(buffer.str(),file,dirname + "/refsolution.csv");
        buffer.str(std::string());
    }
    //! [Export reference point data]


    //! [Export visualization in ParaView]
    if (plot)
    {
        // 4. Plot the mapped spline on the original geometry
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField,dirname + "/" + "Deformation", 1000, mesh);
    }

    if (stress)
    {
        gsMappedSpline<2,real_t> mspline_ori(bb2,solZero);
        gsMappedSpline<2,real_t> mspline_def(bb2,solFull);

        gsFunctionSum<real_t> ori(&geom,&mspline_ori);
        gsFunctionSum<real_t> def(&geom,&mspline_def);

        gsPiecewiseFunction<> membraneStresses;
        assembler.constructStress(ori,def,membraneStresses,stress_type::membrane);
        gsWriteParaview(ori,membraneStresses,dirname + "/" + "MembraneStress",5000);

        gsPiecewiseFunction<> membraneStressesVM;
        assembler.constructStress(ori,def,membraneStressesVM,stress_type::von_mises_membrane);
        gsWriteParaview(ori,membraneStressesVM,dirname + "/" + "MembraneStressVM",5000);

        gsPiecewiseFunction<> flexuralStresses;
        assembler.constructStress(ori,def,flexuralStresses,stress_type::flexural);
        gsWriteParaview(geom,flexuralStresses,dirname + "/" + "FlexuralStress",5000);
    }
    //! [Export visualization in ParaView]
    return EXIT_SUCCESS;
}
#else
int main(int argc, char *argv[])
{
    GISMO_ERROR("G+Smo is not compiled with the gsUnstructuredSplines module.");
    return EXIT_FAILURE;
}
#endif
