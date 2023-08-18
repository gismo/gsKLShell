/** @file TODO

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

#include <gsSpectra/gsSpectra.h>

#include <gsUtils/gsQuasiInterpolate.h>


#include <gsAssembler/gsExprAssembler.h>

#include <gsKLShell/gsThinShellUtils.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot       = false;
    bool mesh       = false;
    bool first      = false;
    bool write      = false;
    bool info       = false;
    bool writeMatrix= false;
    bool writeGeo   = false;
    bool dense      = false;
    index_t numRefine  = 2;
    index_t degree = 3;
    index_t smoothness = 2;
    index_t geometry = 1;
    index_t method = 0;
    index_t nmodes = 10;
    index_t mode   = 0;
    std::string input;

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    real_t shift = 0.01;

    std::string fn1,fn2,fn3,fn4;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";
    std::string out = "ModalResults";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );
    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addString( "o", "out", "Output directory",  out );
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addInt( "m", "method", "Smoothing method to use", method );
    cmd.addInt( "N", "nmodes", "Number of modes", nmodes );
    cmd.addInt( "M", "mode", "Mode number", mode );
    cmd.addReal( "S", "shift", "Set the shift of the solver.",  shift );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("first", "Plot only first mode",first);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("writeGeo", "write geometry",writeGeo);
    cmd.addSwitch("writeMat", "Write projection matrix",writeMatrix);
    cmd.addSwitch( "info", "Print information", info );
    cmd.addSwitch("dense", "Dense eigenvalue computation",dense);
    cmd.addString( "","fullbasis", "File containing all basis details: 1) a multi-patch containing the geometry, 2) a multi-basis containing the local basis, and 3) a sparse matrix being a basis mapper",  fn4 );

    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    gsBoundaryConditions<> bc;

    GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
    GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");
    if (method==3)
        GISMO_ENSURE(smoothness>=1 || smoothness <= degree-2,"Exact C1 method only works for smoothness <= p-2, but smoothness="<<smoothness<<" and p-2="<<degree-2);
    if (method==2 || method==3)
        GISMO_ENSURE(degree > 2,"Degree must be larger than 2 for the approx and exact C1 methods, but it is "<<degree);

    /*
        to do:
        - remove hard-coded IDs from XML reader
     */

    gsFileData<> fd;
    if (!fn4.empty())
    {
        gsInfo<<"Reading geometry from "<<fn4<<"..."<<std::flush;
        gsReadFile<>(fn4, mp);
    }
    else 
    {
        gsInfo<<"Reading geometry from "<<fn1<<"..."<<std::flush;
        gsReadFile<>(fn1, mp);
        if (mp.geoDim()==2)
            mp.embed(3);

        if (method==2 && mp.patch(0).degree(0)<3)
        {
            gsWarn<<"Degree must be larger than 2 for the approx C1 method. Performing degree elevation on the geometry..."<<std::flush;
            mp.degreeIncrease(degree-mp.patch(0).degree(0));
            gsInfo<<"Finished.\n";
        }
    }
    if (mp.nInterfaces()==0 && mp.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it..."<<std::flush;
        mp.computeTopology();
    }
    gsInfo<<"Finished\n";

    fd.read(fn2);
    index_t num = 0;
    gsInfo<<"Reading BCs from "<<fn2<<"..."<<std::flush;
    num = fd.template count<gsBoundaryConditions<>>();
    GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
    fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
    gsInfo<<"Finished\n";

    bc.setGeoMap(mp);

    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    gsInfo<<"Reading thickness from "<<fn2<<" (ID=10) ..."<<std::flush;
    fd.getId(10,t);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Young's Modulus from "<<fn2<<" (ID=11) ..."<<std::flush;
    fd.getId(11,E);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Poisson ratio from "<<fn2<<" (ID=12) ..."<<std::flush;
    fd.getId(12,nu);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading density from "<<fn2<<" (ID=13) ..."<<std::flush;
    fd.getId(13,rho);
    gsInfo<<"Finished\n";

    gsMultiBasis<> dbasis;
    if (fn4.empty())
    {
        gsInfo<<"Making gsMultiBasis..."<<std::flush;
        dbasis = gsMultiBasis<> (mp);
        gsInfo<<"Finished\n";

        gsInfo<<"Setting degree and refinement..."<<std::flush;
        if (method != -1 && method != 2)// && method != 3)
            mp.degreeIncrease(degree-mp.patch(0).degree(0));
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
    }
    
    if (plot || write)
    {
        std::string commands = "mkdir -p " + out;
        const char *command = commands.c_str();
        int systemRet = system(command);
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");
    }
    if (plot)
        gsWriteParaview(mp,"mp",10,true,false);
    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    //     gsDebugVar(mp.patch(p));

    gsMappedBasis<2,real_t> bb2;

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;
    gsMultiPatch<> geom = mp;
    if (!fn4.empty())
    {
        gsInfo<<"Loading the basis data from file: "<<fn4<<"\n";
        fd.read(fn4); //filename: "square_knt"
        fd.getFirst(dbasis);
        fd.getFirst(global2local);//mb.setTopology(mp);
        bb2.init(dbasis,global2local);
    }
    else
    {
        gsInfo<<"Constructing Map..."<<std::flush;
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
            approxC1.options().setSwitch("info",info);
            // approxC1.options().setSwitch("plot",plot);
            approxC1.options().setSwitch("interpolation",true);
            approxC1.options().setInt("gluingDataDegree",-1);
            approxC1.options().setInt("gluingDataSmoothness",-1);
            approxC1.update(bb2);
        }
        else if (method==3) // Andrea
        {
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
            almostC1.options().setSwitch("SharpCorners",false);
            almostC1.compute();
            almostC1.matrix_into(global2local);

            global2local = global2local.transpose();
            geom = almostC1.exportToPatches();
            dbasis = almostC1.localBasis();
            bb2.init(dbasis,global2local);
        }
        else
            GISMO_ERROR("Option "<<method<<" for method does not exist");

        gsInfo<<"Finished\n";

    }

    if (writeMatrix)
    {
        gsWrite(global2local,"mat");
        //gsWrite(geom,"geom");
        //gsWrite(dbasis,"dbasis");
    }
    if (writeGeo)
    {
        gsWrite(geom,"geom");
    }
    if (plot)
        gsWriteParaview(geom,out + "/" + "geom",200,true);

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(geom,t,parameters,rho);

    gsThinShellAssembler<3, real_t, true> assembler;

    //! [Solver loop]
    gsVector<> solVector;

    gsFunctionExpr<> force("0","0","0",3);
    assembler = gsThinShellAssembler<3, real_t, true>(geom,dbasis,bc,force,&materialMatrix);
    // if (method==1)
    assembler.options().setInt("Continuity",-1);
    // else if (method==2)
    //     assembler.options().setInt("Continuity",-1);
    assembler.options().setReal("WeakDirichlet",bcDirichlet);
    assembler.options().setReal("WeakClamped",bcClamped);
    assembler.setSpaceBasis(bb2);
    // gsOptionList options = assembler.options();
    // options.setInt("Continuity",1);
    // assembler.setOptions(options);

    // Initialize the system

    gsInfo<<"Assembling stiffness matrix..."<<std::flush;
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsInfo<<"Finished\n";
    // gsDebugVar(matrix.toDense());
    gsVector<> vector = assembler.rhs();
    gsInfo<<"Assembling mass matrix..."<<std::flush;
    assembler.assembleMass();
    gsSparseMatrix<> mass   = assembler.massMatrix();
    gsInfo<<"Finished\n";
    // gsDebugVar(mass.toDense());

    gsVector<> values;
    gsMatrix<> vectors;

    gsInfo<<"Computing Eigenmodes..."<<std::flush;
    if (dense)
    {
        Eigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  eigSolver;
        eigSolver.compute(matrix-shift*mass,mass);
        values = eigSolver.eigenvalues();
        vectors = eigSolver.eigenvectors();
    }
    else
    {
// #ifdef GISMO_WITH_SPECTRA
#ifdef false
        Spectra::SortRule selectionRule = Spectra::SortRule::LargestMagn;
        Spectra::SortRule sortRule = Spectra::SortRule::SmallestMagn;

        index_t ncvFac = 10;
        index_t number = nmodes;
        // gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> solver(matrix-shift*mass,mass,number,ncvFac*number, shift);
        gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> solver(matrix,mass,number,ncvFac*number, shift);
        solver.init();
        solver.compute(selectionRule,1000,1e-12,sortRule);

        if (solver.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<solver.num_iterations()<<" iterations and with "<<solver.num_operations()<<"operations. \n"; }
        else if (solver.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
        else if (solver.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
        else if (solver.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
        else                                                      { GISMO_ERROR("No error code known"); }

        values  = solver.eigenvalues();
        values.array() += shift;
        vectors = solver.eigenvectors();
#else
        Eigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  eigSolver;
        eigSolver.compute(matrix-shift*mass,mass);
        values = eigSolver.eigenvalues();
        vectors = eigSolver.eigenvectors();
#endif
    }

    gsInfo<<"Finished\n";


    gsDebugVar(values);

    // ! [Solver loop]


    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsMatrix<> modeShape;
        std::string dirname = out;
        std::string output = "modes";
        gsParaviewCollection collection(dirname + "/" + output);

        int N = 1;
        if (!first)
	      N = nmodes;
        //    N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {
            solVector = vectors.col(m).normalized();

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

            gsField<> solField(geom, mspline,true);

            std::string fileName = dirname + "/" + output + util::to_string(m);
            gsWriteParaview<>(solField, fileName, 1000,mesh);
            for (index_t p = 0; p!=geom.nPatches(); p++)
            {
                fileName = output + util::to_string(m);
                collection.addTimestep(fileName,p,m,".vts");
                if (mesh)
                    collection.addTimestep(fileName,p,m,"_mesh.vtp");
            }
        }
        collection.save();
    }
    if (write)
    {
        std::ofstream file;
        file.open(out + "/" + "eigenvalues.csv",std::ofstream::out);
        for (index_t k=0; k!=values.size(); k++)
            file<<std::setprecision(12)<<values.at(k)<<"\n";

        file.close();
    }

    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}
