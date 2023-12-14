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

#include <gsSpectra/gsSpectra.h>

#include <gsUtils/gsQuasiInterpolate.h>


#include <gsAssembler/gsExprAssembler.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

#include <gsKLShell/gsThinShellUtils.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot       = false;
    bool plotGeo    = false;
    bool mesh       = false;
    bool first      = false;
    bool write      = false;
    bool dense      = false;
    index_t method = 0;
    index_t nmodes = 10;
    index_t mode   = 0;
    std::string input;

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    real_t shift = 0.01;

    std::string geomFileName,basisFileName,bvpFileName,optFileName;
    optFileName = "options/solver_options.xml";
    std::string out = "ModalResults";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );
    cmd.addString( "G", "geom","File containing the geometry",  geomFileName );
    cmd.addString( "b", "bas", "File containing the basis (dbasis and global2local)",  basisFileName );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  bvpFileName );
    cmd.addString( "O", "opt", "File containing solver options",  optFileName );
    cmd.addString( "o", "out", "Output directory",  out );
    cmd.addInt( "N", "nmodes", "Number of modes", nmodes );
    cmd.addInt( "M", "mode", "Mode number", mode );
    cmd.addReal( "S", "shift", "Set the shift of the solver.",  shift );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("plotGeo", "plotGeo",plotGeo);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("first", "Plot only first mode",first);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("dense", "Dense eigenvalue computation",dense);

    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(!(geomFileName.empty()) && !(basisFileName.empty()) && !(bvpFileName.empty()),"Not all filenames have been provided.\n"
                        <<"geomFileName  = "<<geomFileName<<"\n"
                        <<"basisFileName = "<<basisFileName<<"\n"
                        <<"bvpFileName   = "<<bvpFileName);

    gsFileData<> fd;

    gsMultiPatch<> geom;
    gsSparseMatrix<> global2local;
    gsMultiBasis<> dbasis;

    gsMappedBasis<2,real_t> bb2;
    gsBoundaryConditions<> bc;

    gsInfo<<"Reading geometry from "<<geomFileName<<"..."<<std::flush;
    gsReadFile<>(geomFileName, geom);
    gsInfo<<"Finished\n";

    // STEP 1: Get curve network with merged linear interfaces
    gsInfo<<"Loading curve network..."<<std::flush;
    geom.computeTopology();
    geom.constructInterfaceRep();
    geom.constructBoundaryRep();
    auto & irep = geom.interfaceRep();
    auto & brep = geom.boundaryRep();
    // gsDebug <<" irep "<< irep.size() <<" \n" ;
    gsDebug <<" brep "<< brep.size() <<" \n" ;

    // outputing...
    gsMultiPatch<> crv_net, iface_net, bnd_net;
    for (auto it = irep.begin(); it!=irep.end(); ++it)
    {
        iface_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }
    for (auto it = brep.begin(); it!=brep.end(); ++it)
    {
        bnd_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }

    if (plot) gsWriteParaview(iface_net,"iface_net",100);
    if (plot) gsWriteParaview(bnd_net,"bnd_net",100);
    // if (plot) gsWriteParaview(crv_net,"crv_net",100);

    gsInfo<<"Reading mapped basis from "<<basisFileName<<"..."<<std::flush;
    fd.read(basisFileName);
    GISMO_ENSURE(fd.template getFirst<gsMultiBasis<>>(dbasis)        ,"No multibasis was read");
    GISMO_ENSURE(fd.template getFirst<gsSparseMatrix<>>(global2local),"No sparse matrix was read");
    gsInfo<<"Finished\n";
    bb2.init(dbasis,global2local);

    gsInfo<<"Reading BVP from "<<basisFileName<<"..."<<std::flush;
    fd.read(bvpFileName);
    index_t num = 0;
    gsInfo<<"\tReading BCs from "<<bvpFileName<<"...\n";
    num = fd.template count<gsBoundaryConditions<>>();
    GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
    fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
    gsInfo<<"Finished\n";

    bc.setGeoMap(geom);

    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    gsInfo<<"\tReading thickness from "<<bvpFileName<<" (ID=10) ..."<<std::flush;
    fd.getId(10,t);
    gsInfo<<"Finished\n";

    gsInfo<<"\tReading Young's Modulus from "<<bvpFileName<<" (ID=11) ..."<<std::flush;
    fd.getId(11,E);
    gsInfo<<"Finished\n";

    gsInfo<<"\tReading Poisson ratio from "<<bvpFileName<<" (ID=12) ..."<<std::flush;
    fd.getId(12,nu);
    gsInfo<<"Finished\n";

    gsInfo<<"\tReading density from "<<bvpFileName<<" (ID=13) ..."<<std::flush;
    fd.getId(13,rho);
    gsInfo<<"Finished\n";
    gsInfo<<"Finished\n";

    if (plot || write)
    {
        std::string commands = "mkdir -p " + out;
        const char *command = commands.c_str();
        int systemRet = system(command);
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");
    }

    if (plotGeo) gsWriteParaview(geom,"geom",1000,true,false);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(geom,t,parameters,rho);

    gsThinShellAssembler<3, real_t, true> assembler;

    //! [Solver loop]
    gsVector<> solVector;

    gsFunctionExpr<> force("0","0","0",3);
    assembler = gsThinShellAssembler<3, real_t, true>(geom,dbasis,bc,force,&materialMatrix);
    assembler.options().setInt("Continuity",-1);
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
#ifdef gsSpectra_ENABLED
        Spectra::SortRule selectionRule = Spectra::SortRule::LargestMagn;
        Spectra::SortRule sortRule = Spectra::SortRule::SmallestMagn;

        index_t ncvFac = 10;
        index_t number = nmodes;
        gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> solver(matrix-shift*mass,mass,number,ncvFac*number, shift);
        // gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> solver(matrix,mass,number,ncvFac*number, shift);
        solver.init();
        solver.compute(selectionRule,1000,1e-6,sortRule);

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
            gsFunctionSum<real_t> def(&geom,&mspline);

            gsField<> solField(geom, mspline,true);

            std::string fileName = dirname + "/" + output + util::to_string(m) + "_";
            gsWriteParaview<>(solField, fileName, 1000,mesh);
            for (index_t p = 0; p!=geom.nPatches(); p++)
            {
                fileName = output + util::to_string(m) + "_" + std::to_string(p);
                collection.addPart(fileName + ".vts",m,"solution",p);
                if (mesh)
                    collection.addPart(fileName + "_mesh.vtp",m,"mesh",p);
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
