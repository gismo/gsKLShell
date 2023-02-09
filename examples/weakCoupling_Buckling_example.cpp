/** @file example_shell3D.cpp

    @brief Simple 3D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsThinShellUtils.h>

#include <gsSpectra/gsSpectra.h>

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool mesh       = false;
    bool first      = false;
    bool write      = false;
    bool info       = false;
    bool dense      = false;

    index_t numRefine  = 2;
    index_t degree = 3;
    index_t smoothness = 2;
    index_t nmodes = 10;
    index_t mode   = 0;

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    real_t ifcDirichlet = 1.0;
    real_t ifcClamped = 1.0;

    real_t shift = 0.01;

    std::string fn1,fn2,fn3;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";
    std::string out = "BucklingResults";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );
    cmd.addReal( "d", "DirIfc", "Dirichlet penalty scalar",  ifcDirichlet );
    cmd.addReal( "c", "ClaIfc", "Clamped penalty scalar",  ifcClamped );
    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addString( "o", "out", "Output directory",  out );
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addInt( "N", "nmodes", "Number of modes", nmodes );
    cmd.addInt( "M", "mode", "Mode number", mode );
    cmd.addReal( "S", "shift", "Set the shift of the solver.",  shift );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("first", "Plot only first mode",first);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("dense", "Dense eigenvalue computation",dense);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Define material parameters and geometry per example]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    gsBoundaryConditions<> bc;

    gsFileData<> fd;
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

    gsInfo<<"Patch 0 has basis: "<<mp.basis(0)<<"\n";

    gsInfo<<"Setting degree and refinement..."<<std::flush;
    GISMO_ENSURE(degree>=mp.patch(0).degree(0),"Degree must be larger than or equal to the degree of the initial geometry, but degree = "<<degree<<" and the original degree = "<<mp.patch(0).degree(0));
    mp.degreeIncrease(degree-mp.patch(0).degree(0));

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine(1,degree-smoothness);
    gsInfo<<"Finished\n";

    gsInfo<<"Patch 0 has basis: "<<mp.basis(0)<<"\n";

    if (plot)
    {
        std::string commands = "mkdir -p " + out;
        const char *command = commands.c_str();
        int systemRet = system(command);
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");
        gsWriteParaview(mp,out + "/" + "mp",10,true,false);
    }

    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    //     gsDebugVar(mp.patch(p));

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters,rho);

    gsThinShellAssembler<3, real_t, true> assembler;

    //! [Solver loop]
    gsVector<> solVector;

    gsMappedBasis<2,real_t> bb2;

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;

    gsInfo<<"Making gsMultiBasis..."<<std::flush;
    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Finished\n";
    // gsMappedSpline<2,real_t> mspline(bb2,coefs);

    assembler = gsThinShellAssembler<3, real_t, true>(mp,dbasis,bc,force,&materialMatrix);
    assembler.options().setReal("WeakDirichlet",bcDirichlet);
    assembler.options().setReal("WeakClamped",bcClamped);
    // Set the penalty parameter for the interface C1 continuity
    assembler.options().setInt("Continuity",-1);
    assembler.options().setReal("IfcDirichlet",ifcDirichlet);
    assembler.options().setReal("IfcClamped",ifcClamped);
    assembler.addWeakC0(mp.topology().interfaces());
    assembler.addWeakC1(mp.topology().interfaces());
    assembler.initInterfaces();
    // gsOptionList options = assembler.options();
    // options.setInt("Continuity",1);
    // assembler.setOptions(options);
    assembler.setPointLoads(pLoads);


    // Initialize the system

    gsInfo<<"Assembling stiffness matrix..."<<std::flush;
    assembler.assemble();
    gsSparseMatrix<> K_L = assembler.matrix();
    gsInfo<<"Finished\n";
    gsVector<> rhs = assembler.rhs();


    gsSparseSolver<>::CGDiagonal solver;
    solver.compute(K_L);
    solVector = solver.solve(rhs);

    gsMultiPatch<> deformation;
    deformation = assembler.constructDisplacement(solVector);
    gsField<> solField(mp,deformation);
    gsWriteParaview(solField,"LinearSolution");

    gsInfo<<"Assembling nonlinear stiffness matrix..."<<std::flush;
    assembler.constructSolution(solVector,mp_def);
    assembler.assembleMatrix(mp_def);
    gsSparseMatrix<real_t> dK = assembler.matrix();
    dK = dK - K_L;
    gsInfo<<"Finished\n";

    gsVector<> values;
    gsMatrix<> vectors;

    gsInfo<<"Computing Eigenmodes..."<<std::flush;
    if (dense)
    {
        Eigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  eigSolver;
        eigSolver.compute(K_L,dK);
        values = eigSolver.eigenvalues();
        vectors = eigSolver.eigenvectors();
    }
    else
    {
#ifdef GISMO_WITH_SPECTRA
        Spectra::SortRule selectionRule = Spectra::SortRule::LargestMagn;
        Spectra::SortRule sortRule = Spectra::SortRule::SmallestMagn;

        index_t ncvFac = 10;
        index_t number = nmodes;
        gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> eigSolver(K_L,dK,number,ncvFac*number, shift);
        // gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> eigSolver(matrix,mass,number,ncvFac*number, shift);
        eigSolver.init();
        eigSolver.compute(selectionRule,1000,1e-6,sortRule);

        if (eigSolver.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<eigSolver.num_iterations()<<" iterations and with "<<eigSolver.num_operations()<<"operations. \n"; }
        else if (eigSolver.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
        else if (eigSolver.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
        else if (eigSolver.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
        else                                                      { GISMO_ERROR("No error code known"); }

        values  = eigSolver.eigenvalues();
        values.array() += shift;
        vectors = eigSolver.eigenvectors();
#else
        Eigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  eigSolver;
        eigSolver.compute(K_L,dK);
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
        std::string output = "modes";
        gsParaviewCollection collection(out + "/" + output);

        int N = 1;
        if (!first)
            N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {
            solVector = vectors.col(m).normalized();
            deformation = assembler.constructDisplacement(solVector);

            // compute the deformation spline
            real_t maxAmpl = std::max(math::abs(deformation.patch(0).coefs().col(2).maxCoeff()),math::abs(deformation.patch(0).coefs().col(2).minCoeff()));
            for (size_t p=1; p!=mp.nPatches(); p++)
            {
                // deformation.patch(p).coefs() -= mp.patch(p).coefs();// assuming 1 patch here
                // Normalize mode shape amplitude in z coordinate
                maxAmpl = std::max(maxAmpl,std::max(math::abs(deformation.patch(p).coefs().col(2).maxCoeff()),math::abs(deformation.patch(p).coefs().col(2).minCoeff())));
            }
            for (size_t p=0; p!=mp.nPatches(); p++)
                if (maxAmpl!=0.0)
                    deformation.patch(p).coefs() = deformation.patch(p).coefs()/maxAmpl;

            gsField<> solField(mp,deformation);

            std::string fileName = out + "/" + output + util::to_string(m) + "_";
            gsWriteParaview<>(solField, fileName, 1000,mesh);
            for (index_t p = 0; p!=mp.nPatches(); p++)
            {
                fileName = output + util::to_string(m) + "_";
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

}// end main
