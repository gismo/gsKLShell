/** @file example_3DBuckling.cpp
    @brief Axially compressed cylinder

    Fig 4 from Oesterle et al 2022

    Oesterle, B., Geiger, F., Forster, D., Fröhlich, M., & Bischoff, M. (2022). A study on the approximation power of NURBS and the significance of exact geometry in isogeometric pre-buckling analyses of shells. Computer Methods in Applied Mechanics and Engineering, 397. https://doi.org/10.1016/j.cma.2022.115144

    Author(s): J. Li

    todo: In the paper they used 60 segments geometry, our geometry only has four patches.
 **/
//! [Include namespace]
#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsStructuralAnalysis/gsBucklingSolver.h>

#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsKLShell/gsFunctionSum.h>

using namespace gismo;

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
            std::string str = std::to_string(matrix(i,j));
            if(j+1 == matrix.cols()){
                file<<std::setprecision(10)<<str;
            }else{
                file<<std::setprecision(10)<<str<<',';
            }
        }
        file<<'\n';
    }
}

int main(int argc, char** argv){

    // Input options
    int numElevate  = 0;
    int numHref     = 4;
    int numElevateL = -1;
    int numHrefL    = -1;
    bool plot       = true;
    bool sparse     = false;
    bool nonlinear  = false;
    bool first  = false;
    bool composite = false;
    int mode = 0;
    int method = 0;

    bool membrane = false;

    real_t E_modulus     = 2.1e5;
    real_t PoissonRatio = 0.0;
//    real_t Density = 1.0;

    real_t thickness = 0.1; // (m) slenderness ratio is 200
    real_t aDim = 20.0;
    real_t bDim = 30.0;
    real_t Load = -1e4;

    index_t Compressibility = 0;
    index_t material = 0;
//    real_t Ratio = 7.0;

    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t fac = 1;

    real_t shift = 0.0;

    int testCase = 0;

    int result = 0;

    bool write = true;

    bool MIP = false;

    index_t nmodes = 1;

    std::string dirname = "Multipatch_cylinder";

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Buckling analysis for thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt( "N", "nmodes", "Number of modes",  nmodes );


    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("R","hRefine2", "Number of dyadic h-refinement (bisection) steps to perform before solving (secondary direction)", numHrefL);
    cmd.addInt("E","degreeElevation2", "Number of degree elevation steps to perform on the Geometry's basis before solving (secondary direction)", numElevateL);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);
    cmd.addInt("m","method", "Smoothing method to use: 0: smoothInterfaces, 1: Almost C1, 2: D-Patch", method);

    cmd.addSwitch("membrane", "Membrane element", membrane);
    cmd.addReal("T","hdim", "thickness of the plate", thickness);
    cmd.addReal("a","adim", "dimension a", aDim);
    cmd.addReal("b","bdim", "dimension b", bDim);

    cmd.addReal("F","fac", "factor linear problem", fac);

    cmd.addReal("s","shift", "eigenvalue shift", shift);

    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("first", "Plot only first", first);
    cmd.addSwitch("write", "Write convergence data to file", write);
    cmd.addSwitch("sparse", "Use sparse solver", sparse);
    cmd.addSwitch("MIP", "Use mixed integration point method", MIP);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Geometry Setup]
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);


    gsMultiPatch<> mp;
    gsFileData<real_t> fn("surfaces/cylinder_4p.xml");
    gsInfo<<"Reading geometry (ID=0) ...";
    fn.getId(0,mp);
    gsDebugVar(mp);
    gsInfo<<"Finished\n";


    // p-refine
    for (size_t p=0; p!=mp.nPatches(); p++)
    {
        for(index_t i = 0; i< numElevate; ++i)
        {
            if (dynamic_cast<gsTensorNurbs<2,real_t> * >(&mp.patch(p)))
            {
                gsWarn<<"Degree elevation applied"<<"\n";
                mp.patch(p).degreeElevate();    // Elevate the degree
            }
            else
                mp.patch(p).degreeIncrease();    // Elevate the degree
        }

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
            mp.patch(p).uniformRefine();
    }

    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";


    // Boundary conditions
    gsBoundaryConditions<> BCs;
    gsInfo<<"Reading boundary conditions (ID=20) ...";
    fn.getId(20,BCs);
    gsInfo<<"Finished\n";
    gsDebugVar(BCs);

    gsFunctionExpr<> forceFun;
    gsFunctionExpr<> pressFun;
//
//    gsInfo<<"Reading force function (ID=21) ...";
//    fd.getId(21,forceFun);
//    gsInfo<<"Finished\n";
    bool pressure = false;
    if ( fn.hasId(22) )
    {
        gsInfo<<"Reading pressure function (ID=22) ...";
        pressure = true;
        fn.getId(22,pressFun);
        gsInfo<<"Finished\n";
    }

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsMatrix<index_t> pid_ploads;
    // Point loads
    gsInfo<<"Reading point load locations (ID=30) ...";
    if ( fn.hasId(30) ) fn.getId(30,points);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load vectors (ID=31) ...";
    if ( fn.hasId(31) ) fn.getId(31,loads);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load patch indices (ID=32) ...";
    if ( fn.hasId(32) ) fn.getId(32,pid_ploads);
    gsInfo<<"Finished\n";


    //Material Matrix


    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");




    gsVector<> tmp(3);
    tmp << 0, 0, 0;



    gsFunctionExpr<> surfForce(tx,ty,tz,3);


    // Initialize solution object
    gsMultiPatch<> mp_def = mp;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp, 3);
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    // Anisotropic material, set the layers to one
    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(), 1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,3);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = &thicks;

    std::vector<gsFunction<>*> parameters; //SvK & Composites & NH & NH_ext
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;
    if      (material==0 && impl==1)
    {
        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            parameters.resize(2);
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);
    }

    gsDebugVar(dbasis.totalSize());

    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refPars, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations (ID=50) ...";
    if ( fd.hasId(50) ) fd.getId(50,refPoints);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches (ID=51) ...";
    if ( fd.hasId(51) ) fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference values (ID=52) ...";
    if ( fd.hasId(52) ) fd.getId(52,refValue);
    gsInfo<<"Finished\n";

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

    // Fix path
    dirname = gsFileManager::getCanonicRepresentation(dirname,true);
    char sep = gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(dirname);

    // plot geometry
    if (plot)
        gsWriteParaview(mp,"mp",1000,true);
    mp.computeTopology();
    gsDebugVar(mp);

    //! [Make unstructured spline]
    gsMultiPatch<> geom;
    gsMappedBasis<2,real_t> bb2;
    gsSparseMatrix<> global2local;

    if (method==0)
    {
        gsSmoothInterfaces<2,real_t> smoothInterfaces(mp);
        smoothInterfaces.options().setSwitch("SharpCorners",false);
        smoothInterfaces.compute();
        smoothInterfaces.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = smoothInterfaces.exportToPatches();
        dbasis = smoothInterfaces.localBasis();
    }
    else if (method==1)
    {
        gsAlmostC1<2,real_t> almostC1(mp);
        almostC1.options().setSwitch("SharpCorners",true);
        almostC1.compute();
        almostC1.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = almostC1.exportToPatches();
        dbasis = almostC1.localBasis();
    }
    else if (method==2)
    {
        gsDPatch<2,real_t> dpatch(mp);
        dpatch.options().setSwitch("SharpCorners",true);
        dpatch.compute();
        dpatch.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = dpatch.exportToPatches();
        dbasis = dpatch.localBasis();
    }
    else
        GISMO_ERROR("Method "<<method<<" unknown");

    bb2.init(dbasis,global2local);
    // Finalize BCs
    BCs.setGeoMap(geom);


    gsDebugVar(geom);

    //Probes
    gsMatrix<> writePoints;
    gsMatrix<index_t> writePatches;

    //Cross-section
    std::vector<gsMatrix<>> writeSection;
    std::vector<index_t> sectionPatches;

    writePoints.resize(2,1);
    writePoints.col(0)<< 1.0,1.0;

    writePatches.resize(1,1);
    writePatches.row(0)<<0;



    gsThinShellAssemblerBase<real_t>* assembler;

    if (membrane)
        assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

    //Construct assembler object
    assembler->setOptions(opts);
    if(method == 0) assembler->options().setInt("Continuity",-1);
    else assembler->options().setInt("Continuity",0);

    if(method == 0)
        assembler->setSpaceBasis(bb2);

    assembler->setPointLoads(pLoads);
    if (pressure)
        assembler->setPressure(pressFun);

    // Assemble linear system to obtain the force vector
    gsDebugVar(assembler->numDofs());
    assembler->assemble();
    gsSparseMatrix<> K_L =  assembler->matrix();
    gsVector<> rhs = assembler->rhs();
    gsDebugVar(rhs.transpose());
//    assembler->assembleMass(true);
//    gsVector<> M = assembler->rhs();
    mp_def = geom;





//    gsSparseSolver<real_t>::CGDiagonal solver;
//    solver.compute(K_L);
//    gsMatrix<> sol = solver.solve(rhs);
//    assembler->constructDisplacement(sol,mp_def);
//    gsField<> solFieldTmp(mp,mp_def);
//    gsWriteParaview<>(solFieldTmp,"tmp",1000);
//    gsDebugVar(sol.transpose());


    gsStructuralAnalysisOps<real_t>::Jacobian_t K_NL = [&assembler,&bb2,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        ThinShellAssemblerStatus status;
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
        size_t d = mp_def.targetDim();
        GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/d,d);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&mp_def,&mspline);
        assembler->constructSolution(solFull,mp_def);
        assembler->assembleMatrix(def);
        status = assembler->assembleMatrix(def);
        m = assembler->matrix();
//        assembler->constructSolution(x,mp_def);
//        status = assembler->assembleMatrix(mp_def);
//        m = assembler->matrix();
        return status == ThinShellAssemblerStatus::Success;
    };
    gsStructuralAnalysisOps<real_t>::dJacobian_t dK_NL = [&assembler,&bb2,&mp_def,&MIP](gsVector<real_t> const &x, gsVector<real_t> const &dx, gsSparseMatrix<real_t> & m)
    {
        ThinShellAssemblerStatus status;
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
        size_t d = mp_def.targetDim();
        GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/d,d);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        if (MIP)
            status = assembler->assembleMatrix(x,solFull-dx);
        else
        {
            gsFunctionSum<real_t> def(&mp_def,&mspline);
//            assembler->constructSolution(x,mp_def);
            status = assembler->assembleMatrix(def);
        }
        m = assembler->matrix();
        return status == ThinShellAssemblerStatus::Success;
    };

    gsBucklingSolver<real_t> buckling(K_L,rhs,dK_NL);
    buckling.options().setInt("solver",2);
    buckling.options().setInt("selectionRule",0);
    buckling.options().setInt("sortRule",4);
    buckling.options().setSwitch("verbose",true);
    buckling.options().setInt("ncvFac",2);


//    buckling.computePower();

    if (!sparse)
        buckling.compute();
    else
        buckling.computeSparse(nmodes);//,2,Spectra::SortRule::LargestMagn,Spectra::SortRule::SmallestMagn);

    gsMatrix<> values = buckling.values();
    gsMatrix<> vectors = buckling.vectors();

    gsMultiPatch<> solution = geom;

    bool mesh = false;

    gsMatrix<> dimensionless_value(values.rows(),values.cols());

    for (index_t k = 0; k<values.rows(); k++){
        dimensionless_value.at(k) = pow(values.at(k),0.5);
    }


    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    for (index_t k = 0; k<10; k++)
    {
        gsInfo<<"\t"<<values.at(k)*Load<<"\n";
    }
    std::string output = "solution";

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        int systemRet = system("mkdir -p Multipatch_cylinder");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        gsMultiPatch<> deformation = solution;
        gsMatrix<> modeShape;
        gsParaviewCollection collection("Multipatch_cylinder/solutions");

        int N = 1;

        bool first = false;


        if (!first)
            N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {
            gsField<> solField;

            std::string fileName = dirname + "/" + output + util::to_string(m);

            /// Make a gsMappedSpline to represent the solution
            // 1. Get all the coefficients (including the ones from the eliminated BCs.)
            gsMatrix<real_t> solFull = assembler->fullSolutionVector(vectors.col(m));

            // 2. Reshape all the coefficients to a Nx3 matrix
            size_t d = mp_def.targetDim();
            GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
            solFull.resize(solFull.rows()/d,d);

            // 3. Make the mapped spline
            gsMappedSpline<2,real_t> mspline(bb2,solFull);

            // 4. Plot the mapped spline on the original geometry
            solField = gsField<>(geom, mspline,true);

            gsWriteParaview<>(solField, fileName, 1000,false,"_");
            // Compute solution based on eigenmode with number 'mode'
            modeShape = vectors.col(m);
            assembler->constructSolution(modeShape, solution);

            // compute the deformation spline
//            deformation = solution;
//            for(int k = 0; k < deformation.nPatches();k++)
//                deformation.patch(k).coefs() -= geom.patch(k).coefs();

            // gsField<> mpField(mp,mp);
            //
            // real_t norm = solField.distanceL2(mpField);

            // Normalize mode shape amplitude in z coordinate
//            gsVector<> maxAmpl(4);
//            for(int k = 0; k < deformation.nPatches();k++) {
//                maxAmpl.row(k) << std::max(math::abs(deformation.patch(k).coefs().col(2).maxCoeff()),
//                                           math::abs(deformation.patch(k).coefs().col(2).minCoeff()));
//                if (maxAmpl.row(k)[0] != 0.0) {
//                    deformation.patch(k).coefs() = deformation.patch(k).coefs() / maxAmpl.row(k)[0];
//                }
//            }
            for (size_t p = 0; p!=geom.nPatches(); p++)
            {
                fileName = output + util::to_string(m) + "_" + util::to_string(p);
                collection.addPart(fileName + ".vts",m);
                if (mesh) collection.addPart(fileName + "_mesh.vtp",m);
            }
        }
        collection.save();
    }

    if (write)
    {
        int systemRet = system("mkdir -p Multipatch_cylinder");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");
        std::string wnM = "Multipatch_cylinder/eigenvalues.txt";
        writeToCSVfile(wnM,values);
    }

    delete materialMatrix;
    delete assembler;

    return result;

}



