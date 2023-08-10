/** @file example_3DBuckling.cpp
    @brief Axially compressed cylinder

    Fig 4 from Oesterle et al 2022

    Oesterle, B., Geiger, F., Forster, D., Fröhlich, M., & Bischoff, M. (2022). A study on the approximation power of NURBS and the significance of exact geometry in isogeometric pre-buckling analyses of shells. Computer Methods in Applied Mechanics and Engineering, 397. https://doi.org/10.1016/j.cma.2022.115144

    Author(s): J. Li
 **/
//! [Include namespace]
#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsStructuralAnalysis/gsBucklingSolver.h>

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
    int numElevate  = 1;
    int numHref     = 1;
    int numElevateL = -1;
    int numHrefL    = -1;
    bool plot       = true;
    bool sparse     = false;
    bool nonlinear  = false;
    bool first  = false;
    bool composite = false;
    int mode = 0;

    real_t E_modulus     = 2.1e5;
    real_t PoissonRatio = 0.0;
//    real_t Density = 1.0;

    real_t thickness = 0.1; // (m) slenderness ratio is 200
    real_t aDim = 20.0;
    real_t bDim = 30.0;

    index_t Compressibility = 0;
    index_t material = 0;
//    real_t Ratio = 7.0;

    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t fac = 1;

    real_t shift = 0.0;

    int testCase = 0;

    int result = 0;

    bool write = false;

    bool MIP = false;

    index_t nmodes = 1;

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

    cmd.addInt("m", "nmode",
               "Mode shape number, starting from 0",
               mode);

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
    gsReadFile<>("surfaces/cylinder.xml", mp);

    gsInfo<< mp.patch(0).coefs().row(0)(2);
    mp.embed(3);
    mp.patch(0).coefs().row(0) << 0.00, -aDim, 0.00;
    mp.patch(0).coefs().row(1) << aDim, -aDim, 0.00;
    mp.patch(0).coefs().row(2) << aDim, 0.00, 0.00;
    mp.patch(0).coefs().row(3) << aDim, aDim, 0.00;
    mp.patch(0).coefs().row(4) << 0.00, aDim, 0.00;
    mp.patch(0).coefs().row(5) << -aDim,aDim, 0.00;
    mp.patch(0).coefs().row(6) << -aDim, 0.00, 0.00;
    mp.patch(0).coefs().row(7) << -aDim, -aDim, 0.00;
    mp.patch(0).coefs().row(8) << 0.00, -aDim, 0.00;
    mp.patch(0).coefs().row(9) << 0.00, -aDim, bDim;
    mp.patch(0).coefs().row(10) << aDim, -aDim, bDim;
    mp.patch(0).coefs().row(11) << aDim, 0.00, bDim;
    mp.patch(0).coefs().row(12) << aDim, aDim, bDim;
    mp.patch(0).coefs().row(13) << 0.00, aDim, bDim;
    mp.patch(0).coefs().row(14) << -aDim, aDim, bDim;
    mp.patch(0).coefs().row(15) << -aDim, 0.00, bDim;
    mp.patch(0).coefs().row(16) << -aDim, -aDim, bDim;
    mp.patch(0).coefs().row(17) << 0.00, -aDim, bDim;

    mp.addInterface(0,boundary::west ,0,boundary::east);
    mp.computeTopology();
    mp.addAutoBoundaries();

    // Elevate the degree in vertical direction
    mp.patch(0).degreeElevate(2,1);
    mp.patch(0).degreeElevate(1,0);



    gsWriteParaview(mp, "mp", 1000, true);



    gsDebugVar(mp);


    if (numHrefL==-1)
        numHrefL = numHref;
    if (numElevateL==-1)
        numElevateL = numElevate;

//    real_t length,width;
    gsMultiBasis<double> dbasis(mp);

    gsInfo<< "Basis (patch 0): "<<mp.patch(0).basis()<<"\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");


    // Boundary conditions
    real_t Load = -1e4;
    gsVector<> neu(3);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    neu << 0, 0, Load;
    gsConstantFunction<> neuData(neu,3);

    //Buckling coefficient
    real_t pressure = 0.0;


    // Apply load to boundary
    // // Pinned-Pinned
    // BCs.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x


    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,0 ); // unknown 1 - y
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,1 ); // unknown 2 - z
    BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false,2 ); // unknown 0 - x
    BCs.addCondition(boundary::north, condition_type::neumann, &neuData); // unknown 0 - x

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

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

    // Construct assembler object
    assembler ->setOptions(opts);

    // Set the penalty parameter for the interface C1 continuity

    real_t ifcDirichlet = 1e3;
    real_t ifContinuity = -1;

    assembler->options().setInt("Continuity",ifContinuity);
    assembler->options().setReal("IfcDirichlet",ifcDirichlet);
    assembler->addWeakC0(mp.topology().interfaces());
    assembler->addWeakC1(mp.topology().interfaces());
    assembler->initInterfaces();
    // Initialise solution object
    gsMultiPatch<> solution = mp;

    assembler->assemble();
    gsSparseMatrix<> K_L =  assembler->matrix();
    gsVector<> rhs = assembler->rhs();

    gsDebugVar(rhs.transpose());


    gsSparseSolver<real_t>::CGDiagonal solver;
    solver.compute(K_L);
    gsMatrix<> sol = solver.solve(rhs);
    assembler->constructDisplacement(sol,mp_def);
    gsField<> solFieldTmp(mp,mp_def);
    gsWriteParaview<>(solFieldTmp,"tmp",1000);
    gsDebugVar(sol.transpose());

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                            Jacobian_t;
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &, gsVector<real_t> const &)>  dJacobian_t;
    Jacobian_t K_NL = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
        assembler->constructSolution(x,mp_def);
        assembler->assemble(mp_def);
        gsSparseMatrix<real_t> m = assembler->matrix();
        return m;
    };
    dJacobian_t dK_NL = [&assembler,&mp_def,&MIP](gsVector<real_t> const &x, gsVector<real_t> const &dx)
    {
        if (MIP)
            assembler->assembleMatrix(x,x-dx);
        else
        {
            assembler->constructSolution(x,mp_def);
            assembler->assembleMatrix(mp_def);
        }

        gsSparseMatrix<real_t> m = assembler->matrix();
        return m;
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
        buckling.computeSparse(shift,nmodes);//,2,Spectra::SortRule::LargestMagn,Spectra::SortRule::SmallestMagn);

    gsMatrix<> values = buckling.values();
    gsMatrix<> vectors = buckling.vectors();

    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    for (index_t k = 0; k<10; k++)
    {
        gsInfo<<"\t"<<values.at(k)*Load<<"\n";
    }
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        int systemRet = system("mkdir -p CylinderBucklingResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        gsMultiPatch<> deformation = solution;
        gsMatrix<> modeShape;
        gsParaviewCollection collection("CylinderBucklingResults/modes");

        int N = 1;
        if (!first)
            N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {

            // Compute solution based on eigenmode with number 'mode'
            modeShape = vectors.col(m);
            assembler->constructSolution(modeShape, solution);

            // compute the deformation spline
            deformation = solution;
            deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

            real_t maxAmpl = std::max(math::abs(deformation.patch(0).coefs().col(2).maxCoeff()),math::abs(deformation.patch(0).coefs().col(2).minCoeff()));
            if (maxAmpl!=0.0)
            {
                deformation.patch(0).coefs() = deformation.patch(0).coefs()/maxAmpl;
            }

            gsField<> solField(mp,deformation);
            std::string fileName = "CylinderBucklingResults/modes" + util::to_string(m);
            gsWriteParaview<>(solField, fileName, 5000);
            fileName = "modes" + util::to_string(m) + "0";
            collection.addTimestep(fileName,m,".vts");
        }
        collection.save();
    }

    if (write)
    {
        int systemRet = system("mkdir -p CylinderBucklingResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");
        std::string wnM = "CylinderBucklingResults/eigenvalues.txt";
        writeToCSVfile(wnM,values);
    }

    delete materialMatrix;
    delete assembler;
    return result;

}








