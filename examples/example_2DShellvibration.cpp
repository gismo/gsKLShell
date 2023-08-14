/** @file example_2DShellvibration.cpp

    @brief 2D shell vibration

    Example 26 from Liu et al 2006

    Liu, B., Xing, Y., Wang, Z., Lu, X., & Sun, H. (2017) Non-uniform rational Lagrange functions and its applications to isogeometric analysis of in-plane and flexural vibration of thin plates. Computer Methods in Applied Mechanics and Engineering, 321, 173-208. http://dx.doi.org/10.1016/j.cma.2017.04.007

    Our results:
    -----------------------------------------------------------
    Omega 0 = 13.264506012621019693
    Omega 1 = 25.948329127222432788
    Omega 2 = 46.868225986961668639
    Omega 3 = 68.842673494422726321
    Omega 4 = 74.844691337606633397

    Liu et al
    ------------------------------------------------------------
    Omega 0 = 7.98598284677940
    Omega 1 = 22.7690915637245
    Omega 2 = 47.5081920801538
    Omega 3 = 69.6137907904339
    Omega 4 = 80.9901698316950
    Author(s): J. Li
 **/

//! [Include namespace]
#include <gismo.h>

#include <iostream>
#include <fstream>
#include <gismo.h>
#include <gsKLShell/gsThinShellUtils.h>
#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

#include <gsStructuralAnalysis/gsModalSolver.h>

using namespace gismo;
//! [Include namespace]
void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
            std::string str = std::to_string(matrix(i,j));
            if(j+1 == matrix.cols()){
                file<<std::setprecision(20)<<matrix(i,j);
            }else{
                file<<std::setprecision(20)<<matrix(i,j)<<",";
            }
        }
        file<<'\n';
    }
}


int main(int argc, char *argv[])
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    int numKref     = 1;
    bool plot       = true;
    bool nonlinear  = false;
    bool first  = false;
    int mode = 1;

    real_t E_modulus = 71e9; //Pa
    real_t PoissonRatio = 0.3;
    real_t Density = 2700.0; //kg/m^3
    real_t thickness = 0.05; //mm


//    real_t rho = 1e0;

    int result = 0;
    bool write = true;
    bool plot_eigens = true;



    gsCmdLine cmd("2D shell vibration problem");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("m", "nmode",
               "Mode shape number, starting from 0",
               mode);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("first", "Plot only first", first);
    cmd.addSwitch("write", "Write convergence data to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Generate Geometry]
    real_t EA,EI,r,D;
    real_t PI = 3.1415926535;
    real_t r0 = 0.5;
    real_t r1 = 1.0;
    real_t Area = 0.25*PI*(pow(r1,2) - pow(r0,2));
    EI = 1.0/12.0*((r1-r0)*math::pow(thickness*2,3))*E_modulus;
    std::string fn;
    fn = "surfaces/quarterannulus.xml";
    EA = E_modulus*Area;
    r = math::sqrt(EI/EA);
    gsInfo<<"EI = "<<EI<<"; EA = "<<EA<<"; r = "<<r<<"\n";

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::NurbsQuarterAnnulus(r0,r1));
    mp.addAutoBoundaries();
    mp.embed(3);

    //! [Generate Geometry]

    for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";
    gsInfo<<"Coefficients (patch 0): "<< mp.patch(0).coefs().size()<< "\n";

    gsWriteParaview<>(mp, "mp", 500,true,false);

    // Initiate eigenfrequency
    real_t omega1 = 0.0;
    real_t poisson_ratio = 0.3;
    real_t poisson_ratio_new = poisson_ratio/(1-poisson_ratio);
    real_t E_modulus_new = E_modulus/(1-poisson_ratio * poisson_ratio);

    // Boundary conditions
    std::vector< std::pair<patchSide,int> > clamped;
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsPointLoads<real_t> pMass = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 );
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 1 );
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 2 );
    BCs.addCondition(boundary::west, condition_type::clamped,0,0,true,-1);
    // Right
//    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 );
//    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 1 );
//    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 2 );
//    BCs.addCondition(boundary::east, condition_type::clamped,0,0,true,-1);
    // Top
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0 );
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 1 );
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 2 );
    BCs.addCondition(boundary::north, condition_type::clamped,0,0,true,-1);
    // Bottom
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 );
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 1 );
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 2 );
    BCs.addCondition(boundary::south, condition_type::clamped,0,0,true,-1);




    omega1 = 6.8138;


    //! [Refinement]

    for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();

    for(index_t i = 0; i< numKref; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    gsMultiBasis<> dbasis(mp);


    // Initialize solution object

    gsFunctionExpr<> surfForce(tx,ty,tz,3);

    // Construct assembler object

    gsOptionList options;

    index_t material = 0;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    // Linear isotropic material model
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    real_t G_shear = E_modulus_new/2*(1+poisson_ratio_new);
    real_t D_shear = E_modulus * pow(thickness,3)/(12 * pow(1-poisson_ratio,2));
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(poisson_ratio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
//    gsConstantFunction<> ratio(Ratio,3);
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

    gsMaterialMatrixBase<real_t>* materialMatrix;


    options.addInt("Materal", "Material model: (0): SvK | (1): NH ", material);
    options.addInt("Implementation", "Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral", 1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp, t, parameters, rho, options);

    // Construct assembler object

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

    std::string assemberOptionsFile("options/solver_options.xml");
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    assembler->setOptions(opts);

    // Initialise solution and object
    gsMultiPatch<> solution = mp;

    // assemble system solve the system Kw + Mw'' = 0
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    assembler->assembleMass();
    gsSparseMatrix<> M = assembler->massMatrix();
    gsMatrix<real_t> Minv = M.toDense().inverse();

    gsDebug<<"K = \n"<<K.toDense()<<"\n";
    gsDebug<<"M = \n"<<M.toDense()<<"\n";
    gsDebug<<"Determinant of M = "<<M.toDense().determinant()<<"\n";
    gsDebug<<"Inverse of M = \n"<<M.toDense().inverse()<<"\n";

    gsModalSolver<real_t> modal(K,M);
    modal.options().setInt("solver",2);
    modal.options().setInt("selectionRule",0);
    modal.options().setInt("sortRule",4);
    modal.options().setSwitch("verbose",true);
    modal.options().setInt("ncvFac",2);

    modal.compute();

    gsMatrix<> values = modal.values();
    gsMatrix<> vectors = modal.vectors();

    gsMatrix<> dimensionless_value(values.rows(),values.cols());

    for (index_t k = 0; k<values.rows(); k++){
        dimensionless_value.at(k) = pow(values.at(k),0.5);
    }

    dimensionless_value = dimensionless_value * pow(Density * thickness/D_shear,0.5);

    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<dimensionless_value.at(k)<<"\n";
    gsInfo<<"\n";

    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        int systemRet = system("mkdir -p ModalResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        gsMultiPatch<> deformation = solution;
        gsMatrix<> modeShape;
        gsParaviewCollection collection("ModalResults/modes");

        int N = 1;
        if (!first)
            N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {

            // Compute solution based on eigenmode with number 'mode'
            modeShape = modal.vector(m);
            assembler->constructSolution(modeShape, solution);

            // compute the deformation spline
            deformation = solution;
            deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

            // Normalize mode shape amplitude in z coordinate
            real_t maxAmpl = std::max(math::abs(deformation.patch(0).coefs().col(2).maxCoeff()),math::abs(deformation.patch(0).coefs().col(2).minCoeff()));
            if (maxAmpl!=0.0)
            {
                deformation.patch(0).coefs() = deformation.patch(0).coefs()/maxAmpl;
            }

            gsField<> solField(mp,deformation);
            std::string fileName = "ModalResults/modes" + util::to_string(m);
            gsWriteParaview<>(solField, fileName, 5000);
            fileName = "modes" + util::to_string(m) + "0";
            collection.addTimestep(fileName,m,".vts");

        }

        gsFunctionExpr<> analytical("0","0","sin(3.1415926535*x)*sin(3.1415926535*y)",3);
        gsPiecewiseFunction<> func(analytical);
        gsField<> an(mp,func);
        gsWriteParaview(an,"analytical");
        collection.save();
    }

    if (write)
    {
        int systemRet = system("mkdir -p ModalResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        std::string wnM = "ModalResults/eigenvalues.txt";

        writeToCSVfile(wnM,dimensionless_value);
    }

    delete materialMatrix;
    delete assembler;
    return result;

}// end main
