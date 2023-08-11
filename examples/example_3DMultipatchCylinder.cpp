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

    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool nonlinear = false;

    bool membrane = false;
    bool composite = false;

    real_t E_modulus     = 2.1e5;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 0.01;

    real_t aDim = 20.0;
    real_t bDim = 30.0;

    real_t ifcDirichlet = 1.0;
    real_t ifcClamped = 1.0;

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Buckling analysis for multipatch thin shell.");
    cmd.addReal( "D", "Dir", "Dirichlet penalty scalar",  ifcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped penalty scalar",  ifcClamped );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 0 = square plate with pressure; 1 = Scordelis Lo Roof; 2 = quarter hemisphere; 3 = pinched cylinder",  testCase );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("composite", "Composite material", composite);
    cmd.addSwitch( "nl", "Print information", nonlinear );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Geometry Setup]
    gsMultiPatch<> mp;
    std::string fn;

    gsReadFile<>("surfaces/cylinder_4p.xml", mp);

    gsDebug<<mp;

    //! [Apply Boundary Condition]
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    // Initialize surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    real_t Load = -1e4;
    gsVector<> neu(3);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    neu << 0, 0, Load;
    gsConstantFunction<> neuData(neu,3);

    //Buckling coefficient
    real_t pressure = 0.0;

    // Apply pressure to boundary


    BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
    BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
    BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z



    BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false,0 ); // unknown 1 - y
    BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false,1 ); // unknown 2 - z
    BCs.addCondition(0,boundary::north, condition_type::collapsed, 0, 0, false,2 ); // unknown 0 - x
    BCs.addCondition(0,boundary::north, condition_type::neumann, &neuData); // unknown 0 - x

    BCs.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
    BCs.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
    BCs.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z


    BCs.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false,0 ); // unknown 1 - y
    BCs.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false,1 ); // unknown 2 - z
    BCs.addCondition(1,boundary::north, condition_type::collapsed, 0, 0, false,2 ); // unknown 0 - x
    BCs.addCondition(1,boundary::north, condition_type::neumann, &neuData); // unknown 0 - x

    BCs.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
    BCs.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
    BCs.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z



    BCs.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false,0 ); // unknown 1 - y
    BCs.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false,1 ); // unknown 2 - z
    BCs.addCondition(2,boundary::north, condition_type::collapsed, 0, 0, false,2 ); // unknown 0 - x
    BCs.addCondition(2,boundary::north, condition_type::neumann, &neuData); // unknown 0 - x

    BCs.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
    BCs.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
    BCs.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z



    BCs.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false,0 ); // unknown 1 - y
    BCs.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false,1 ); // unknown 2 - z
    BCs.addCondition(3,boundary::north, condition_type::collapsed, 0, 0, false,2 ); // unknown 0 - x
    BCs.addCondition(3,boundary::north, condition_type::neumann, &neuData); // unknown 0 - x

    gsFunctionExpr<> surfForce(tx,ty,tz,3);

    // Initialize solution obejct
    gsMultiPatch<> mp_def = mp;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp, 3);
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    std::vector<gsFunction<>*> parameters; //SvK & Composites & NH & NH_ext
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;

    index_t Compressibility = 0;

    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);

    gsBucklingSolver<real_t> buckling(K_L,rhs,dK_NL);
    buckling.options().setInt("solver",2);
    buckling.options().setInt("selectionRule",0);
    buckling.options().setInt("sortRule",4);
    buckling.options().setSwitch("verbose",true);
    buckling.options().setInt("ncvFac",2);






//    delete materialMatrix;
//    delete assembler;
//    return result;

}








