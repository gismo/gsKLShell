/** @file example_pointload.cpp

    @brief Asymmetric loading applied with an offset angle (pi/50)
    from the crown of a semicircular arch.

    Example 13 from Liu et al 2006

    Liu, B., Xing, Y., Wang, Z., Lu, X., & Sun, H. (2017) Non-uniform rational Lagrange functions and its applications to isogeometric analysis of in-plane and flexural vibration of thin plates. Computer Methods in Applied Mechanics and Engineering, 321, 173-208. http://dx.doi.org/10.1016/j.cma.2017.04.007

    Author(s): J. Li
 **/

//![Include namespace]
#include <gismo.h>
#include <gsKLShell/gsThinShellUtils.h>
#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

using namespace gismo;

int main(int argc, char *argv[]){

    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool Compressibility = false;
    index_t material = 0;
    bool verbose = false;
    std::string fn;

    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

//    real_t Ratio = 7.0;

    gsCmdLine cmd("2D clamped open sectiorial membrane.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
//    cmd.addReal( "R", "Ratio", "Mooney Rivlin Ratio",  Ratio );
//    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("comp", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("composite", "Composite material", composite);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Material properties]
    real_t E_modulus = 7.1e10; //Pa
    real_t PoissonRatio = 0.3;
    real_t Density = 2700.0; //kg/m^3
    real_t thickness = 0.005; //mm
    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    real_t omega = 6.8138; // omega2 = 8.2666 omega3 = 12.855 omega4 = 13.742
    gsInfo<<"Eigen frequency: "<< omega;
    //! [Material properties]

    //! [Make geometry, refine and evaluate]
    gsMultiPatch<> mp, mp_vib;
    real_t r0 = 0.5; //m
    real_t r1 = 1.0;
    gsTensorNurbs<2> geo = *gsNurbsCreator<>::NurbsQuarterAnnulus(r0,r1);
    mp.addPatch(geo);
    mp.addAutoBoundaries();
    real_t PI = 3.1415926535;
    real_t Area = 0.25*PI*(pow(r1,2) - pow(r0,2));


    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();
    // Set the deformed configuration
    mp_vib = mp;
    if (plot) gsWriteParaview<>( mp_vib  , "mp", 1000, true);
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";
    //! [Make geometry, refine and evaluate]

    //! [Set boundary conditions]
    real_t U_inf = 1;
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    gsFunctionExpr<> phi_in("0",2);
    gsFunctionExpr<> phi_out(std::to_string(U_inf*r0*r1/(r1*r1)) + "*x",2);
    gsFunctionExpr<> symmetry("0",2);
    gsFunctionExpr<> wall(std::to_string(r0) + "*1/(sqrt(1+y^2/x^2))",2);

    BCs.addCondition(0, boundary::north, condition_type::dirichlet, phi_in);
    BCs.addCondition(0, boundary::east, condition_type::dirichlet, phi_out);
    BCs.addCondition(0, boundary::south, condition_type::neumann, symmetry);
    BCs.addCondition(0, boundary::west, condition_type::neumann, wall);



    //! [Set boundary conditions]



}