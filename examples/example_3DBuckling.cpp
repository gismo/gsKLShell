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
    bool plot       = false;
    bool sparse     = false;
    bool nonlinear  = false;
    bool first  = false;
    int mode = 0;

    real_t E_modulus     = 1e8;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;

    real_t thickness = 1e-3;
    real_t aDim = 1.0;
    real_t bDim = 1.0;

    index_t Compressibility = 0;
    index_t material = 0;
    real_t Ratio = 7.0;
    bool composite = false;
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

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp;
    gsReadFile<>("surfaces/cylinder.xml");
}