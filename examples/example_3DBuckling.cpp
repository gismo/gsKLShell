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

template <class T>
gsMultiPatch<T>