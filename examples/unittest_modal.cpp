/** @file unittest_modal.cpp

    @brief Simple unit test for shells using modal analysis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

using namespace gismo;

gsVector<real_t> numerical(bool composite)
{
    // Input options
    int numElevate  = 2;
    int numHref     = 4;

    real_t thickness = 0.01;
    real_t E_modulus = 1e5;
    real_t Density = 1e0;
    real_t PoissonRatio = 0.3;

    gsMultiPatch<> mp;

    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);

    for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();

    gsMultiBasis<> dbasis(mp);

    // Boundary conditions
    gsBoundaryConditions<> BCs;

    // Plate
    // Pinned-Pinned-Pinned-Pinned
        // Left
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        // Right
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        // Top
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        // Bottom
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    // Initialise solution object
    gsMultiPatch<> mp_def = mp;

    // Linear isotropic material model
    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsConstantFunction<> rho(Density,3);

    // Linear anisotropic material model
    real_t pi = math::atan(1)*4;
    index_t kmax = 5;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Rs(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    Rs[0] = Rs[1] = Rs[2] = Rs[3] = Rs[4] = &rho;


    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = Gs[1] = Gs[2] = Gs[3] = Gs[4] = &Gfun;

    gsConstantFunction<> phi1, phi2, phi3, phi4, phi5;
    phi1.setValue(0/kmax * pi / 2.0,3);
    phi2.setValue(1/kmax * pi / 2.0,3);
    phi3.setValue(2/kmax * pi / 2.0,3);
    phi4.setValue(3/kmax * pi / 2.0,3);
    phi5.setValue(4/kmax * pi / 2.0,3);

    Phis[0] = &phi1;
    Phis[1] = &phi2;
    Phis[2] = &phi3;
    Phis[3] = &phi4;
    Phis[4] = &phi5;

    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = Ts[1] = Ts[2] = Ts[3] = Ts[4] = &thicks;

    std::vector<gsFunction<>*> parameters;
    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;

    if (composite)
    {
        materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis,Rs);
    }
    else
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    // options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    // options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",static_cast<int>(!composite));
    // materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

    assembler->assemble();
    gsSparseMatrix<> K =  assembler->matrix();
    assembler->assembleMass();
    gsSparseMatrix<> M =  assembler->matrix();

    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;
    eigSolver.compute(K,M);
    gsMatrix<> values  = eigSolver.eigenvalues();
    gsMatrix<> vectors = eigSolver.eigenvectors();

    values = values.cwiseSqrt();
    values = values.col(0).head(10);

    delete assembler;
    return values;
}

gsVector<real_t> analytical()
{
    real_t thickness = 0.01;
    real_t E_modulus = 1e5;
    real_t Density = 1e0;
    real_t PoissonRatio = 0.3;

    real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));

    std::vector<real_t> omegas;
    for (index_t m=1; m!=10; m++)
      for (index_t n=1; n!=10; n++)
        omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

    std::sort(omegas.begin(),omegas.end());
    omegas.resize(10);
    gsAsVector<> analytical(omegas);

    return analytical;
}

int main (int argc, char** argv)
{
    gsVector<> an = analytical();
    gsVector<> num, relError;
    real_t error;

    std::vector<bool> composite { true, false };

    for (std::vector<bool>::iterator comp = composite.begin(); comp!=composite.end(); comp++)
    {
        num = numerical(*comp);
        relError = (num - an).array()/an.array();
        error = relError.norm();

        if (error < 1e-3)
            gsInfo<<"Passed\n";
        else
            gsInfo<<"Failed\n";
    }
    return 1;
}
