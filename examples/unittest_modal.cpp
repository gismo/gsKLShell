/** @file gsThinShell_Buckling.cpp

    @brief Example to compute eigenvalues and eigenmodes of a buckled shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

using namespace gismo;


int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 2;
    int numHref     = 4;

    real_t thickness     = 1;
    real_t E_modulus     = 1e0;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;

    gsMultiPatch<> mp;


    gsCmdLine cmd("Modal analysis for thin shells.");

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


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

    std::vector<real_t> omegas;

    thickness = 0.01;
    E_modulus = 1e5;
    Density = 1e0;
    PoissonRatio = 0.3;
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


    real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
    for (index_t m=1; m!=10; m++)
      for (index_t n=1; n!=10; n++)
        omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

    std::sort(omegas.begin(),omegas.end());
    omegas.resize(10);
    gsAsVector<> analytical(omegas);

    // Initialise solution object
    gsMultiPatch<> mp_def = mp;

    // Linear isotropic material model
    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,mp_def,t,parameters,rho,options);

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
    analytical = analytical.head(10);
    gsVector<> relError = (values - analytical).array()/values.array();

    real_t error = relError.norm();
    gsDebugVar(error);

    if (error < 1e-3)
        gsInfo<<"Passed";
    else
        gsInfo<<"Failed";


    return (error < 1e-3);
}
