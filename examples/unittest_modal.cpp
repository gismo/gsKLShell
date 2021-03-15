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
    gsFunctionExpr<> rho(std::to_string(Density),3);

    // Linear anisotropic material model
    real_t pi = math::atan(1)*4;
    index_t kmax = 5;
    gsVector<> E11(kmax), E22(kmax), G12(kmax), nu12(kmax), nu21(kmax), thick(kmax), phi(kmax);
    E11.setZero(); E22.setZero(); G12.setZero(); nu12.setZero(); nu21.setZero(); thick.setZero(); phi.setZero();
    for (index_t k=0; k != kmax; ++k)
    {
        E11.at(k) = E22.at(k) = E_modulus;
        nu12.at(k) = nu21.at(k) = PoissonRatio;
        G12.at(k) = 0.5 * E_modulus / (1+PoissonRatio);
        thick.at(k) = thickness/kmax;
        phi.at(k) = static_cast<real_t>(k) / kmax * pi/2.0;
    }

    gsConstantFunction<> E11fun(E11,3);
    gsConstantFunction<> E22fun(E22,3);
    gsConstantFunction<> G12fun(G12,3);
    gsConstantFunction<> nu12fun(nu12,3);
    gsConstantFunction<> nu21fun(nu21,3);
    gsConstantFunction<> thickfun(thick,3);
    gsConstantFunction<> phifun(phi,3);

    std::vector<gsFunction<>*> parameters;
    if (!composite)
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
    }
    else
    {
        parameters.resize(6);
        parameters[0] = &E11fun;
        parameters[1] = &E22fun;
        parameters[2] = &G12fun;
        parameters[3] = &nu12fun;
        parameters[4] = &nu21fun;
        parameters[5] = &phifun;
    }

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;

    if (composite)
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",0);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,mp_def,t,parameters,rho,options);
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,mp_def,t,parameters,rho,options);
    }

    // options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    // options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",static_cast<int>(!composite));
    // materialMatrix = getMaterialMatrix<3,real_t>(mp,mp_def,t,parameters,rho,options);

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
