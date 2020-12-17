/** @file gsMaterialMatrix_test.cpp

    @brief Code for the arc-length method of a shell based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsKLShell/gsMaterialMatrix.h>

using namespace gismo;

int main (int argc, char** argv)
{

	int material = 2;
	int Compressibility = 0;
	gsCmdLine cmd("Thin shell plate example.");

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsMultiPatch<> mp, mp_def;

    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);
    mp_def = mp;


    gsMatrix<> ones;
    // translate second patch up
    ones = gsMatrix<>::Ones(mp.patch(0).coefs().rows(),1); // patch 1
    mp_def.patch(0).coefs().col(1) += ones; // patch 1

    real_t thickness = 1.0;
    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.499;
    real_t Ratio = 2;
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsConstantFunction<> ratio(Ratio,3);


    std::vector<gsFunction<>*> parameters(3);
    parameters[0] = &E;
    parameters[1] = &nu;
    parameters[2] = &ratio;

    gsMaterialMatrix mm(mp,mp_def,t,parameters);
	mm.options().setInt("MaterialLaw",material);
    mm.options().setInt("Compressibility",Compressibility);

    gsVector<> pt(2);
    pt.setConstant(0.25);

    gsMatrix<> result;
    mm.makeMatrix(0);
    mm.eval_into(pt,result);
    gsInfo<<"matrix 0 = \n"<<result.reshape(3,3)<<"\n";

    mm.makeMatrix(1);
    mm.eval_into(pt,result);
    gsInfo<<"matrix 1 = \n"<<result.reshape(3,3)<<"\n";

    mm.makeMatrix(2);
    mm.eval_into(pt,result);
    gsInfo<<"matrix 2 = \n"<<result.reshape(3,3)<<"\n";

    mm.makeVector(0);
    mm.eval_into(pt,result);
    gsInfo<<"vector 0 = \n"<<result.reshape(3,1)<<"\n";

    mm.makeVector(1);
    mm.eval_into(pt,result);
    gsInfo<<"vector 1 = \n"<<result.reshape(3,1)<<"\n";

	return 0;
}