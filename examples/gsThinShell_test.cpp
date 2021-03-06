/** @file gsThinShell_test2.cpp

    @brief Example testing and debugging thin shell solver. Based on gsThinShell_test

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixComposite.h>
#include <gsKLShell/gsMaterialMatrixBase.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool Compressibility = 0;
    index_t material = 0;
    bool verbose = false;
    std::string fn;
    bool membrane = false;


    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;
    real_t Ratio = 7.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addReal( "R", "Ratio", "Mooney Rivlin Ratio",  Ratio );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addSwitch( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    if (testCase == 1 )
    {
        thickness = 0.25;
        E_modulus = 4.32E8;
        fn = "surface/scordelis_lo_roof.xml";
        gsReadFile<>(fn, mp);
        PoissonRatio = 0.0;
    }
    else if (testCase == 2)
    {
        thickness = 0.04;
        E_modulus = 6.825E7;
        PoissonRatio = 0.3;
        gsReadFile<>("surface/quarter_hemisphere.xml", mp);
    }
    else if (testCase == 3)
    {
        thickness = 3;
        E_modulus = 3E6;
        PoissonRatio = 0.3;
        gsReadFile<>("surface/pinched_cylinder.xml", mp);
    }
    else if (testCase == 4)
    {
        real_t L = 1.0;
        real_t B = 1.0;
        real_t mu = 1.5e6;
        thickness = 0.001;
        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;
        E_modulus = 2*mu*(1+PoissonRatio);
        // PoissonRatio = 0;
        mp = Rectangle(L, B);

        gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";
    }
    else if (testCase == 5)
    {
        thickness = 1;
        E_modulus = 1;
        PoissonRatio = 0.5;
        gsReadFile<>("planar/unitcircle.xml", mp);
    }
    else if (testCase == 6)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        PoissonRatio = 0.3;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("surface/quarter_sphere.xml", mp);
    }
    else if (testCase == 7)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        PoissonRatio = 0.3;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("surface/quarter_frustrum.xml", mp);
    }
    // else if (testCase == 8)
    // {
    //     thickness = 1;
    //     PoissonRatio = 0.4998999999;
    //     E_modulus = 240.5653612;
    //     gsReadFile<>("surfaces/cooks_membrane.xml", mp);
    //     mp.embed(3);

    //     real_t mu = E_modulus/(2*(1+PoissonRatio));
    //     real_t K = E_modulus/(3-6*PoissonRatio);
    //     gsDebug<<"K = "<<K<<"\tmu = "<<mu<<"\n";
    // }
    else if (testCase == 8 || testCase == 9)
    {
        thickness = 0.1;
        real_t mu = 4.225;
        PoissonRatio = 0.5;
        E_modulus = 2*mu*(1+PoissonRatio);
        // gsReadFile<>("quarter_frustrum.xml", mp);

        // R1 is radius on bottom, R2 is radius on top
        mp = FrustrumDomain(numRefine,numElevate+2,2.0,1.0,1.0);
    }
    else if (testCase == 10)
    {
        E_modulus = 1;
        // thickness = 0.15;
        thickness = 1;
        if (!Compressibility)
          PoissonRatio = 0.499;
        else
          PoissonRatio = 0.45;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);

    }
    else if (testCase == 17)
    {
        // Unit square
        mp = RectangularDomain(0,2,1.0,1.0);
        E_modulus = 4.497000000e6;
        thickness = 0.001;
        PoissonRatio = 0.4999;
    }
    else if (testCase == 18)
    {
        // Unit square
        // gsReadFile<>("planar/annulus_4p.xml", mp);
        gsReadFile<>("surface/frustrum.xml", mp);
        mp.computeTopology();

        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0;
    }
    else if (testCase == 20)
    {
        // Unit square
        // gsReadFile<>("planar/annulus_4p.xml", mp);
        gsMultiPatch<> mp_old;
        mp_old.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp = mp_old.uniformSplit();
        mp.computeTopology();
        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0;
    }
    else
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.embed(3);
        mp.addAutoBoundaries();
        E_modulus = 1e0;
        thickness = 1e0;
        // PoissonRatio = 0.5;
        // PoissonRatio = 0.499;
        PoissonRatio = 0.0;
        // PoissonRatio = 0.5;
    }
    //! [Read input file]
    // p-refine
    if (testCase != 8 && testCase != 9)
    {
        if (numElevate!=0)
            mp.degreeElevate(numElevate);

        // h-refine
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine();
    }
    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);


    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    gsConstantFunction<> displx(0.1,3);
    gsConstantFunction<> disply(0.25,3);

    gsFunctionExpr<> neuDataFun1;
    gsConstantFunction<> neuData(neu,3);
    real_t pressure = 0.0;
    if (testCase == 0)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,i);
        }

        // bc.addCondition(boundary::north, condition_type::clamped, 0, 0 ,false,2);
        // bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
        // bc.addCondition(boundary::south, condition_type::clamped, 0, 0 ,false,2);
        // bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        // tmp << 0,0,0;
        tmp << 0,0,-1;

        // Point loads
        // gsVector<> point(2);
        // gsVector<> load (3);
        // point<< 0.5, 0.5 ; load << 0.0, 1.0, 0.0 ;
        // pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 1)
    {
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z

        bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z

        // Surface forces
        tmp << 0, 0, -90;
    }
    else if (testCase == 2)
    {
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Surface forces
        tmp.setZero();

        // Point loads
        gsVector<> point(2);
        gsVector<> load (3);
        point<< 0.0, 0.0 ; load << 1.0, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
        point<< 1.0, 0.0 ; load << 0.0, -1.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 3)
    {
        // Symmetry in y-direction for back side
        bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

        // Diaphragm conditions for left side
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction: for right side
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in z-direction:for the front side
        bc.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 );

        // Surface forces
        tmp.setZero();

        // Point loads
        gsVector<> point(2); point<< 1.0, 1.0 ;
        gsVector<> load (3); load << 0.0, 0.0, -0.25 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 4) // Uniaxial tension; use with hyperelastic material model!
    {
        neu << 2625, 0, 0;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::neumann, &neuData ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.


        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

        gsVector<> point(2);
        gsVector<> load (3);
        point<< 1.0, 0.5 ;
        load << 1.0,0.0, 0.0;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 5)
    {
        for (index_t i = 0; i!=3; i++)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
        }

        // Point loads
        gsVector<> point(2); point<< 0.5,0.5 ;
        gsVector<> load (3); load << 0.0, 0.0, -1 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 6) // balloon
    {
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        // bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Pressure
        pressure = 5e3;
    }
    else if (testCase == 7)
    {
        neu << 0, 0, -100.0;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::north, condition_type::neumann, &neuData );
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );
    }
    // else if (testCase == 8)
    // {
    //     bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    //     bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

    //     // Symmetry in x-direction:
    //     neu << 0, 6.25, 0;
    //     neuData.setValue(neu,3);
    //     bc.addCondition(boundary::east, condition_type::neumann, &neuData );

    //     bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
    //     bc.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, 2 );
    //     bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
    //     bc.addCondition(boundary::west,  condition_type::dirichlet, 0, 0, false, 2 );
    // }
    else if (testCase == 8)
    {
    //     real_t Load = -0.0168288;
    //     neu << 0, 0, Load;
    //     neuData.setValue(neu,3);

    //     bc.addCondition(boundary::north, condition_type::neumann, &neuData );
    //     bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
    //     bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
    //     bc.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

        displx.setValue(-0.027815,3);
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
        bc.addCondition(boundary::north, condition_type::dirichlet, &displx, 0, false, 2 ); // unknown 1 - y

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

    }
    else if (testCase == 9)
    {
        real_t  Load = -0.0168288;
        neu << 0, 0, Load;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::north, condition_type::neumann, &neuData );
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
        bc.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

    }
    else if (testCase == 10)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        // bc.addCondition(boundary::east, condition_type::dirichlet, &displx, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (3); load << 0.25, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 11)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }

        bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        tmp<<0,0,-1;
    }
    else if (testCase == 12)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }

        bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        neu << 0, 0, -0.1;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
    }
    else if (testCase == 13)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }

        bc.addCondition(boundary::north, condition_type::clamped, 0, 0 ,false,2);
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
        bc.addCondition(boundary::south, condition_type::clamped, 0, 0 ,false,2);
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        tmp << 0,0,-1;
    }
    else if (testCase == 14)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        tmp << 0,0,-1;

        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 2 );
    }
    else if (testCase == 14)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        tmp << 0,0,-1;

        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 2 );
    }
    else if (testCase == 15)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z


        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (3); load << 0.1, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 16)
    {
      bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        real_t Load = 1e-2;
        gsVector<> point(2);
        gsVector<> load (3);
        point<< 1.0, 0.5 ;
        load << 0.0, 0.0, Load ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 17)
    {
        real_t Load = 1e-1;
        neu << -Load, 0, 0;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::collapsed, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    }
    else if (testCase == 20)
    {
        real_t Load = 1e-1;
        neu << 0, 0, -Load;
        neuData.setValue(neu,3);

        real_t T1 = 0.5;
        std::string nx = std::to_string(-T1) + "*cos(atan2(y,x))";
        std::string ny = std::to_string(-T1) + "*sin(atan2(y,x))";
        std::string nz = "0";
        neuDataFun1 = gsFunctionExpr<>(nx,ny,nz,3);

        // for (index_t p=0; p!=mp.nPatches(); p++)
        // {
        //     // bc.addCondition(p,boundary::west, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        //     bc.addCondition(p,boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        //     bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 1 - y
        //     bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        //     bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        // }

        bc.addCondition(0,boundary::south, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        // bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        // bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        // bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x

        // bc.addCondition(1,boundary::west, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x

        // bc.addCondition(2,boundary::south, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x

        // bc.addCondition(3,boundary::east, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x


        tmp << 0,0,-1e-6;
    }

    gsDebugVar(bc);
    //! [Refinement]

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsConstantFunction<> pressFun(pressure,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(Ratio,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));

    gsConstantFunction<> alpha1(2.0,3);
    gsConstantFunction<> mu1(7.0*mu/8.0,3);
    gsConstantFunction<> alpha2(-2.0,3);
    gsConstantFunction<> mu2(-mu/8.0,3);
    // gsConstantFunction<> alpha3()
    // gsConstantFunction<> mu3()

    // gsMaterialMatrix materialMatrixNonlinear(mp,mp_def,t,E,nu,rho);
    std::vector<gsFunction<>*> parameters(3);
    parameters[0] = &E;
    parameters[1] = &nu;
    parameters[2] = &ratio;
    gsMaterialMatrixBase<real_t>* materialMatrix;

    // Linear anisotropic material model
    real_t pi = math::atan(1)*4;
    index_t kmax = 1;
    gsVector<> E11(kmax), E22(kmax), G12(kmax), nu12(kmax), nu21(kmax), thick(kmax), phi(kmax);
    E11.setZero(); E22.setZero(); G12.setZero(); nu12.setZero(); nu21.setZero(); thick.setZero(); phi.setZero();
    for (index_t k=0; k != kmax; ++k)
    {
        E11.at(k) = E22.at(k) = E_modulus;
        nu12.at(k) = nu21.at(k) = PoissonRatio;
        G12.at(k) = 0.5 * E_modulus / (1+PoissonRatio);
        thick.at(k) = thickness/kmax;
        phi.at(kmax) = k / kmax * pi/2.0;
    }

    gsConstantFunction<> E11fun(E11,3);
    gsConstantFunction<> E22fun(E22,3);
    gsConstantFunction<> G12fun(G12,3);
    gsConstantFunction<> nu12fun(nu12,3);
    gsConstantFunction<> nu21fun(nu21,3);
    gsConstantFunction<> thickfun(thick,3);
    gsConstantFunction<> phifun(phi,3);

    if (material==0)
    {
        materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,mp_def,t,parameters,rho);
    }
    if (material==1)
    {
        // materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,mp_def,t,parameters,rho);
        materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,mp_def,E11fun,E22fun,G12fun,nu12fun,nu21fun,thickfun,phifun);
    }
    else if ((material==2) && (!Compressibility))
    {
        constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
        materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    }
    else if ((material==2) && (Compressibility))
    {
        constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
        materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    }
    else if ((material==12) && (!Compressibility))
    {
        constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
        materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    }
    else if ((material==12) && (Compressibility))
    {
        constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
        materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    }
    else if ((material==22) && (!Compressibility))
    {
        constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
        materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    }
    else if ((material==22) && (Compressibility))
    {
        constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
        materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    }

    // else if ((material==3) && (!Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==3) && (Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==13) && (!Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==13) && (Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==23) && (!Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==23) && (Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    // }

    // else if ((material==14) && (!Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==14) && (Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    // }

    // else if ((material==5) && (!Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Analytical>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==5) && (Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Analytical>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==15) && (!Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Spectral>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==15) && (Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Spectral>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==25) && (!Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Generalized>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,false>(mp,mp_def,t,parameters,rho);
    // }
    // else if ((material==25) && (Compressibility))
    // {
    //     constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Generalized>::id;
    //     materialMatrix = new gsMaterialMatrix<3,real_t,id,true>(mp,mp_def,t,parameters,rho);
    // }


    std::vector<gsFunction<>*> parameters2(6);
    if (material==14)
    {
        parameters2[0] = &E;
        parameters2[1] = &nu;
        parameters2[2] = &mu1;
        parameters2[3] = &alpha1;

        parameters2[4] = &mu2;
        parameters2[5] = &alpha2;

        // parameters[6] = ;
        // parameters[7] = ;
        materialMatrix->setParameters(parameters2);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    if(membrane)
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,bc,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

    assembler->setPointLoads(pLoads);
    if (pressure!= 0.0)
        assembler->setPressure(pressFun);

    // gsVector<> found_vec(3);
    // found_vec<<0,0,2;
    // gsConstantFunction<> found(found_vec,3);
    // assembler->setFoundation(found);

    gsStopwatch stopwatch,stopwatch2;
    real_t time = 0.0;
    real_t totaltime = 0.0;

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      time += stopwatch.stop();
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      time += stopwatch.stop();
      return assembler->rhs();
    };

    // Define Matrices
    stopwatch.restart();
    stopwatch2.restart();
    assembler->assemble();
    time += stopwatch.stop();

    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();

    gsDebugVar(assembler->matrix().toDense());
    gsDebugVar(assembler->rhs().transpose());

    // Solve linear problem
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    gsDebugVar(solVector.transpose());

    real_t residual = vector.norm();
    real_t residual0 = residual;
    real_t residualOld = residual;
    gsVector<real_t> updateVector = solVector;
    gsVector<real_t> resVec = Residual(solVector);
    gsSparseMatrix<real_t> jacMat;
    for (index_t it = 0; it != 100; ++it)
    {
        jacMat = Jacobian(solVector);
        solver.compute(jacMat);
        updateVector = solver.solve(resVec); // this is the UPDATE
        solVector += updateVector;

        gsDebugVar(assembler->matrix().toDense());
        gsDebugVar(assembler->rhs().transpose());
        gsDebugVar(updateVector.transpose());

        resVec = Residual(solVector);
        residual = resVec.norm();

        gsInfo<<"Iteration: "<< it
           <<", residue: "<< residual
           <<", update norm: "<<updateVector.norm()
           <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
           <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
           <<"\n";

        residualOld = residual;

        if (updateVector.norm() < 1e-6)
            break;
        else if (it+1 == it)
            gsWarn<<"Maximum iterations reached!\n";


            // ADD DIRICHLET HOMOGENIZE
    }

    totaltime += stopwatch2.stop();

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    if (testCase==4)
    {
        gsVector<> pt(2);
        pt<<1,0;
      gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);
      real_t S = 2625 / 1e-3 / lambdas(0) / lambdas(2);
      real_t San = mu * (math::pow(lambdas(1),2)-1/lambdas(1));
      gsInfo<<"S = \t"<<S<<"\t San = \t"<<San<<"\t |S-San| = \t"<<abs(S-San)<<"\n";
      gsInfo<<"lambda = \t"<<lambdas(1)<<"\t 1/lambda = \t"<<1/lambdas(1)<<"\t lambda_0 = \t"<<lambdas(0)<<"\t lambda_2 = \t"<<lambdas(2)<<"\n";
    }


    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( u_sol   , G, "solution");

        // gsFileManager::open("solution.pvd");

        gsInfo <<"Maximum deformation coef: "
               << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
        gsInfo <<"Minimum deformation coef: "
               << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    }
    if (stress)
    {

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsPiecewiseFunction<> stretches;
        assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
        gsField<> Stretches(mp_def,stretches, true);

        // gsPiecewiseFunction<> membraneStresses_p;
        // assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
        // gsField<> membraneStress_p(mp_def,membraneStresses_p, true);

        // gsPiecewiseFunction<> flexuralStresses_p;
        // assembler->constructStress(mp_def,flexuralStresses_p,stress_type::principal_stress_flexural);
        // gsField<> flexuralStress_p(mp_def,flexuralStresses_p, true);

        gsPiecewiseFunction<> stretch1;
        assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
        gsField<> stretchDir1(mp_def,stretch1, true);

        gsPiecewiseFunction<> stretch2;
        assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
        gsField<> stretchDir2(mp_def,stretch2, true);

        gsPiecewiseFunction<> stretch3;
        assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
        gsField<> stretchDir3(mp_def,stretch3, true);


        gsField<> solutionField(mp,deformation, true);


        // gsField<> stressField = assembler->constructStress(mp_def,stress_type::membrane_strain);

        gsWriteParaview(solutionField,"Deformation");
        gsWriteParaview(membraneStress,"Membrane Stress");
        gsWriteParaview(flexuralStress,"Flexural Stress");
        gsWriteParaview(stretchDir1,"Principal Direction 1");
        gsWriteParaview(stretchDir2,"Principal Direction 2");
        gsWriteParaview(stretchDir3,"Principal Direction 3");
    }
    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    return EXIT_SUCCESS;

}// end main

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomain(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> Rectangle(T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 2; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,2,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,2,1);

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}


template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomainVert(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains z-coordinates (height)
    coefs.col(2).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomain90(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(1).segment(k*len0,len0) = coefvec0;
    // Second column contains z-coordinates (height)
    coefs.col(0).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  // n = number of uniform refinements over the height; n = 0, only top and bottom part

  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,3,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,3,1);

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      kv1.uniformRefine();

  gsDebug<<kv1;

  // Make basis
  // gsTensorNurbsBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  index_t N = math::pow(2,n)+2;
  gsMatrix<> coefs(3*N,dim);
  gsMatrix<> tmp(3,3);
  T R,H;

  gsMatrix<> weights(3*N,1);
  for (index_t k=0; k!= N; k++)
  {
    R = k*(R2-R1)/(N-1) + R1;
    H = k*h/(N-1);
    tmp<< R,0,H,
          R,R,H,
          0,R,H;

    coefs.block(3*k,0,3,3) = tmp;

    weights.block(3*k,0,3,1) << 1,0.70711,1;
  }

  // Create gsGeometry-derived object for the patch
  gsTensorNurbs<2,real_t> shape(kv0,kv1,coefs,weights);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  // Elevate up to order p
  if (p>2)
  {
    for(index_t i = 2; i< p; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree
  }

  // // Refine n times
  // for(index_t i = 0; i< n; ++i)
  //     mp.patch(0).uniformRefine();

  return mp;
}


/*
    to do:
    =  make function for construction of the solution given the space and the mp
*/



/*
template<class T>
void gsShellAssembler<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<unsigned> acts,globalActs;

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        if ( m_pLoads[i].parametric )
        {
            m_bases.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_bases.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else
        {
            gsMatrix<> forcePoint;
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, acts );
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, bVals);
        }

        // translate patch-local indices to global dof indices
        for (size_t j = 0; j< 3; ++j)
        {
            if (m_pLoads[i].value[j] != 0.0)
            {
                u.dofMappers[j].localToGlobal(acts, m_pLoads[i].patch, globalActs);

                for (index_t k=0; k < globalActs.rows(); ++k)
                {
                    if (int(globalActs(k,0)) < m_dofs)
                        m_rhs(globalActs(k,0), 0) += bVals(k,0) * m_pLoads[i].value[j];
                }
            }
        }
    }
}
*/
