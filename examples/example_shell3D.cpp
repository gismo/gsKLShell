/** @file example_shell3D.cpp

    @brief Simple 3D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixContainer.h>
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsCore/gsPiecewiseFunction.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> Plate(T Lp, T Wp, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> Strip(T Lb, T Hw, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> IBeam(T Lb, T Hw, T Wf, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> TBeam(T Lb, T Hw, T Wf, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> LBeam(T Lb, T Hw, T Wf, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> PanelT(T Lp, T Wp, T Hw, T Wf, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> PanelL(T Lp, T Wp, T Hw, T Wf, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> PlateGirderL(T PanelLength, T PanelWidth, T GirderHeight, T GirderFlangeWidth, T WebHeight, T WebFlangeWidth, T x = 0, T y = 0, T z = 0);


// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool nonlinear = false;

    bool membrane = false;
    bool composite = false;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;

    real_t ifcDirichlet = 1.0;
    real_t ifcClamped = 1.0;


    gsCmdLine cmd("2D shell example.");
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
    //! [Parse command line]

    //! [Define material parameters and geometry per example]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    std::string fn;
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
        thickness = 1;
        E_modulus = 1;
        PoissonRatio = 0.0;

        // MULTIPATCH CASE plate
        gsReadFile<>("planar/two_squares.xml", mp);
        mp.embed(3);
    }
    else if (testCase == 5)
    {
        thickness = 1;
        E_modulus = 1;
        PoissonRatio = 0.3;

        // MULTIPATCH CASE plate
        // mp = TBeam(10.,1.,0.5);

        gsReadFile<>("multipatches/T-beam.xml", mp);
    }
    else if (testCase == 6)
    {
        thickness = 1;
        E_modulus = 1;
        PoissonRatio = 0.3;

        // MULTIPATCH CASE plate
        gsReadFile<>("multipatches/I-beam.xml", mp);
    }
    else if (testCase == 7)
    {
        thickness = 1;
        E_modulus = 1;
        PoissonRatio = 0.3;

        mp = PanelT(10.,5.,1.,1.);
    }
    else if (testCase == 8)
    {
        thickness = 1;
        E_modulus = 1;
        PoissonRatio = 0.0;

        // MULTIPATCH CASE plate
        gsReadFile<>("planar/two_squares.xml", mp);
        mp.embed(3);
    }
    else if (testCase == 9)
    {
        thickness = 10.; // plate
        E_modulus = 206000;
        PoissonRatio = 0.3;
        mp = PanelL(3500.,3500.,200.,38.);
        gsWrite<>(mp,"PanelL");
    }
    else if (testCase == 10)
    {
        thickness = 10.; // plate
        E_modulus = 206000;
        PoissonRatio = 0.3;
        mp = PlateGirderL(7020.,3500.,580.,250.,200.,38.);
        gsWrite<>(mp,"PlateGirderL");

        // Plate
        // t = 10
        // Girder
        // tw = 15
        // Hw = 580
        // tf = 30
        // Wf = 250
        // Stiffener
        // tw = 10
        // Hw = 200
        // tf = 15
        // Wf = 38

    }
    else
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.embed(3);
        mp.addAutoBoundaries();
        E_modulus = 1e0;
        thickness = 1e0;
        PoissonRatio = 0.3;
    }
    //! [Define material parameters and geometry per example]

    //! [Refine and elevate]
    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);
    //! [Refine and elevate]

    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";


    //! [Set boundary conditions]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsPiecewiseFunction<> force(mp.nPatches());
    gsPiecewiseFunction<> t(mp.nPatches());
    // gsPiecewiseFunction<> nu(mp.nPatches());

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

    gsVector<> refPoint(2);
    real_t refPatch = 0;
    if (testCase == 0)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,i);
        }

        // pressure = -1;
        refPoint<<0.5,0.5;
        tmp << 0,0,-1;

        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);

        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);

        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);
    }
    else if (testCase == 1)
    {
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 1 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 2 );

        bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 1 );
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 2 );

        // Surface forces
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);

        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);

        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);

        refPoint<<0.5,1;
    }
    else if (testCase == 2)
    {
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 );

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Surface forces
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);

        // // material parameters
        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);

        // Point loads
        gsVector<> point(2);
        gsVector<> load (3);
        point<< 0.0, 0.0 ; load << 1.0, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
        point<< 1.0, 0.0 ; load << 0.0, -1.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );

        refPoint = point;
    }
    else if (testCase == 3)
    {
        // Symmetry in y-direction for back side
        bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

        // Diaphragm conditions for left side
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 );

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
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);

        // // material parameters
        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);

        // Point loads
        gsVector<> point(2); point<< 1.0, 1.0 ;
        gsVector<> load (3); load << 0.0, 0.0, -0.25 ;
        pLoads.addLoad(point, load, 0 );

        refPoint = point;
    }
    else if (testCase == 4)
    {
        for (index_t d = 0; d!=3; d++)
        {
            bc.addCondition(0, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::west, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(0, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::south, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::north, condition_type::dirichlet, 0, 0, false, d);
        }

        // Surface forces
        tmp << 0,0,-1e-3;
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> force1(tmp,3);
        force.addPiece(force1);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);
        gsFunctionExpr<> t1(std::to_string(thickness), 3);
        t.addPiece(t1);

        // // material parameters
        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);
        // gsFunctionExpr<> nu1(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu1);

        // Point loads
        gsVector<> point(2); point<< 0.0, 0.5 ;
        refPoint = point;
    }
    else if (testCase == 5)
    {
        for (index_t d = 0; d!=3; d++)
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::dirichlet, 0, 0, false, d);
        for (size_t p=0; p!=mp.nPatches(); ++p)
            bc.addCondition(p, boundary::east, condition_type::clamped, 0, 0, false, 2);


        // Surface forces
        tmp << 0,0,0;
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> force1(tmp,3);
        force.addPiece(force1);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> force2(tmp,3);
        force.addPiece(force2);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);
        gsFunctionExpr<> t1(std::to_string(thickness), 3);
        t.addPiece(t1);
        gsFunctionExpr<> t2(std::to_string(thickness), 3);
        t.addPiece(t2);

        // // material parameters
        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);
        // gsFunctionExpr<> nu1(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu1);
        // gsFunctionExpr<> nu2(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu2);

        // Point loads
        gsVector<> point(2); point<< 1.0, 1.0 ;
        refPoint = point;
    }
    else if (testCase == 6)
    {
        for (index_t d = 0; d!=3; d++)
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::dirichlet, 0, 0, false, d);
        for (size_t p=0; p!=mp.nPatches(); ++p)
            bc.addCondition(p, boundary::east, condition_type::clamped, 0, 0, false, 2);

        // Point loads
        tmp << 0,0,0;
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> force1(tmp,3);
        force.addPiece(force1);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> force2(tmp,3);
        force.addPiece(force2);
        tmp << 0,0,0;
        gsConstantFunction<> force3(tmp,3);
        force.addPiece(force3);
        tmp << 0,0,0;
        gsConstantFunction<> force4(tmp,3);
        force.addPiece(force4);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);
        gsFunctionExpr<> t1(std::to_string(thickness), 3);
        t.addPiece(t1);
        gsFunctionExpr<> t2(std::to_string(thickness), 3);
        t.addPiece(t2);
        gsFunctionExpr<> t3(std::to_string(thickness), 3);
        t.addPiece(t3);
        gsFunctionExpr<> t4(std::to_string(thickness), 3);
        t.addPiece(t4);

        // // material parameters
        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);
        // gsFunctionExpr<> nu1(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu1);
        // gsFunctionExpr<> nu2(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu2);
        // gsFunctionExpr<> nu3(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu3);
        // gsFunctionExpr<> nu4(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu4);

        gsVector<> point(2); point<< 1.0, 1.0 ;
        refPoint = point;
    }
    else if (testCase == 7)
    {
        for (index_t d = 0; d!=3; d++)
        {
            bc.addCondition(0, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(0, boundary::east, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(1, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::south, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::west, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(2, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(2, boundary::west, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(3, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(3, boundary::south, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(4, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(4, boundary::south, condition_type::dirichlet, 0, 0, false, d);

        }

        // Point loads
        tmp << 0,0,1e-2;
        gsConstantFunction<> piece0(tmp,3);
        force.addPiece(piece0);
        tmp << 0,0,1e-2;
        gsConstantFunction<> piece1(tmp,3);
        force.addPiece(piece1);
        tmp << 0,0,0;
        gsConstantFunction<> piece2(tmp,3);
        force.addPiece(piece2);
        tmp << 0,0,0;
        gsConstantFunction<> piece3(tmp,3);
        force.addPiece(piece3);
        tmp << 0,0,0;
        gsConstantFunction<> piece4(tmp,3);
        force.addPiece(piece4);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);
        gsFunctionExpr<> t1(std::to_string(thickness), 3);
        t.addPiece(t1);
        gsFunctionExpr<> t2(std::to_string(thickness), 3);
        t.addPiece(t2);
        gsFunctionExpr<> t3(std::to_string(thickness), 3);
        t.addPiece(t3);
        gsFunctionExpr<> t4(std::to_string(thickness), 3);
        t.addPiece(t4);

        // // material parameters
        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);
        // gsFunctionExpr<> nu1(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu1);
        // gsFunctionExpr<> nu2(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu2);
        // gsFunctionExpr<> nu3(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu3);
        // gsFunctionExpr<> nu4(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu4);

        gsVector<> point(2); point<< 1.0, 1.0 ;
        refPoint = point;
    }
    else if (testCase == 8)
    {
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0 ,false,0);

        bc.addCondition(0,boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 ,false,2);

        bc.addCondition(1,boundary::west, condition_type::clamped  , 0, 0, false,2);
        bc.addCondition(0,boundary::east, condition_type::clamped  , 0, 0, false,2);

        bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );

        // Surface forces
        tmp << 0,0,0;
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);
        tmp << 0,0,0;
        gsConstantFunction<> force1(tmp,3);
        force.addPiece(force1);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);
        gsFunctionExpr<> t1(std::to_string(thickness), 3);
        t.addPiece(t1);

        // // material parameters
        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);
        // gsFunctionExpr<> nu1(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu1);

        real_t Load = 1e0;
        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (3); load << Load,0.0, 0.0;
        pLoads.addLoad(point, load,0);
    }
    else if (testCase == 9)
    {
        // Horizontal displacements
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 1 );      // stiffener web
        bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 1 );      // stiffener web
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 1 );     // stiffener flange
        bc.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false, 1 );     // stiffener flange

        // Vertical displacements
        bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );      // plate
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );      // plate
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );      // stiffener web
        bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );      // stiffener web
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );     // stiffener flange
        bc.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );     // stiffener flange


        // Clampings
        bc.addCondition(0,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::east, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(2,boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(3,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(3,boundary::south, condition_type::clamped, 0, 0, false, 2 );

        std::vector<index_t> plateIDs{0,1};
        std::vector<index_t> stiffenerWebIDs{2};
        std::vector<index_t> stiffenerFlangeIDs{3};

        real_t pressure = 1;
        // Surface loads
        gsVector<> pressvec(3);
        pressvec << 0,0,pressure;
        tmp << 0,0,0;
        gsConstantFunction<> forcep(pressvec,3);
        gsConstantFunction<> force0(tmp,3);

        typedef typename std::vector<gsFunction<real_t>*>       FunctionContainer;
        FunctionContainer forceContainer(mp.nPatches());
        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            forceContainer.at(*it) = &forcep;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            forceContainer.at(*it) = &force0;

        for (size_t p = 0; p != forceContainer.size(); p++)
            force.addPiece(*forceContainer.at(p));


        FunctionContainer tContainer(mp.nPatches());

        gsVector<> tmp2(1);
        gsFunctionExpr<> tp(std::to_string(15.), 3);
        gsFunctionExpr<> tsw(std::to_string(5.), 3);
        gsFunctionExpr<> tsf(std::to_string(10.), 3);

        gsDebugVar(tsf.targetDim());

        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            tContainer.at(*it) = &tp;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            tContainer.at(*it) = &tsw;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            tContainer.at(*it) = &tsf;

        for (size_t p = 0; p != tContainer.size(); p++)
            t.addPiece(*tContainer.at(p));
    }
    else if (testCase == 10)
    {
        // Plate
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(3,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(4,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(5,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(5,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(0,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(3,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(3,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(4,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(5,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(5,boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Girder
        // Web
        bc.addCondition(13,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(20,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(16,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(21,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(13,boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(20,boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(16,boundary::east, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(21,boundary::east, condition_type::clamped, 0, 0, false, 0 );

        //Flange
        bc.addCondition(14,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(15,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(17,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(18,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(14,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(15,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(17,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(18,boundary::south, condition_type::clamped, 0, 0, false, 2 );

        // Stiffener
        // Web
        bc.addCondition(8,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(6,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(8,boundary::west, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(6,boundary::east, condition_type::clamped, 0, 0, false, 1 );
        //Flange
        bc.addCondition(9,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(7,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(9,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(7,boundary::north, condition_type::clamped, 0, 0, false, 2 );

        std::vector<index_t> plateIDs{0,1,2,3,4,5};
        std::vector<index_t> girderWebIDs{10,13,16,19,20,21};
        std::vector<index_t> girderFlangeIDs{11,12,14,15,17,18};
        std::vector<index_t> stiffenerWebIDs{6,8};
        std::vector<index_t> stiffenerFlangeIDs{7,9};

        real_t pressure = 1e0;
        // Surface loads
        gsVector<> pressvec(3);
        pressvec << 0,0,pressure;
        tmp << 0,0,0;
        gsConstantFunction<> forcep(pressvec,3);
        gsConstantFunction<> force0(tmp,3);

        typedef typename std::vector<gsFunction<real_t>*>       FunctionContainer;
        FunctionContainer forceContainer(mp.nPatches());
        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            forceContainer.at(*it) = &forcep;
        // Girder webs
        for (typename std::vector<index_t>::iterator it = girderWebIDs.begin(); it!=girderWebIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Girder flanges
        for (typename std::vector<index_t>::iterator it = girderFlangeIDs.begin(); it!=girderFlangeIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            forceContainer.at(*it) = &force0;

        for (size_t p = 0; p != forceContainer.size(); p++)
            force.addPiece(*forceContainer.at(p));

        // force = gsPiecewiseFunction<>(forceContainer);

        // Plate (patches 0, 1, 4, 5)
        // tp = 10
        // Girder web (patches 0, 1, 4, 5)
        // tgw = 15
        // Girder flange (patches 0, 1, 4, 5)
        // tgf = 30
        // Stiffener web (patches 0, 1, 4, 5)
        // tsw = 10
        // Stiffener flange (patches 0, 1, 4, 5)
        // tsf = 15

        FunctionContainer tContainer(mp.nPatches());

        gsVector<> tmp2(1);
        gsFunctionExpr<> tp(std::to_string(15.), 3);
        gsFunctionExpr<> tgw(std::to_string(1.), 3);
        gsFunctionExpr<> tgf(std::to_string(2.), 3);
        gsFunctionExpr<> tsw(std::to_string(1.), 3);
        gsFunctionExpr<> tsf(std::to_string(1.), 3);

        gsDebugVar(tsf.targetDim());

        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            tContainer.at(*it) = &tp;
        // Girder webs
        for (typename std::vector<index_t>::iterator it = girderWebIDs.begin(); it!=girderWebIDs.end(); it++)
            tContainer.at(*it) = &tgw;
        // Girder flanges
        for (typename std::vector<index_t>::iterator it = girderFlangeIDs.begin(); it!=girderFlangeIDs.end(); it++)
            tContainer.at(*it) = &tgf;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            tContainer.at(*it) = &tsw;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            tContainer.at(*it) = &tsf;

        for (size_t p = 0; p != tContainer.size(); p++)
            t.addPiece(*tContainer.at(p));

        // t = gsPiecewiseFunction<>(tContainer);
    }
    else
        GISMO_ERROR("Test case not known");
    //! [Set boundary conditions]


    //! [Make material functions]
    // Linear isotropic material model
    gsConstantFunction<> pressFun(pressure,3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    // gsConstantFunction<> t(thickness,3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> nu(PoissonRatio,3);
    // gsFunctionExpr<> force("0","0","0",3);

    // Linear anisotropic material model (only one layer for example purposes)
    index_t kmax = 1; // number of layers
    std::vector<gsFunctionSet<> * > Gs(kmax); // Material matrices
    std::vector<gsFunctionSet<> * > Ts(kmax); // Thickness per layer
    std::vector<gsFunctionSet<> * > Phis(kmax); // Fiber angle per layer

    // Make material matrix
    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = &Gfun;

    // Define fiber angle
    gsConstantFunction<> phi;
    phi.setValue(0,3);
    Phis[0] = &phi;

    // Define thickness
    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = &thicks;

    //! [Make assembler]
    std::vector<gsFunction<>*> parameters;
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    // Make gsMaterialMatrix depending on the user-defined choices
    if (composite)
    {
        materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
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

    gsMaterialMatrixContainer<real_t> materialMats(mp.nPatches());
    for (size_t p = 0; p!=mp.nPatches(); p++)
        materialMats.add(materialMatrix);


    // gsMaterialMatrixContainer<real_t> materialMatsSingle(mp.nPatches());
    // gsMaterialMatrixBase<real_t> * mmtmp;
    // for (size_t p = 0; p!=mp.nPatches(); p++)
    // {
    //     parameters.resize(2);
    //     parameters[0] = const_cast<gsFunction<> *>(&(E.function(p)));
    //     parameters[1] = const_cast<gsFunction<> *>(&(nu.function(p)));
    //     options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    //     options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    //     mmtmp = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    //     materialMatsSingle.add(mmtmp);
    // }



    // gsDebugVar(materialMats);
    // gsDebugVar(materialMatsSingle);

    // gsMatrix<> z(1,1);
    // z.setZero();
    // // gsMaterialMatrixEval2<real_t,MaterialOutput::VectorN> vectorN(materialMats,&mp,z);

    // gsMaterialMatrixEval<real_t,MaterialOutput::Thickness> Thickness(materialMatrix,&mp,z);
    // gsMaterialMatrixEval<real_t,MaterialOutput::Parameters> Parameters(materialMatrix,&mp,z);


    // gsMatrix<> result;
    // gsVector<> pt(2); pt.setConstant(0.25);

    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    // {
    //     gsDebug<<"-----------Patch "<<p<<"\n";
    //     Thickness.piece(p).eval_into(pt,result);
    //     gsDebug<<"Thickness: "<<result.transpose()<<"\n";

    //     Parameters.piece(p).eval_into(pt,result);
    //     gsDebug<<"Parameters: "<<result.transpose()<<"\n";
    // }


    // gsMaterialMatrixEval<real_t,MaterialOutput::Thickness> Thickness2(materialMatsSingle,&mp,z);
    // gsMaterialMatrixEval<real_t,MaterialOutput::Parameters> Parameters2(materialMatsSingle,&mp,z);


    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    // {
    //     gsDebug<<"-----------Patch "<<p<<"\n";
    //     Thickness2.piece(p).eval_into(pt,result);
    //     gsDebug<<"Thickness: "<<result.transpose()<<"\n";

    //     Parameters2.piece(p).eval_into(pt,result);
    //     gsDebug<<"Parameters: "<<result.transpose()<<"\n";
    // }

    // Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t>* assembler;
    // if(membrane) // no bending term
    //     assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,bc,force,materialMats);
    // else
    //     assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMats);

    if(membrane) // no bending term
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,bc,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);


    // Set the penalty parameter for the interface C1 continuity
    assembler->options().setReal("IfcDirichlet",ifcDirichlet);
    assembler->options().setReal("IfcClamped",ifcClamped);
    assembler->addWeakC0(mp.topology().interfaces());
    assembler->addWeakC1(mp.topology().interfaces());
    assembler->initInterfaces();

    assembler->setPointLoads(pLoads);
    if (pressure!= 0.0)
        assembler->setPressure(pressFun);
    //! [Make assembler]

    // Set stopwatch
    gsStopwatch stopwatch,stopwatch2;
    real_t time = 0.0;
    real_t totaltime = 0.0;

    //! [Define jacobian and residual]
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      ThinShellAssemblerStatus status;
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
      GISMO_ASSERT(status==ThinShellAssemblerStatus::Success,"Assembly failed");
      time += stopwatch.stop();
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      ThinShellAssemblerStatus status;
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleVector(mp_def);
      GISMO_ASSERT(status==ThinShellAssemblerStatus::Success,"Assembly failed");
      time += stopwatch.stop();
      return assembler->rhs();
    };
    //! [Define jacobian and residual]

    stopwatch.restart();
    stopwatch2.restart();
    assembler->assemble();
    time += stopwatch.stop();

    //! [Assemble linear part]
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler->numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]


    //! [Solve non-linear problem]
    if (nonlinear)
    {
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
        }
    }
    //! [Solve non-linear problem]

    totaltime += stopwatch2.stop();

    //! [Construct and evaluate solution]
    mp_def = assembler->constructSolution(solVector);
    gsMultiPatch<> deformation = assembler->constructDisplacement(solVector);
    //! [Construct and evaluate solution]


    //! [Construct and evaluate solution]
    gsVector<> refVals = deformation.patch(refPatch).eval(refPoint);
    // real_t numVal;
    // if      (testCase == 0 || testCase == 1 || testCase == 3)
    //     numVal = refVals.at(2);
    // else
    //     numVal = refVals.at(1);

    // gsInfo << "Displacement at reference point: "<<numVal<<"\n";
    gsInfo << "Displacement at reference point: "<<refVals<<"\n";
    //! [Construct and evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        // gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( solField, "Deformation", 1000, true);
        gsWriteParaview<>( solField, "Deformation", 1000, false);

        if (testCase==3)
            gsWarn<<"Deformations are possibly zero in Paraview, due to the default precision (1e-5).\n";
    }
    if (stress)
    {
        if (testCase==2)
        {
            gsWarn<<"Stresses cannot be plotted for this case due to the singularity at the top of the geometry.\n";
        }
        else
        {
            gsPiecewiseFunction<> membraneStresses;
            assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
            gsField<> membraneStress(mp_def,membraneStresses, true);

            gsPiecewiseFunction<> flexuralStresses;
            assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
            gsField<> flexuralStress(mp_def,flexuralStresses, true);

            gsWriteParaview(membraneStress,"MembraneStress",1000);
            gsWriteParaview(flexuralStress,"FlexuralStress",1000);
        }
    }
    // ! [Export visualization in ParaView]

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;

}// end main

template <class T>
gsMultiPatch<T> Plate(T Lp, T Wp, T x, T y, T z)
{
    gsMultiPatch<T> result, tmp;

    // Web
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< Lp,0,0;
    result.patch(0).coefs().row(2)<< 0,Wp,0;
    result.patch(0).coefs().row(3)<< Lp,Wp,0;

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> Strip(T Lb, T Hw, T x, T y, T z)
{
    gsMultiPatch<T> result, tmp;

    // Web
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< 0,Lb,0;
    result.patch(0).coefs().row(2)<< 0,0,Hw;
    result.patch(0).coefs().row(3)<< 0,Lb,Hw;

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> IBeam(T Lb, T Hw, T Wf, T x, T y, T z)
{
    gsMultiPatch<T> result, tmp;

    // Web
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< 0,Lb,0;
    result.patch(0).coefs().row(2)<< 0,0,Hw/2.;
    result.patch(0).coefs().row(3)<< 0,Lb,Hw/2.;

    // Flange, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< 0,0,Hw/2.;
    result.patch(1).coefs().row(1)<< Wf/2,0,Hw/2.;
    result.patch(1).coefs().row(2)<< 0,Lb,Hw;
    result.patch(1).coefs().row(3)<< Wf/2,Lb,Hw/2.;

    // Flange, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(2).embed(3);
    result.patch(2).coefs().row(0)<< -Wf/2,0,Hw/2.;
    result.patch(2).coefs().row(1)<< 0,0,Hw/2.;
    result.patch(2).coefs().row(2)<< -Wf/2,Lb,Hw/2.;
    result.patch(2).coefs().row(3)<< 0,Lb,Hw/2.;

    // Flange, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(3).embed(3);
    result.patch(3).coefs().row(0)<< 0,0,-Hw/2.;
    result.patch(3).coefs().row(1)<< Wf/2,0,-Hw/2.;
    result.patch(3).coefs().row(2)<< 0,Lb,-Hw;
    result.patch(3).coefs().row(3)<< Wf/2,Lb,-Hw/2.;

    // Flange, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(4).embed(3);
    result.patch(4).coefs().row(0)<< -Wf/2,0,-Hw/2.;
    result.patch(4).coefs().row(1)<< 0,0,-Hw/2.;
    result.patch(4).coefs().row(2)<< -Wf/2,Lb,-Hw/2.;
    result.patch(4).coefs().row(3)<< 0,Lb,-Hw/2.;

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.computeTopology();

    // GISMO_ERROR("Interfaces not yet configured");

    // result.addInterface(&result.patch(0),4,&result.patch(1),1);
    // result.addInterface(&result.patch(0),4,&result.patch(2),2);

    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> TBeam(T Lb, T Hw, T Wf, T x, T y, T z)
{
    gsMultiPatch<T> result, tmp;

    // Web
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< 0,Lb,0;
    result.patch(0).coefs().row(2)<< 0,0,Hw;
    result.patch(0).coefs().row(3)<< 0,Lb,Hw;

    // Flange, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< 0,0,Hw;
    result.patch(1).coefs().row(1)<< Wf/2,0,Hw;
    result.patch(1).coefs().row(2)<< 0,Lb,Hw;
    result.patch(1).coefs().row(3)<< Wf/2,Lb,Hw;

    // Flange, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(2).embed(3);
    result.patch(2).coefs().row(0)<< -Wf/2,0,Hw;
    result.patch(2).coefs().row(1)<< 0,0,Hw;
    result.patch(2).coefs().row(2)<< -Wf/2,Lb,Hw;
    result.patch(2).coefs().row(3)<< 0,Lb,Hw;

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.computeTopology();

    // result.addInterface(&result.patch(0),4,&result.patch(1),1);
    // result.addInterface(&result.patch(0),4,&result.patch(2),2);
    // result.addInterface(&result.patch(1),1,&result.patch(2),2);

    result.addAutoBoundaries();
    gsWrite(result,"result");

    return result;
}

template <class T>
gsMultiPatch<T> LBeam(T Lb, T Hw, T Wf, T x, T y, T z)
{
    gsMultiPatch<T> result, tmp;

    // Web
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< 0,Lb,0;
    result.patch(0).coefs().row(2)<< 0,0,Hw;
    result.patch(0).coefs().row(3)<< 0,Lb,Hw;

    // Flange, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< 0,0,Hw;
    result.patch(1).coefs().row(1)<< Wf,0,Hw;
    result.patch(1).coefs().row(2)<< 0,Lb,Hw;
    result.patch(1).coefs().row(3)<< Wf,Lb,Hw;

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.addInterface(&result.patch(0),4,&result.patch(1),1);

    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> PanelT(T Lp, T Wp, T Hw, T Wf, T x, T y, T z)

{
    gsMultiPatch<T> result, tmp;

    // Base plate, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< Wp/2,0,0;
    result.patch(0).coefs().row(2)<< 0,Lp,0;
    result.patch(0).coefs().row(3)<< Wp/2,Lp,0;

    // Base plate, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< -Wp/2,0,0;
    result.patch(1).coefs().row(1)<< 0,0,0;
    result.patch(1).coefs().row(2)<< -Wp/2,Lp,0;
    result.patch(1).coefs().row(3)<< 0,Lp,0;

    // T-Beam
    gsMultiPatch<> beam = TBeam(Lp,Hw,Wf);

    for (size_t p=0; p!=beam.nPatches(); p++)
        result.addPatch(beam.patch(p));

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    // result.addInterface(&result.patch(1),2,&result.patch(0),1);
    // result.addInterface(&result.patch(2),3,&result.patch(0),1);
    // result.addInterface(&result.patch(2),4,&result.patch(3),1);
    // result.addInterface(&result.patch(2),4,&result.patch(4),2);

    result.computeTopology();
    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> PanelL(T Lp, T Wp, T Hw, T Wf, T x, T y, T z)

{
    gsMultiPatch<T> result, tmp;

    // Base plate, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< Wp/2,0,0;
    result.patch(0).coefs().row(2)<< 0,Lp,0;
    result.patch(0).coefs().row(3)<< Wp/2,Lp,0;

    // Base plate, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< -Wp/2,0,0;
    result.patch(1).coefs().row(1)<< 0,0,0;
    result.patch(1).coefs().row(2)<< -Wp/2,Lp,0;
    result.patch(1).coefs().row(3)<< 0,Lp,0;

    // T-Beam
    gsMultiPatch<> beam = LBeam(Lp,Hw,Wf);

    for (size_t p=0; p!=beam.nPatches(); p++)
        result.addPatch(beam.patch(p));

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    // result.addInterface(&result.patch(1),2,&result.patch(0),1);
    // result.addInterface(&result.patch(2),3,&result.patch(0),1);
    // result.addInterface(&result.patch(2),4,&result.patch(3),1);

    result.computeTopology();
    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> PlateGirderL(T PanelLength, T PanelWidth, T GirderHeight, T GirderFlangeWidth, T WebFlangeHeight, T WebFlangeWidth, T x, T y, T z)

{
    gsMultiPatch<T> result, tmp;

    // make sub panels
    std::vector<gsMultiPatch<T>> panels(14);
    panels.at(0) = Plate(PanelLength/2.,PanelWidth/2.,                  0.,                 -PanelWidth/2., 0.);
    panels.at(1) = Plate(PanelLength/2.,WebFlangeWidth,                 0.,                 0.,             0.);
    panels.at(2) = Plate(PanelLength/2.,PanelWidth/2.-WebFlangeWidth,   0.,                 WebFlangeWidth, 0.);
    panels.at(3) = Plate(PanelLength/2.,PanelWidth/2.,                  -PanelLength/2.,    -PanelWidth/2., 0.);
    panels.at(4) = Plate(PanelLength/2.,WebFlangeWidth,                 -PanelLength/2.,    0.,             0.);
    panels.at(5) = Plate(PanelLength/2.,PanelWidth/2.-WebFlangeWidth,   -PanelLength/2.,    WebFlangeWidth, 0.);

    panels.at(6) = LBeam(PanelLength/2.,WebFlangeHeight,WebFlangeWidth);
    panels.at(7) = LBeam(PanelLength/2.,WebFlangeHeight,WebFlangeWidth);

    for (size_t p = 0; p!=panels[6].nPatches(); p++)
        panels[6].patch(p).coefs().col(0).swap(panels[6].patch(p).coefs().col(1));
    for (size_t p = 0; p!=panels[7].nPatches(); p++)
    {
        panels[7].patch(p).coefs().col(0).swap(panels[7].patch(p).coefs().col(1));
        panels[7].patch(p).coefs().col(0).array() -= PanelLength / 2.;
    }

    panels.at(8) = TBeam(WebFlangeWidth,              GirderHeight,GirderFlangeWidth,0.,0.,WebFlangeHeight);
    panels.at(9) = TBeam(PanelWidth/2.-WebFlangeWidth,GirderHeight,GirderFlangeWidth,0.,WebFlangeWidth,WebFlangeHeight);
    panels.at(10) = TBeam(PanelWidth/2.,               GirderHeight,GirderFlangeWidth,0.,-PanelWidth/2.,WebFlangeHeight);

    panels.at(11) = Strip(WebFlangeWidth,              WebFlangeHeight,0.,0.);
    panels.at(12) = Strip(PanelWidth/2.-WebFlangeWidth,WebFlangeHeight,0.,WebFlangeWidth);
    panels.at(13) = Strip(PanelWidth/2.,               WebFlangeHeight,0.,-PanelWidth/2.);

    for (typename std::vector<gsMultiPatch<T>>::iterator it = panels.begin(); it!=panels.end(); it++)
        for (size_t p = 0; p!=it->nPatches(); p++)
           result.addPatch(it->patch(p));

    result.computeTopology();
    result.addAutoBoundaries();

    return result;
}
