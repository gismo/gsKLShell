/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESSpline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>
#include <gsStructuralAnalysis/gsALMConsistentCrisfield.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{
// #ifndef GISMO_STRUCTURALANALYSIS
//     GISMO_ERROR("This code should be compiled with cmake flag -DGISMO_STRUCTURALANALYSIS=ON");
// #else

    bool plot       = false;
    bool stress     = false;
    bool write      = false;
    bool mesh       = false;
    bool last       = false;
    bool info       = false;
    bool writeMatrix= false;
    bool nonlinear  = false;
    bool SingularPoint = false;
    bool quasiNewton = false;
    index_t quasiNewtonInt = -1;
    index_t numRefine  = 2;
    index_t degree = 3;
    index_t smoothness = 2;
    index_t geometry = 1;
    index_t method = 0;
    std::string input;

    index_t step          = 10;
    index_t ALMmethod        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method)
    real_t dL         = 0; // Arc length
    real_t dLb        = 0.5; // Arc length to find bifurcation
    real_t Lmax       = std::numeric_limits<real_t>::max();
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e5;
    real_t tau = 1e4;
    real_t shift = -1e2;

    real_t Radius = 0.1;

    std::string fn1,fn2,fn3;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";
    std::string dirname = "ArcLengthResults";

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    index_t nmodes = 15;

    std::vector<index_t> modevec;

    index_t Compressibility = 0;
    index_t material = 3;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral


    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );

    cmd.addReal( "R", "Radius", "Relative radius compared to the width",  Radius );


    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );

    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addString( "o", "out", "Dir name of the output",  dirname );
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addInt( "m", "method", "Smoothing method to use", method );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("last", "last case only",last);
    cmd.addSwitch("writeMat", "Write projection matrix",writeMatrix);
    cmd.addSwitch( "info", "Print information", info );

    cmd.addReal("S","shift", "Shift for stability eigenvalue computation", shift);
    cmd.addReal("F","factor", "factor for bifurcation perturbation", tau);
    cmd.addInt("A","ALMmethod", "Arc length method; 1: Crisfield's method; 2: RIks' method.", ALMmethod);
    cmd.addReal("L","dLb", "arc length", dLb);
    cmd.addReal("l","dL", "arc length after bifurcation", dL);
    cmd.addReal("b","Lmax", "Maximum L", Lmax);
    cmd.addInt("n", "maxmodes", "Number of modes to be computed", nmodes);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addSwitch("bifurcation", "Compute singular points and bifurcation paths", SingularPoint);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addMultiInt("i", "modes", "Modes to select", modevec);

    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
    GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");
    if (method==3)
        GISMO_ENSURE(smoothness>=1 || smoothness <= degree-2,"Exact C1 method only works for smoothness <= p-2, but smoothness="<<smoothness<<" and p-2="<<degree-2);
    if (method==2 || method==3)
        GISMO_ENSURE(degree > 2,"Degree must be larger than 2 for the approx and exact C1 methods, but it is "<<degree);

    if (dL==0)
    {
      dL = dLb;
    }

    real_t aDim,bDim;
    real_t thickness = 0.14e-3;
    real_t E_modulus = 1;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    real_t Ratio = 7.0;

    if ((!Compressibility) && (material!=0))
      PoissonRatio = 0.5;
    else
      PoissonRatio = 0.499;

    real_t mu, C01,C10;
    if (material==3)
    {
      C10 = 6.21485502e4; // c1/2
      C01 = 15.8114570e4; // c2/2
      Ratio = C10/C01;
      mu = 2*(C01+C10);
    }
    else
    {
      C10 = 19.1010178e4;
      mu = 2*C10;
    }
    E_modulus = 2*mu*(1+PoissonRatio);
    gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"; mu = "<<mu<<"; ratio = "<<Ratio<<"\n";

    gsMultiPatch<> mp,mp_def,geom;

    bDim = 0.14; aDim = 2*bDim;

    /////////////////////////////////////////////////////////////////////////////////////
    gsMultiPatch<> mp_nurbs;
    GISMO_ASSERT(Radius<1,"Radius cannot be larger than 1");
    real_t Ri = bDim*Radius;
    real_t Ro = Ri+(bDim - Ri)/2; // middle between Ri an bDim
    gsKnotVector<> KVx (0,1,0,2) ;
    gsKnotVector<> KVy (0,1,0,3) ;
    gsMatrix<> C(6,2) ;
    C <<  Ri , 0  ,  Ro, 0
        , Ri , Ri ,  Ro, Ro
        , 0 ,  Ri ,  0 , Ro ;

    // Set weights
    gsMatrix<> ww(6,1) ;
    ww.setOnes();
    ww.at(2)= 0.707106781186548 ;
    ww.at(3)= 0.707106781186548 ;
    gsTensorNurbs<2,real_t> annulus(KVx,KVy,C,ww);
    annulus.swapDirections(0,1);
    std::vector<gsGeometry<> *> pieces = annulus.uniformSplit(0);

    gsTensorNurbs<2,real_t> bottom = static_cast<gsTensorNurbs<2,real_t> &>(*pieces[0]);
    gsTensorNurbs<2,real_t> top    = static_cast<gsTensorNurbs<2,real_t> &>(*pieces[1]);
    mp_nurbs.addPatch(bottom);
    mp_nurbs.addPatch(top);

    bottom.coefs().block(0,0,3,2)   = bottom.coefs().block(3,0,3,2);
    bottom.weights().block(0,0,3,1) = bottom.weights().block(3,0,3,1);
    bottom.coefs().block(3,0,3,1).setConstant(aDim);
    bottom.weights().block(3,0,3,1).setOnes();
    mp_nurbs.addPatch(bottom);

    top.coefs().block(0,0,3,2)   = top.coefs().block(3,0,3,2);
    top.weights().block(0,0,3,1) = top.weights().block(3,0,3,1);
    top.coefs().block(3,1,3,1).setConstant(bDim);
    top.weights().block(3,0,3,1).setOnes();
    mp_nurbs.addPatch(top);

    // we use C to store the coefs of the new patch
    C.resize(4,2);
    C.row(0) = top.coefs().row(0); // the first coefficient of top is the EV
    C.row(1)<<aDim,top.coefs()(0,1);
    C.row(2)<<top.coefs()(0,1),bDim; // the fourth coefficient of top is the second control point we need
    C.row(3)<<aDim,bDim; // the fourth coefficient of top is the second control point we need
    gsTensorBSpline<2,real_t> patch(KVx,KVx,C);
    mp_nurbs.addPatch(patch);

    // mp_nurbs.computeTopology();

    gsWriteParaview<>(mp_nurbs,"mp_nurbs",1000,true,true);
    /////////////////////////////////////////////////////////////////////////////////////

    // gsTensorBSplineBasis<2,real_t> tbasis(KVx,KVx);
    // tbasis.degreeElevate(degree-tbasis.degree(0));
    // for (int r =0; r < numRefine; ++r)
    //     tbasis.uniformRefine(1,degree-smoothness);

    // gsDebugVar(tbasis);
    // for (index_t p = 0; p != mp_nurbs.nPatches(); p++)
    // {
    //     gsQuasiInterpolate<real_t>::localIntpl(tbasis,mp_nurbs.patch(p),C);

    //     // gsL2Projection<real_t>::projectGeometry(mp_nurbs.basis(p),mp_nurbs.patch(p),C);
    //     gsDebugVar(C);
    //     // gsL2Projection<real_t>::projectGeometry(tbasis,mp_nurbs.patch(p),C);
    //     // C.resize(C.rows()/mp_nurbs.geoDim(),mp_nurbs.geoDim());
    //     gsDebugVar(tbasis);
    //     gsDebugVar(C);
    //     mp.addPatch(gsTensorBSpline<2,real_t>(tbasis,C));
    // }

    // Construct mp and project mp_nurbs on it
    gsTensorBSpline<2,real_t> * tb;
    gsTensorNurbs<2,real_t> * tn;

    for (size_t p=0; p!=mp_nurbs.nPatches(); p++)
    {
        if ((tn = dynamic_cast<gsTensorNurbs<2,real_t> *>(&mp_nurbs.patch(p))))
            mp.addPatch(gsTensorBSpline<2,real_t>(tn->basis().knots(0),tn->basis().knots(1),tn->coefs()));
        else if ((tb = dynamic_cast<gsTensorBSpline<2,real_t> *>(&mp_nurbs.patch(p))))
            mp.addPatch(*tb);
        else
            GISMO_ERROR("Cannot construct multipatch");

        tb = dynamic_cast<gsTensorBSpline<2,real_t> *>(&mp.patch(p));
        if (tb->degree(0) > tb->degree(1))
            tb->degreeElevate(1,1);
        else if (tb->degree(0) < tb->degree(1))
            tb->degreeElevate(1,0);

        if (tb->basis().maxDegree() < degree)
            tb->degreeElevate(degree-tb->basis().maxDegree());

    }

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine(1,degree-smoothness);

    gsMultiBasis<> basis(mp);


    mp.embed(3);
    mp.fixOrientation();
    mp.computeTopology();

    gsWrite(mp,"mp");

    mp_def = mp;

    gsInfo<<"alpha = "<<aDim/bDim<<"; beta = "<<bDim/thickness<<"\n";

    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<< 0.0,0.5;
    writePoints.col(1)<< 0.5,0.5;
    writePoints.col(2)<< 1.0,0.5;

    BCs.addCondition(2,boundary::north, condition_type::collapsed, 0, 0 ,false,0);
    BCs.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0 ,false,1);
    BCs.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0 ,false,2);
    BCs.addCondition(2,boundary::north, condition_type::clamped  , 0, 0, false,2);

    BCs.addCondition(4,boundary::east, condition_type::collapsed, 0, 0 ,false,0);
    BCs.addCondition(4,boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
    BCs.addCondition(4,boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
    BCs.addCondition(4,boundary::east, condition_type::clamped  , 0, 0, false,2);


    BCs.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0 ,false,0);
    BCs.addCondition(1,boundary::east, condition_type::clamped  , 0, 0, false,2);
    BCs.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0 ,false,0);
    BCs.addCondition(3,boundary::east, condition_type::clamped  , 0, 0, false,2);


    BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
    BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
    BCs.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
    BCs.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
    BCs.setGeoMap(mp);

    // ADD DIFFERENT MATERIALS!

    real_t Load = 1e0;
    gsVector<> point(2); point<< 1.0, 1.0 ;
    gsVector<> load (3); load << Load,0.0, 0.0;
    pLoads.addLoad(point, load, 4 );

    dirname = dirname + "/QuarterSheetHole_-r" + std::to_string(numRefine) + "-p" + std::to_string(degree) + "-s" + std::to_string(smoothness) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility) + "-alpha" + std::to_string(aDim/bDim) + "-beta" + std::to_string(bDim/thickness);
    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    int systemRet = system(command);
    GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

    SingularPoint = true;

    index_t cross_coordinate = 0;
    real_t cross_val = 0.0;

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);
    gsConstantFunction<> t(thickness,3);
    gsConstantFunction<> E(E_modulus,3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsConstantFunction<> rho(Density,3);
    gsConstantFunction<> ratio(Ratio,3);

    mu = E_modulus / (2 * (1 + PoissonRatio));

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
    }
    else if (material==1 || material==2) // NH & NH_ext
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==3) // MR
    {
      parameters.resize(3);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &ratio;
    }

    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    if      (material==0 && impl==1)
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    //! [Solver loop]
    gsVector<> solVector;

    gsMappedBasis<2,real_t> bb2;

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;

    gsMultiBasis<> dbasis(mp);

    if (method==-1)
    {
        // identity map
        global2local.resize(dbasis.totalSize(),dbasis.totalSize());
        for (size_t k=0; k!=dbasis.totalSize(); ++k)
            global2local.coeffRef(k,k) = 1;
        geom = mp;
        bb2.init(dbasis,global2local);
    }
    else if (method==0)
    {
        gsMPBESSpline<2,real_t> cgeom(mp,3);
        gsMappedBasis<2,real_t> basis = cgeom.getMappedBasis();

        global2local = basis.getMapper().asMatrix();
        geom = cgeom.exportToPatches();
        auto container = basis.getBasesCopy();
        dbasis = gsMultiBasis<>(container,mp.topology());
        bb2.init(dbasis,global2local);
    }
    else if (method==1)
    {
        geom = mp;
        gsDPatch<2,real_t> dpatch(geom);
        dpatch.options().setSwitch("SharpCorners",false);
        dpatch.compute();
        dpatch.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = dpatch.exportToPatches();
        dbasis = dpatch.localBasis();
        bb2.init(dbasis,global2local);
    }
    else if (method==2) // Pascal
    {
        // The approx. C1 space
        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        // approxC1.options().setSwitch("info",info);
        // approxC1.options().setSwitch("plot",plot);
        approxC1.options().setSwitch("interpolation",true);
        approxC1.options().setInt("gluingDataDegree",-1);
        approxC1.options().setInt("gluingDataSmoothness",-1);
        approxC1.update(bb2);
        geom = mp;
    }
    else if (method==3) // Andrea
    {
        dbasis = gsMultiBasis<>(mp);
        gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
        smoothC1.init();
        smoothC1.compute();

        global2local = smoothC1.getSystem();
        global2local = global2local.transpose();
        smoothC1.getMultiBasis(dbasis);
        bb2.init(dbasis,global2local);
        geom = mp;
    }
    else if (method==4)
    {
        gsAlmostC1<2,real_t> almostC1(mp);
        almostC1.compute();
        almostC1.matrix_into(global2local);

        global2local = global2local.transpose();
        geom = almostC1.exportToPatches();
        dbasis = almostC1.localBasis();
        bb2.init(dbasis,global2local);
    }
    else
        GISMO_ERROR("Option "<<method<<" for method does not exist");

    if (writeMatrix)
    {
        gsWrite(global2local,"mat");
        //gsWrite(geom,"geom");
        //gsWrite(dbasis,"dbasis");
    }

    if (plot) gsWriteParaview(geom,"geom",1000,true,false);

    gsThinShellAssembler<3, real_t, true> assembler;
    assembler = gsThinShellAssembler<3, real_t, true>(geom,dbasis,BCs,force,materialMatrix);
    assembler.options().setInt("Continuity",-1);
    assembler.options().setReal("WeakDirichlet",bcDirichlet);
    assembler.options().setReal("WeakClamped",bcClamped);
    assembler.setSpaceBasis(bb2);
    assembler.setPointLoads(pLoads);

    // Initialize the system
    // Linear
    assembler.assemble();
    gsVector<> Force = assembler.rhs();
    gsSparseMatrix<> K_L = assembler.matrix();

    // Nonlinear
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >                                     Residual_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &force) >                             ALResidual_t;
    Jacobian_t Jacobian = [&geom,&bb2,&assembler](gsVector<real_t> const &x)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleMatrix(def);
        // gsSparseMatrix<real_t> m =
        return assembler.matrix();
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&geom,&bb2,&assembler /*,&Force_const*/](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleVector(def);
        gsVector<real_t> Fint = -(assembler.rhs() - force);
        gsVector<real_t> result = Fint - lam * force/*- Force_const*/;
        return result; // - lam * force;
    };

    gsALMBase<real_t> * arcLength;
    if (ALMmethod==0)
      arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (ALMmethod==1)
      arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (ALMmethod==2)
      arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else if (ALMmethod==3)
      arcLength = new gsALMConsistentCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
      GISMO_ERROR("Method unknown");

#ifdef GISMO_WITH_PARDISO
    arcLength->options().setString("Solver","PardisoLU"); // LDLT solver
#else
    arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
#endif

    arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    arcLength->options().setReal("Length",dLb);
    arcLength->options().setInt("MaxIter",10);
    // arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
    arcLength->options().setInt("AdaptiveIterations",5);
    arcLength->options().setReal("Perturbation",tau);
    // arcLength->options().setReal("Scaling",0.0);
    arcLength->options().setReal("Tol",tol);
    arcLength->options().setReal("TolU",tolU);
    arcLength->options().setReal("TolF",tolF);
    arcLength->options().setSwitch("Verbose",true);
    if (quasiNewtonInt>0)
    {
      quasiNewton = true;
      arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength->options().setSwitch("Quasi",quasiNewton);
    arcLength->options().setReal("Shift",shift);

    gsInfo<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR
    bool bisected = false;
    real_t dLb0 = dLb;

    std::string output = "solution";
    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection membraneStressCollection(dirname + "/MembraneStress");
    gsParaviewCollection flexuralStressCollection(dirname + "/FlexuralStress");

    index_t k=0;
    real_t L=0;

    while (k < step && L < Lmax)
    // while (false)
    {
        gsInfo<<"Load step "<< k<<"; Lold = "<<L<<"\n";
        // assembler.constructSolution(solVector,solution);
        arcLength->step();

        // gsInfo<<"m_U = "<<arcLength->solutionU()<<"\n";
        if (!(arcLength->converged()))
        {
            gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
            dLb = dLb / 2.;
            arcLength->setLength(dLb);
            arcLength->setSolution(Uold,Lold);
            bisected = true;
            k -= 1;
            continue;
        }
        arcLength->computeStability(arcLength->solutionU(),quasiNewton);

        if (arcLength->stabilityChange() && SingularPoint)
        {
            gsInfo<<"Bifurcation spotted!"<<"\n";
            arcLength->computeSingularPoint(1e-4, 5, Uold, Lold, 1e-4, 0, false);
            arcLength->switchBranch();
            dLb0 = dLb = dL;
            arcLength->setLength(dLb);
            SingularPoint = false;
        }
        // else if (k == 10)
        // {
        //     dLb0 = dLb = dL;
        //     arcLength->setLength(dLb);
        // }
        indicator = arcLength->indicator();

        solVector = arcLength->solutionU();
        Uold = solVector;
        Lold = arcLength->solutionL();

        if (plot || write || stress)
        {
            /// Make a gsMappedSpline to represent the solution
            // 1. Get all the coefficients (including the ones from the eliminated BCs.)
            gsMatrix<real_t> solFull = assembler.fullSolutionVector(solVector);

            // 2. Reshape all the coefficients to a Nx3 matrix
            GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
            solFull.resize(solFull.rows()/3,3);

            // 3. Make the mapped spline
            gsMappedSpline<2,real_t> mspline(bb2,solFull);

            gsFunctionSum<real_t> def(&mp,&mspline);

            // 4. Plot the mapped spline on the original geometry
            gsField<> solField(geom, mspline,true);

            if (plot)
            {
                std::string fileName = dirname + "/" + output + util::to_string(k) + "_";
                gsWriteParaview<>(solField, fileName, 1000,mesh);
                for (index_t p = 0; p!=mp.nPatches(); p++)
                {
                    fileName = output + util::to_string(k) + "_";
                    collection.addPart(fileName + std::to_string(p) + ".vts",k,"",p);
                    if (mesh)
                        collection.addPart(fileName + "_mesh.vtp",k,"",p);
                }
            }
        }

        if (!bisected)
        {
          dLb = dLb0;
          arcLength->setLength(dLb);
        }
        bisected = false;

        L = arcLength->solutionL();
        k++;
    }

    if (plot)
    {
        collection.save();
    }
    if (stress)
    {
        membraneStressCollection.save();
        flexuralStressCollection.save();
    }

    delete arcLength;

    return EXIT_SUCCESS;
// #endif
}
