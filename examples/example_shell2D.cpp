/** @file example_shell2D.cpp

    @brief Simple 2D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixEval.h>

#include <gsKLShell/gsMaterialMatrixTFT.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>

#ifdef GISMO_WITH_IPOPT
#include <gsIpopt/gsOptProblem.h>
#endif
//#include <gsThinShell/gsNewtonIterator.h>


using namespace gismo;

/**
 * @brief
 * Simple optimization example, to demonstrate the definition of an
 * optimization problem using the base class gsOptProblem.
 *
 *  This class implements the following NLP.
 *
 * min_x f(x) = (... )^2     (objetive function) --> square minimize it to zero
 *  s.t.
 *       beta > 0
 *       t sigma t > 0
 *       0.5 * pi <= theta <= 0.5 * pi
 *
 */
#ifdef GISMO_WITH_IPOPT

template <typename T>
class gsOptProblemExample : public gsOptProblem<T>
{
public:

    gsOptProblemExample(const gsMatrix<T> & C, const gsMatrix<T> & e)
    :
    m_C(C),
    m_e(e)
    {
        m_numDesignVars  = 1;
        m_numConstraints = 1;
        m_numConJacNonZero = 1;

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        // theta has a lower bound of -1 and an upper bound of 1
    const T inf = std::numeric_limits<T>::infinity();
        m_desLowerBounds[0] = -0.5*M_PI;
        m_desUpperBounds[0] =  0.5*M_PI;

        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);

        // we have one equality constraint, so we set the bounds on
        // this constraint to be equal (and zero).
        m_conLowerBounds[0] = 0;
        m_conUpperBounds[0] = inf;

        // m_conLowerBounds[0] = -1;
        // m_conUpperBounds[0] = 1;

        // we initialize x in bounds, in the upper right quadrant
        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign(0,0) = 0.0;

        //
        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);
    }


public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        T theta = u(0,0);
        T n1 = math::cos(theta);
        T n2 = math::sin(theta);
        T m1 = -math::sin(theta);
        T m2 = math::cos(theta);

        gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<T,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<T,1,1> gamma = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );

        gsMatrix<T,1,1> result = n2_vec.transpose() * m_C * m_e + gamma * n2_vec.transpose() * m_C * n1_vec;
        GISMO_ASSERT(result.rows()==1 && result.cols() ==1,"f is not scalar!");

        return result(0,0) * result(0,0);
    }

    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        T theta = u(0,0);
        result.resize(m_numConstraints,1);
        // return the value of the constraints: g(x)
        real_t n1 = math::cos(theta);
        real_t n2 = math::sin(theta);
        real_t m1 = -math::sin(theta);
        real_t m2 = math::cos(theta);

        gsVector<real_t,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<real_t,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<real_t,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<real_t,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<real_t,1,1> gamma = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );
        result(0) = gamma(0,0);
    }

    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // at the moment only a full finite difference matrix is returned.

        gsVector<T> uu = u;
        gsVector<T> e1( this->m_numConstraints );
        gsVector<T> e2( this->m_numConstraints );
        gsAsVector<T> ee1( e1.data() , e1.rows() );
        gsAsVector<T> ee2( e2.data() , e2.rows() );

        index_t lastDesginVar = -1;

        // TODO: Replace by a better value or use AD...
        const T h = T(0.00001);

        for( index_t i = 0 ; i < this->m_numConJacNonZero; ++i )
        {
            index_t row = this->m_conJacRows[i];  // constrains
            index_t col = this->m_conJacCols[i];  // designVariables

            if( lastDesginVar != col )
            {
                gsAsConstVector<T> uuMap( uu.data() , uu.rows() );

                uu(col) -= h/2.;
                evalCon_into( uuMap, ee1);
                uu(col) += h;
                evalCon_into( uuMap, ee2 );
                uu(col) = u(col);

                lastDesginVar = col;
            }

            result(i) = (0.5*e1(row) - 0.5*e2(row)) / h;

        }

    }

private:

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;

    const gsMatrix<T> m_C;
    const gsMatrix<T> m_e;
};
#endif

template<class T>
T findMod(T a, T b)
{
    T mod;
    // Handling negative values
    if (a < 0)
        mod = -a;
    else
        mod =  a;
    if (b < 0)
        b = -b;

    // Finding mod by repeated subtraction

    while (mod >= b)
        mod = mod - b;

    // Sign of result typically depends
    // on sign of a.
    if (a < 0)
        return -mod;

    return mod;
}

template <class T>
class gsScalarRootFinder
{
protected:
    T m_x;
    T m_xmin;
    T m_xmax;
    T m_error;
    T m_tolerance;

    index_t m_maxIterations;
    index_t m_iteration;

    typedef std::function < T ( T const &) > function_t;
    function_t m_function;

public:
    gsScalarRootFinder(T xmin, T xmax, function_t fun)
    :
    m_xmin(xmin),
    m_xmax(xmax),
    m_function(fun),
    m_iteration(0),
    m_maxIterations(15),
    m_error(1.0),
    m_tolerance(1e-15)
    {

    }

    void compute()
    {
        // From: https://lemesurierb.github.io/elementary-numerical-analysis-python/notebooks/root-finding-without-derivatives-python.html
        T x_older = m_xmin;
        T x_more_recent = m_xmax;
        T x_new;
        T f_x_older = m_function(x_older);
        T f_x_more_recent = m_function(x_more_recent);
        T f_x_new;
        for (m_iteration = 0; m_iteration <= m_maxIterations; m_iteration++)
        {
            x_new = (x_older * f_x_more_recent - f_x_older * x_more_recent)/(f_x_more_recent - f_x_older);
            x_new = findMod(x_new , M_PI);

            f_x_new = m_function(x_new);

            gsDebugVar(x_new);
            gsDebugVar(x_more_recent);

            x_older = x_more_recent;
            x_more_recent = x_new;
            f_x_older = f_x_more_recent;
            f_x_more_recent = f_x_new;

            m_error = std::abs(x_older - x_more_recent);

            gsDebug<<m_iteration<<"\t"<<m_error<<"\n";

            if (m_error < m_tolerance)
                break;
        }
        if (m_iteration < m_maxIterations)
        {
            m_x = x_new;
        }
        else
        {
            gsWarn<<"Did not converge in "<<m_iteration<<"iterations\n";
            m_x = NAN;
        }


        // real_t a = m_xmin;
        // real_t b = m_xmax;
        // real_t fa = m_function(a);
        // real_t fb = m_function(b);
        // real_t c;
        // for (m_iteration = 0; m_iteration!= m_maxIterations; m_iteration++)
        // {
        //         c = (a * fb - fa * b)/(fb - fa);
        //         real_t fc = m_function(c);
        //         if (fa * fc < 0)
        //         {
        //             b = c;
        //             fb = fc;
        //         }
        //         else
        //         {
        //             a = c;
        //             fa = fc;
        //         }

        //         m_error = std::abs(b-a);

        //         gsDebug<<m_iteration<<"\t"<<m_error<<"\n";

        //     if (m_error < m_tolerance)
        //         break;
        // }
        // if (m_iteration < m_maxIterations)
        // {
        //     m_x = c;
        // }
        // else
        // {
        //     gsWarn<<"Did not converge in "<<m_iteration<<"iterations\n";
        //     m_x = NAN;
        // }

    }

    T result() { return m_x; }
    T error()  { return m_error; }
};

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
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
    bool membrane = false;

    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

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
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("comp", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("composite", "Composite material", composite);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    real_t length = 1;
    real_t width = 1;
    bool nonlinear = true;
    if (testCase == 0)
    {
        real_t mu = 1.5e6;
        thickness = 0.001;
        if (!Compressibility)
          PoissonRatio = 0.499;
        else
          PoissonRatio = 0.45;
        E_modulus = 2*mu*(1+PoissonRatio);

        length= 1;
        width = 1;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= length;
        mp.patch(0).coefs().col(1) *= width;
        mp.addAutoBoundaries();
    }
    else if (testCase == 1)
    {
        E_modulus = 1;
        thickness = 1;
        if (!Compressibility)
          PoissonRatio = 0.499;
        else
          PoissonRatio = 0.45;

      length = 2;
      width = 1;

      mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
      mp.patch(0).coefs().col(0) *= length;
      mp.patch(0).coefs().col(1) *= width;
      mp.addAutoBoundaries();
    }
    else if (testCase == 2)
    {
        nonlinear = false;
        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0.3;

        length = 3;
        width = 1;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= length;
        mp.patch(0).coefs().col(1) *= width;
        mp.patch(0).coefs()(0,0) = width;

    }
    else if (testCase == 3)
    {
        nonlinear = true;
        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0.3;

        length = 3;
        width = 1;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= length;
        mp.patch(0).coefs().col(1) *= width;
    }
    else if (testCase == 4)
    {
        mp.addPatch( gsNurbsCreator<>::NurbsAnnulus(1,2) ); // degree
        mp.patch(0).coefs().col(0) *= length;
        mp.patch(0).coefs().col(1) *= width;

        mp.addInterface(0,boundary::east, 0, boundary::west);

    }


    gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";
    gsDebugVar(PoissonRatio);

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true,true);
    gsWrite(mp_def,"mp");

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(2);
    tmp << 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    real_t pressure = 0.0;


    std::string fx = "0";
    std::string fy = "0";

    if (testCase == 0)
        fx = "2625";
    else if (testCase == 2)
    {
        real_t sigmax = 1e-1;
        char buffer[2000];
        sprintf(buffer,"%e ( 1 - y/%e)",sigmax,width);
        fx = buffer;
    }
    else if (testCase == 3)
        fx = "0.1";

    gsFunctionExpr<> neuData(fx,fy,2);

    if (testCase == 0) // Uniaxial tension; use with hyperelastic material model!
    {
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
    }
    else if (testCase == 1)
    {
        for (index_t i=0; i!=2; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
        }

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (2); load << 0.25, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 2)
    {
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );

        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
    }
    else if (testCase == 3)
    {
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );

        bc.addCondition(boundary::north, condition_type::collapsed, 0, 0 ,false, 0);
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );

        gsVector<> point(2); point<< 0.5,1.0;
        gsVector<> load (2); load << 0.25, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 4)
    {
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false, 0);
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );

        gsVector<> point(2);
        gsVector<> load (2);

        point<< 0.0,0.0;
        load << 0.0, 1.0 ;
        pLoads.addLoad(point, load, 0 );

        point<< 0.25,0.0;
        load << -1.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );

        point<< 0.5,0.0;
        load << 0.0,-1.0 ;
        pLoads.addLoad(point, load, 0 );

        point<< 0.75,0.0;
        load << 1.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else
        GISMO_ERROR("Test case not known");

    //! [Refinement]

    // Linear isotropic material model and Neo-Hookean material
    gsConstantFunction<> force(tmp,2);
    gsConstantFunction<> pressFun(pressure,2);
    gsFunctionExpr<> t(std::to_string(thickness),2);
    gsFunctionExpr<> E(std::to_string(E_modulus),2);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
    gsFunctionExpr<> rho(std::to_string(Density),2);

    // Mooney-Rivlin material
    gsConstantFunction<> ratio(Ratio,2);

    // Ogden material
    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,2);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,2);
    gsConstantFunction<> alpha2(5.0,2);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,2);
    gsConstantFunction<> alpha3(-2.0,2);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,2);
    gsMaterialMatrixBase<real_t>* materialMatrix;

    // Linear anisotropic material model
    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,2);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,2);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,2);
    Ts[0] = &thicks;

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK & Composites
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
    else if (material==4) // OG
    {
      parameters.resize(8);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
      parameters[4] = &mu2;
      parameters[5] = &alpha2;
      parameters[6] = &mu3;
      parameters[7] = &alpha3;
    }

    gsOptionList options;
    if      (material==0)
    {
        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<2,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,materialMatrix);

    assembler->setPointLoads(pLoads);

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

    // Solve linear problem
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);

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

    totaltime += stopwatch2.stop();

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);
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

        gsPiecewiseFunction<> pstress_m;
        assembler->constructStress(mp_def,pstress_m,stress_type::principal_stress_membrane);
        gsField<> pstressM(mp_def,pstress_m, true);

        gsPiecewiseFunction<> pstress_f;
        assembler->constructStress(mp_def,pstress_f,stress_type::principal_stress_flexural);
        gsField<> pstressF(mp_def,pstress_f, true);

        gsPiecewiseFunction<> stretch1;
        assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
        gsField<> stretchDir1(mp_def,stretch1, true);

        gsPiecewiseFunction<> stretch2;
        assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
        gsField<> stretchDir2(mp_def,stretch2, true);

        gsPiecewiseFunction<> stretch3;
        assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
        gsField<> stretchDir3(mp_def,stretch3, true);

        gsPiecewiseFunction<> VMStresses;
        assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
        gsField<> VMStress(mp_def,VMStresses, true);



        gsWriteParaview(membraneStress,"MembraneStress",5000);
        gsWriteParaview(VMStress,"MembraneStressVM",5000);
        gsWriteParaview(flexuralStress,"FlexuralStress",5000);
        gsWriteParaview(Stretches,"PrincipalStretch",5000);
        gsWriteParaview(pstressM,"PrincipalMembraneStress",5000);
        gsWriteParaview(pstressF,"PrincipalFlexuralStress",5000);
        gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
        gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
        gsWriteParaview(stretchDir2,"PrincipalDirection2",5000);
        gsWriteParaview(stretchDir3,"PrincipalDirection3",5000);


    }

    if (testCase==2)
    {
        gsPiecewiseFunction<> VMStresses;
        assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
        gsField<> VMStress(mp_def,VMStresses, true);
        gsVector<> pt(2);
        pt<<0,0;
        gsInfo<<"Stress in corner: "<<VMStresses.piece(0).eval(pt)<<"\n";
    }

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////



    gsMaterialMatrixLinear<2,real_t> * materialMatrixLinear = new gsMaterialMatrixLinear<2,real_t>(mp,t,parameters,rho);
    gsMaterialMatrixTFT<2,real_t> * materialMatrixTFT = new gsMaterialMatrixTFT<2,real_t>(materialMatrixLinear);

    gsMatrix<> z(1,1);
    z.setZero();
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> matA(materialMatrixTFT,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> matB(materialMatrixTFT,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> matC(materialMatrixTFT,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> matD(materialMatrixTFT,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> vecN(materialMatrixTFT,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> vecM(materialMatrixTFT,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::PStrainN> pstrain(materialMatrix,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::PStressN> pstress(materialMatrix,mp_def,z);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> matAL(materialMatrixLinear,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixB> matBL(materialMatrixLinear,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixC> matCL(materialMatrixLinear,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixD> matDL(materialMatrixLinear,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> vecNL(materialMatrixLinear,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::VectorM> vecML(materialMatrixLinear,mp_def,z);


    gsMaterialMatrixEval<real_t,MaterialOutput::TensionField> tensionfield(materialMatrixLinear,mp_def,z);


    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> mat(materialMatrixLinear,mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir> pdir(materialMatrixLinear,mp_def,z);



    gsField<> tensionField(mp_def, tensionfield, true);
    gsInfo<<"Plotting in Paraview...\n";
    gsWriteParaview<>( tensionField, "TensionField", 10000, true);


    gsVector<> pt(2);
    pt<<0.5,0.5;

    gsMatrix<> result;

    matA.eval_into(pt,result);
    gsMatrix<real_t,3,3> Amat = result.reshape(3,3);
    gsDebugVar(Amat);

    matB.eval_into(pt,result);
    gsMatrix<real_t,3,3> Bmat = result.reshape(3,3);
    gsDebugVar(Bmat);

    matC.eval_into(pt,result);
    gsMatrix<real_t,3,3> Cmat = result.reshape(3,3);
    gsDebugVar(Cmat);

    matD.eval_into(pt,result);
    gsMatrix<real_t,3,3> Dmat = result.reshape(3,3);
    gsDebugVar(Dmat);

    vecN.eval_into(pt,result);
    gsMatrix<real_t,3,1> Nvec = result;
    gsDebugVar(Nvec);

    vecM.eval_into(pt,result);
    gsMatrix<real_t,3,1> Mvec = result;
    gsDebugVar(Mvec);

///////////////////////

    matAL.eval_into(pt,result);
    gsMatrix<real_t,3,3> AmatL = result.reshape(3,3);
    gsDebugVar(AmatL);

    matBL.eval_into(pt,result);
    gsMatrix<real_t,3,3> BmatL = result.reshape(3,3);
    gsDebugVar(BmatL);

    matCL.eval_into(pt,result);
    gsMatrix<real_t,3,3> CmatL = result.reshape(3,3);
    gsDebugVar(CmatL);

    matDL.eval_into(pt,result);
    gsMatrix<real_t,3,3> DmatL = result.reshape(3,3);
    gsDebugVar(DmatL);

    vecNL.eval_into(pt,result);
    gsMatrix<real_t,3,1> NvecL = result;
    gsDebugVar(NvecL);

    vecML.eval_into(pt,result);
    gsMatrix<real_t,3,1> MvecL = result;
    gsDebugVar(MvecL);


// //     pstrain.eval_into(pt,result);
// //     gsDebugVar(result);


//     pdir.eval_into(pt,result);
//     gsDebugVar(result);

//     mp_def.patch(0).eval_into(pt,result);
//     gsDebugVar(result);

//     tensionfield.eval_into(pt,result);
//     gsDebugVar(result);

//     gsWarn<<"This should not be pstrain.\n";
//     pstrain.eval_into(pt,result);
//     gsMatrix<real_t,3> e = result; // THIS IS ACTUALLY S!!
//     gsDebugVar(e);

//     mat.eval_into(pt,result);
//     gsMatrix<real_t,3,3> C = result.reshape(3,3);
//     gsDebugVar(Amat);


//     typedef std::function < real_t ( real_t const &) > function_t;
//     function_t f_fun = [&C,&e](real_t const & theta)
//     {
//         real_t n1 = math::cos(theta);
//         real_t n2 = math::sin(theta);
//         real_t m1 = -math::sin(theta);
//         real_t m2 = math::cos(theta);

//         gsVector<real_t,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
//         gsVector<real_t,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
//         gsVector<real_t,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
//         gsVector<real_t,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

//         gsMatrix<real_t,1,1> gamma = - ( n1_vec.transpose() * C * e ) / ( n1_vec.transpose() * C * n1_vec );
//         gsMatrix<real_t,1,1> result = n2_vec.transpose() * C * e + gamma * n2_vec.transpose() * C * n1_vec;
//         GISMO_ASSERT(result.rows()==1 && result.cols() ==1,"f is not scalar!");
//         return result(0,0);
//     };

//     typedef std::function < real_t ( real_t const &) > function_t;
//     function_t gamma_fun = [&C,&e](real_t const & theta)
//     {
//         real_t n1 = math::cos(theta);
//         real_t n2 = math::sin(theta);
//         real_t m1 = -math::sin(theta);
//         real_t m2 = math::cos(theta);

//         gsVector<real_t,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
//         gsVector<real_t,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
//         gsVector<real_t,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
//         gsVector<real_t,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

//         gsMatrix<real_t,1,1> gamma = - ( n1_vec.transpose() * C * e ) / ( n1_vec.transpose() * C * n1_vec );
//         return gamma(0,0);
//     };


//     gsVector<> x = gsVector<>::LinSpaced(100,-0.5 * M_PI,0.5 * M_PI);

// #ifdef GISMO_WITH_IPOPT
//     gsOptProblemExample<real_t> opt(C,e);

//     for (index_t k = 0; k!=x.size(); k++)
//     {
//         std::vector<real_t> X,dX,C,dC;
//         X.push_back(x[k]);
//         dX.resize(1);
//         gsAsVector<> der(dX);

//         C.resize(1);
//         gsAsVector<> con(C);

//         dC.resize(1);
//         gsAsVector<> dcon(dC);


//         opt.gradObj_into(gsAsConstVector<>(X),der);
//         opt.evalCon_into(gsAsConstVector<>(X),con);
//         opt.jacobCon_into(gsAsConstVector<>(X),dcon);
//         gsInfo<<x[k]<<","<<opt.evalObj(gsAsConstVector<>(X))<<","<<der<<","<<con<<","<<dcon<<"\n";
//     }


//     // Run optimizer
//     opt.solve();

//     // Print some details in the output
//     gsDebugVar(opt.objective());
//     gsDebugVar(opt.iterations());
//     gsDebugVar(opt.currentDesign());
//     real_t theta = opt.currentDesign()(0,0);

//     // gsScalarRootFinder<real_t> rf(0.5*M_PI,0.99*M_PI,f_fun);
//     // rf.compute();
//     // real_t theta = rf.result();
//     // gsDebugVar(rf.result());
//     // gsDebugVar(rf.error());

//     real_t n1 = math::cos(theta);
//     real_t n2 = math::sin(theta);
//     real_t m1 = -math::sin(theta);
//     real_t m2 = math::cos(theta);
//     gsVector<real_t,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
//     gsVector<real_t,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
//     gsVector<real_t,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

//     gsMatrix<real_t,1,1> denum = n1_vec.transpose() * C * n1_vec;

//     gsMatrix<real_t,3,3> C_I = C - 1 / (  n1_vec.transpose() * C * n1_vec ) * C * ( n1_vec * n1_vec.transpose() ) * C;



//     gsMatrix<real_t,1,1> gamma = - ( n1_vec.transpose() * C * e ) / ( n1_vec.transpose() * C * n1_vec );

//     gsMatrix<real_t,1,1> tmp2 = (n1_vec.transpose() * C * n2_vec);

//     GISMO_ASSERT(tmp2.rows()==1 && tmp2.cols()==1,"Must be scalar");

//     gsMatrix<real_t,1,1> df = n4_vec.transpose() * C * (e + gamma(0,0) * n1_vec)
//                             + 2 * gamma * ( n2_vec.transpose() * C * n2_vec
//                             - math::pow(tmp2(0,0),2) / (n1_vec.transpose() * C * n1_vec) );


//     gsMatrix<real_t,3,1> b = n2_vec - ( (n1_vec.transpose() * C * n2_vec)(0,0) / ( n1_vec.transpose() * C * n1_vec )(0,0)) * n1_vec;


//     gsMatrix<real_t,3,3> C_II = C_I + 2 * gamma(0,0) / df(0,0) * (C * b * b.transpose() * C);

//     gsDebugVar(C_I);
//     gsDebugVar(C_II);
//     gsDebugVar(C);

//     gsDebugVar(C_I * e);
//     gsDebugVar(C_II*e);

//     gsDebugVar(n1_vec.transpose() * C_I * e);
//     gsDebugVar(n1_vec.transpose() * C_II* e);
//     gsDebugVar(n2_vec.transpose() * C_I * e);
//     gsDebugVar(n2_vec.transpose() * C_II* e);


//     gsDebugVar(gamma);
// #else
//     gsInfo<<"GISMO_WITH_IPOPT is not enabled\n";
// #endif


    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////


    delete assembler;
    return EXIT_SUCCESS;

}// end main
