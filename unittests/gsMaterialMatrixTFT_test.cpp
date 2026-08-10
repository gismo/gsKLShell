/** @file example_shell2D.cpp

    @brief Simple 2D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

using namespace gismo;


template<class T>
class STFTfun : public gsFunction<T>
{
public:
    STFTfun(const gsMatrix<T> & gori, gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_gori(gori),
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        m_mmTFT = new gsMaterialMatrixTFT<2,T>(mm);
    }

    ~STFTfun() { delete m_mmTFT; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        for (index_t k=0; k!=u.cols(); k++)
        {
            gsVector<T> E = u.col(k);
            E(2,0) *= 0.5;
            gsMatrix<T> Cmat = 2*E+m_gori;

            gsMatrix<T> C = m_mm->eval3D_matrix_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);
            gsMatrix<T> S = m_mm->eval3D_vector_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);

            gsMatrix<T> THETA = m_mmTFT->eval_theta(C,S,E);

            T n1 = math::cos(THETA(0,0));
            T n2 = math::sin(THETA(0,0));
            T m1 = -math::sin(THETA(0,0));
            T m2 = math::cos(THETA(0,0));
            gsVector<T> n1_vec(3); n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            gsVector<T> n2_vec(3); n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            T gamma = - ( n1_vec.transpose() * S).value() / ( n1_vec.transpose() * C.reshape(3,3) * n1_vec ).value();
            gsVector<T> Ew = gamma * n1_vec;
            result.col(k) = S + C.reshape(3,3) * Ew;
        }
    }

    short_t domainDim() const {return 3;}
    short_t targetDim() const {return 3;}

private:
    const gsMatrix<T> m_gori;
    gsMaterialMatrixBase<T> * m_mm;
    gsMaterialMatrixTFT<2,T> * m_mmTFT;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};


template<class T>
class Sfun : public gsFunction<T>
{
public:
    Sfun(const gsMatrix<T> & gori, gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_gori(gori),
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
    }


    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        for (index_t k=0; k!=u.cols(); k++)
        {
            gsVector<T> E = u.col(k);
            E(2,0) *= 0.5;
            gsMatrix<T> Cmat = 2*E+m_gori;

            result.col(k) = m_mm->eval3D_vector_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);
        }
    }

    short_t domainDim() const {return 3;}
    short_t targetDim() const {return 3;}

private:
    const gsMatrix<T> m_gori;
    gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template<class T>
class Cfun : public gsFunction<T>
{
public:
    Cfun(const gsMatrix<T> & gori, gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_gori(gori),
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        m_mmTFT = new gsMaterialMatrixTFT<2,T>(mm);
    }

    ~Cfun() { delete m_mmTFT; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(9,u.cols());
        for (index_t k=0; k!=u.cols(); k++)
        {
            gsVector<T> E = u.col(k);
            E(2,0) *= 0.5;
            gsMatrix<T> Cmat = 2*E+m_gori;

            result.col(k) = m_mm->eval3D_matrix_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);
        }
    }

    short_t domainDim() const {return 3;}
    short_t targetDim() const {return 9;}

private:
    const gsMatrix<T> m_gori;
    gsMaterialMatrixBase<T> * m_mm;
    gsMaterialMatrixTFT<2,T> * m_mmTFT;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

SUITE(gsMaterialMatrixTFT_test)                 // The suite should have the same name as the file
{
    void MM_CHECK(const index_t material, const index_t impl, const bool Compressibility);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(MM_SvK)
    {
     index_t mat = 0, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }

    TEST(MM_NH_Incomp_Analytical)
    {
     index_t mat = 1, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Comp_Analytical)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Incomp_Generic)
    {
     index_t mat = 1, impl = 2;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Comp_Generic)
    {
     index_t mat = 1, impl = 2;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Incomp_Spectral)
    {
     index_t mat = 1, impl = 3;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Comp_Spectral)
    {
     index_t mat = 1, impl = 3;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }

    TEST(MM_MR_Incomp_Analytical)
    {
     index_t mat = 3, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Comp_Analytical)
    {
     index_t mat = 3, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Incomp_Generic)
    {
     index_t mat = 3, impl = 2;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Comp_Generic)
    {
     index_t mat = 3, impl = 2;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Incomp_Spectral)
    {
     index_t mat = 3, impl = 3;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Comp_Spectral)
    {
     index_t mat = 3, impl = 3;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }

    TEST(MM_OG_Incomp_Spectral)
    {
     index_t mat = 4, impl = 3;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_OG_Comp_Spectral)
    {
     index_t mat = 4, impl = 3;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }

    /*  (task 65) RUNTIME identity gate -- see the twin in gsThinShellAssembler_test.cpp.
        Twelve of the fourteen non-SvK MM_* tests were named for materials they did not
        run (every body passed mat = 1, impl = 1). The arguments are therefore not
        trusted: gsMaterialMatrixNonlinear::print() reports the object's own TEMPLATE
        parameters <matId, comp> (gsMaterialMatrixNonlinear.hpp:153-192), i.e. what
        getMaterialMatrix ACTUALLY constructed.
    */
    void CHECK_material_identity(const gsMaterialMatrixBase<real_t> & mm,
                                 const index_t material, const index_t impl,
                                 const bool Compressibility)
    {
        /*  THE EXPECTATION COMES FROM THE TEST NAME, NOT FROM THE ARGUMENTS.

            An earlier version of this helper compared the constructed object against
            the `material`/`impl` ARGUMENTS -- and it was MEASURED not to fire when two
            tests were reverted to the original fictional `mat = 1, impl = 1` (task 65,
            poison round 1: both suites stayed green). Of course: reverting the argument
            moves BOTH sides of that comparison. The defect this task exists to prevent
            is a mismatch between the test NAME and what the test RUNS, so the name is
            the only admissible anchor. With the name on one side and the constructed
            template instantiation on the other, the argument appears nowhere in the
            chain and the original defect becomes unrepresentable.
        */
        const std::string tn = UnitTest::CurrentTest::Details()->testName;

        std::string matName, implName, compName;
        if      (tn.find("_NH_")!=std::string::npos) matName = "Neo-Hookean\n";
        else if (tn.find("_MR_")!=std::string::npos) matName = "Mooney-Rivlin\n";
        else if (tn.find("_OG_")!=std::string::npos) matName = "Ogden\n";
        // the trailing newline matters: "Neo-Hookean" is a PREFIX of the NH_ext name

        if      (tn.find("_Analytical") !=std::string::npos) implName = "Analytical implementation";
        else if (tn.find("_Generic")    !=std::string::npos ||
                 tn.find("_Generalized")!=std::string::npos) implName = "Generalized implementation";
        else if (tn.find("_Spectral")   !=std::string::npos) implName = "Spectral implementation";

        // "_Incomp" must be tested BEFORE "_Comp"
        if      (tn.find("_Incomp")!=std::string::npos) compName = "\tIncompressible ";
        else if (tn.find("_Comp")  !=std::string::npos) compName = "\tCompressible ";

        if (matName.empty() || implName.empty() || compName.empty())
        {
            // A test that reaches the material path but whose name does not say which
            // material it runs is exactly the condition this gate exists to forbid.
            gsInfo << "[MATERIAL] FAIL: test name '"<<tn<<"' does not encode "
                      "material / implementation / compressibility\n";
            CHECK(false);
            return;
        }

        // gsMaterialMatrixNonlinear::print() reports the object's own TEMPLATE
        // parameters <matId, comp> (gsMaterialMatrixNonlinear.hpp:153-192), i.e. what
        // getMaterialMatrix ACTUALLY constructed -- not what was requested.
        std::ostringstream oss;
        mm.print(oss);
        const std::string s = oss.str();

        const bool okMat  = (s.find(matName)  != std::string::npos);
        const bool okImpl = (s.find(implName) != std::string::npos);
        const bool okComp = (s.find(compName) != std::string::npos);
        gsInfo << "[MATERIAL] "<<tn<<" : name wants "
               << compName.substr(1) << matName.substr(0,matName.size()-1) << " / " << implName
               << " ; args (mat "<<material<<", impl "<<impl<<", comp "<<Compressibility
               << ") CONSTRUCTED "
               << (okMat&&okImpl&&okComp ? "MATCH" : "MISMATCH") << "\n";
        if (!(okMat&&okImpl&&okComp))
            gsInfo << "[MATERIAL] constructed object reports:\n"<<s;
        CHECK(okMat);
        CHECK(okImpl);
        CHECK(okComp);
    }

    void MM_CHECK(const index_t material, const index_t impl, const bool Compressibility)
    {
        // (task 65) THE SKIP MUST RETURN -- the same defect task 61 fixed in
        // balloon_CHECK / UAT_CHECK, left behind here (task 61 open item 3). Without
        // the return, the vacuous CHECK(true) was recorded and the body ran on into
        // getMaterialMatrix.h's GISMO_ERROR for Ogden at a non-Spectral implementation.
        // Still dead code after this task: the only mat = 4 registrations are
        // MM_OG_{Incomp,Comp}_Spectral at impl = 3, so nothing reaches it (verified
        // against the dispatch in getMaterialMatrix.h:230-255, where OG x Analytical
        // and OG x Generalized are the only unsupported combinations used here).
        if (material==4 && impl!=3)
        {
            CHECK(true);
            return;
        }

        real_t E_modulus;
        real_t PoissonRatio;
        real_t thickness = 0.01;
        real_t Ratio = 7.0;

        real_t mu = 4.225e5;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        if (!Compressibility && !(material==0))
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;

        E_modulus = 2*mu*(1+PoissonRatio);

        //! [Read input file]
        gsMultiPatch<> mp;
        gsMultiPatch<> mp_def;

        real_t L = 1;
        real_t B = 1;
        real_t Delta = 0.1;
        mp_def.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp_def.patch(0).coefs().col(0) *= L;
        mp_def.patch(0).coefs().col(1) *= B;
        mp_def.addAutoBoundaries();

        mp = mp_def;
        mp.patch(0).coefs()(2,0) += Delta;
        // mp_def.patch(0).coefs()(3,0) += Delta;

        index_t numRefine  = 2;
        index_t numElevate = 3;

        // p-refine
        if (numElevate!=0)
            mp.degreeElevate(numElevate);

        // h-refine
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine();

        //! [Refinement]
        gsMultiBasis<> dbasis(mp);

        // Linear isotropic material model and Neo-Hookean material
        gsFunctionExpr<> t(std::to_string(thickness),2);
        gsFunctionExpr<> E(std::to_string(E_modulus),2);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
        // Mooney-Rivlin material
        gsConstantFunction<> ratio(Ratio,2);
        // Ogden material
        gsConstantFunction<> alpha1fun(alpha1,2);
        gsConstantFunction<> mu1fun(mu1,2);
        gsConstantFunction<> alpha2fun(alpha2,2);
        gsConstantFunction<> mu2fun(mu2,2);
        gsConstantFunction<> alpha3fun(alpha3,2);
        gsConstantFunction<> mu3fun(mu3,2);

        std::vector<gsFunctionSet<>*> parameters;
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
            parameters[2] = &mu1fun;
            parameters[3] = &alpha1fun;
            parameters[4] = &mu2fun;
            parameters[5] = &alpha2fun;
            parameters[6] = &mu3fun;
            parameters[7] = &alpha3fun;
        }

        gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
        gsMaterialMatrixBase<real_t>::uPtr materialMatrixTFT;
        gsOptionList options;
        if      (material==0)
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
            materialMatrixTFT = memory::make_unique(new gsMaterialMatrixTFT<2,real_t,true>(*materialMatrix));
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
            materialMatrixTFT = memory::make_unique(new gsMaterialMatrixTFT<2,real_t,false>(*materialMatrix));
            // (task 65) what was CONSTRUCTED, not what was asked for. Checked on the
            // BASE material matrix, not on the TFT wrapper (which has its own print).
            CHECK_material_identity(*materialMatrix,material,impl,Compressibility);
        }

        gsVector<> testpt(2);
        testpt.setConstant(0.25);
        gsMatrix<> g;

        // Compute deformation tensor
        gsMapData<> map;
        map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
        map.points = testpt;
        static_cast<const gsFunction<>&>(mp_def.patch(0)).computeMap(map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

        g = map.jacobian(0);
        g.resize(2,2);
        gsMatrix<> G(3,3);
        G.setZero();
        G.block(0,0,2,2) = g.transpose() * g;
        G(2,2) = 1;

        gsMatrix<> Cmat(3,1);
        Cmat(0,0) = G(0,0);
        Cmat(1,0) = G(1,1);
        Cmat(2,0) = G(0,1);

        static_cast<const gsFunction<>&>(mp.patch(0)).computeMap(map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
        g = map.jacobian(0);
        g.resize(2,2);
        G.setZero();
        G.block(0,0,2,2) = g.transpose() * g;
        G(2,2) = 1;

        gsMatrix<> gori(3,1);
        gori(0,0) = G(0,0);
        gori(1,0) = G(1,1);
        gori(2,0) = G(0,1);

        gsMatrix<> z(1,1); z.setZero();

        materialMatrix->setDeformed(&mp_def);
        materialMatrixTFT->setDeformed(&mp_def);

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////

        gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> matTFT(materialMatrixTFT,&mp_def,z);
        gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> vecTFT(materialMatrixTFT,&mp_def,z);

        gsVector<> pt(2);
        pt = testpt;

        gsMatrix<> resss;

        /// STRAIN
        gsMatrix<> e = 1./2. * (Cmat - gori);
        e(2,0) *= 2;
        gsMatrix<> strain_MM = materialMatrix->eval3D_strain(0,pt,z);
        CHECK_MATRIX_CLOSE(e,strain_MM,1e-3);

        /// MATRIX
        Sfun<real_t> Sfun(gori,materialMatrix.get(),0,pt,z);
        Sfun.deriv_into(e,resss);
        gsMatrix<> C_FD = resss.reshape(3,3); // Finite Differences
        gsMatrix<> C_MM = materialMatrix->eval3D_matrix(0,pt,z,MaterialOutput::Generic).reshape(3,3); // implemented
        CHECK_MATRIX_CLOSE(C_FD,C_MM,1e-3);

        /// dMATRIX
        Cfun<real_t> Cfun(gori,materialMatrix.get(),0,pt,z);
        Cfun.deriv_into(e,resss);
        gsMatrix<> dC_FD = resss.reshape(9,3);
        gsMatrix<> dC_MM = materialMatrix->eval3D_dmatrix(0,pt,z,MaterialOutput::Generic).reshape(9,3);
        dC_MM *= 2; // NOTE: dmatrix returns d(mm)/dC, and dC_FD is d(mm)/dE = 2*d(mm)/dC
        /*  (task 65) THE TOLERANCE HERE WAS SCALE-BLIND, AND RE-ARMING THE ARGUMENTS
            EXPOSED IT. This read
                CHECK_MATRIX_CLOSE(dC_FD,dC_MM,1e-3);
            and CHECK_MATRIX_CLOSE is ENTRYWISE ABSOLUTE (gsUnitTest/gs/CheckMatrix.h:37).
            The entries of dC are O(1e7) here, so "1e-3" silently demanded ~1e-10
            RELATIVE agreement between a finite-difference derivative and an analytic
            tangent. That was attainable only on the Analytical/Generalized paths the
            file used to run; it failed on all three newly-armed *_Comp_Spectral tests.

            ADJUDICATED as neither a library defect nor an oracle defect, on three
            independent measurements (task 65 report):
              1. FD STEP SWEEP. gsFunction::deriv_into (gsFunction.hpp:93-126) hardcodes
                 a 4-point stencil at h = 1e-5. Composing the function with a domain
                 rescaling to sweep the EFFECTIVE step gives, for NH/Spectral/comp,
                 max|dC_FD - dC_MM| = 3.2e-2, 7.6e-3, 1.3e-3, 1.2e-3, 4.5e-3 at
                 h = 1e-7 ... 1e-3: a textbook roundoff/truncation V-curve whose MINIMUM
                 (1.1e-3 for Ogden) already exceeds the old 1e-3 gate. No FD step could
                 have passed it. A wrong analytic tangent would instead show an
                 h-independent floor.
              2. CROSS-IMPLEMENTATION FINGERPRINT. For NH compressible the analytic
                 |dC_MM| is 14035337.358630737 (Analytical), ...358662579 (Generalized),
                 ...360463997 (Spectral) -- the three independent implementations of the
                 SAME tangent disagree among themselves by 1.3e-10 relative, exactly the
                 size of the FD discrepancy. No implementation could pass a 1e-10 gate.
              3. The relative deviation is 5.6e-12 ... 1.9e-10 over all fifteen
                 registered combinations, i.e. 10 correct significant digits everywhere.

            THE GATE IS THEREFORE RESCALED, NOT LOOSENED IN SUBSTANCE -- but stated
            plainly: in RELATIVE terms this IS looser than the ~1.2e-10 the absolute
            1e-3 happened to impose on the paths that used to pass. Nobody chose 1.2e-10;
            it was an accident of the entry magnitudes, and it sits below the accuracy
            floor of the FD oracle itself. The new gate scales with real_t precision
            (float/double/multiprecision) and with the magnitude actually compared:
            1e8*eps ~ 2.2e-8 relative in double, which is ~120x the worst floor measured
            above and still demands EIGHT correct significant digits of the tangent.
        */
        // The scale is max(|dC|,|C|) and not |dC| alone: for SvK dC vanishes identically,
        // which would make the gate exactly zero and leave the check hostage to the last
        // bit of an FD difference that only happens to cancel exactly today.
        const real_t dC_scale = math::max(dC_MM.array().abs().maxCoeff(),
                                          C_MM.array().abs().maxCoeff());
        const real_t dC_tol   = 1e8*std::numeric_limits<real_t>::epsilon()*dC_scale;
        gsInfo << "[MM_DIAG] mat "<<material<<" impl "<<impl
               << (Compressibility ? " comp" : " incomp")
               << " : |C|max "<<C_MM.array().abs().maxCoeff()
               << " dC_absdiff "<<(dC_FD-dC_MM).array().abs().maxCoeff()
               << " dC_reldiff "<<(dC_FD-dC_MM).array().abs().maxCoeff()/dC_scale
               << " (gate "<<dC_tol<<") ; C_reldiff "
               << (C_FD-C_MM).array().abs().maxCoeff()/C_MM.array().abs().maxCoeff()
               << "\n";
        CHECK_MATRIX_CLOSE(dC_FD,dC_MM,dC_tol);

        /// STRESS TFT
        STFTfun<real_t> STFT(gori,materialMatrix.get(),0,pt,z);
        STFT.eval_into(e,resss);
        gsMatrix<> STFT_test = resss;
        gsMatrix<> STFT_MM   = vecTFT.piece(0).eval(pt).reshape(3,1);
        /*  (task 65) SECOND INSTANCE OF THE VERY DEFECT THIS TASK WAS SENT TO FIX.
            This line read
                CHECK_MATRIX_CLOSE(STFT_test,STFT_test,1e-3);
            i.e. it compared the finite-difference-free reference to ITSELF. STFT_MM
            was computed and dropped, so the TFT stress path (gsMaterialMatrixTFT's
            VectorN output) was gated by NOTHING in any of the fifteen MM_* tests --
            fictional coverage of the same shape as the fictional material arguments,
            and in the same file.
        */
        gsInfo << "[MM_STFT] mat "<<material<<" impl "<<impl
               << (Compressibility ? " comp" : " incomp")
               << " : ref ["<<STFT_test.transpose()<<" ] vs TFT ["<<STFT_MM.transpose()
               << " ] , max|diff| = "<<(STFT_test-STFT_MM).array().abs().maxCoeff()<<"\n";
        CHECK_MATRIX_CLOSE(STFT_test,STFT_MM,1e-3);

        /// MATRIX TFT
        STFT.deriv_into(e,resss);
        gsMatrix<> CTFT_FD = resss.reshape(3,3);
        gsMatrix<> CTFT_MM = matTFT.piece(0).eval(pt).reshape(3,3);
        CHECK_MATRIX_CLOSE(CTFT_FD,CTFT_MM,1e-3);
    }


}
