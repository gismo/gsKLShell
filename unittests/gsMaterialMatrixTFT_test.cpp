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
            gsVector<T> Eeng = u.col(k); // engineering Voigt strain
            gsVector<T> E = Eeng;
            E(2,0) *= 0.5;
            gsMatrix<T> Cmat = 2*E+m_gori;

            gsMatrix<T> C = m_mm->eval3D_matrix_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);
            gsMatrix<T> S = m_mm->eval3D_vector_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);

            gsMatrix<T> THETA = m_mmTFT->eval_theta(C,S,Eeng);

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

/// Tolerance relative to the largest entry of \a ref
static inline real_t TFT_relTol(const gsMatrix<real_t> & ref, const real_t relTol)
{
    return relTol * ref.cwiseAbs().maxCoeff();
}

/// Truncated (tension-field) stress S + C*gamma*n1 on the geometry-free "_C" path, as a function
/// of the engineering Voigt strain. Records, for every evaluated column, whether the state is
/// wrinkled (S has one negative and one positive principal value, E has a positive principal
/// value, gamma > 0, and n2^T (S + gamma C n1) = 0 holds at the root theta).
template<class T>
class TFTStressStencilFun : public gsFunction<T>
{
public:
    TFTStressStencilFun(const gsMatrix<T> & gori, gsMaterialMatrixBase<T> * mm, const gsMaterialMatrixTFT<2,T,false> & mmTFT, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_gori(gori), m_mm(mm), m_mmTFT(&mmTFT), m_patch(patch), m_u(u), m_z(z),
    m_columns(0), m_allWrinkled(true)
    {
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        for (index_t k=0; k!=u.cols(); k++)
        {
            gsVector<T> Eeng = u.col(k);
            gsVector<T> E = Eeng;
            E(2,0) *= 0.5;
            gsMatrix<T> Cmat = 2*E+m_gori;

            gsMatrix<T> C = m_mm->eval3D_matrix_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);
            gsMatrix<T> S = m_mm->eval3D_vector_C(Cmat,m_patch,m_u.col(0),m_z(0,0),MaterialOutput::Generic);

            gsMatrix<T> THETA = m_mmTFT->eval_theta(C,S,Eeng);

            T n1 = math::cos(THETA(0,0));
            T n2 = math::sin(THETA(0,0));
            T m1 = -math::sin(THETA(0,0));
            T m2 = math::cos(THETA(0,0));
            gsVector<T> n1_vec(3); n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            gsVector<T> n2_vec(3); n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            gsMatrix<T> C33 = C.reshape(3,3);
            T gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C33 * n1_vec ).value();
            gsVector<T> Sp = S + C33 * (gamma * n1_vec);
            result.col(k) = Sp;

            // Principal values of the 2x2 symmetric component matrices (signs are metric independent)
            T f = ( n2_vec.transpose() * Sp ).value();
            T trS = S(0)+S(1), detS = S(0)*S(1)-S(2)*S(2);
            T minS = trS/2 - math::sqrt(trS*trS/4 - detS);
            T maxS = trS/2 + math::sqrt(trS*trS/4 - detS);
            T trE = Eeng(0)+Eeng(1), detE = Eeng(0)*Eeng(1)-Eeng(2)*Eeng(2)/4;
            T maxE = trE/2 + math::sqrt(trE*trE/4 - detE);

            const bool wrinkled = minS < 0 && maxS > 0 && maxE > 0 && gamma > 0 && math::abs(f) <= 1e-8 * S.norm();
            m_allWrinkled = m_allWrinkled && wrinkled;
            m_columns++;
        }
    }

    short_t domainDim() const {return 3;}
    short_t targetDim() const {return 3;}

    index_t columns() const { return m_columns; }
    bool allWrinkled() const { return m_allWrinkled; }

private:
    const gsMatrix<T> m_gori;
    gsMaterialMatrixBase<T> * m_mm;
    const gsMaterialMatrixTFT<2,T,false> * m_mmTFT;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
    mutable index_t m_columns;
    mutable bool m_allWrinkled;
};

/// Material matrix C (column-major 3x3 in 9 entries) of the wrapped material as a function of the
/// engineering Voigt strain.
template<class T>
class TFTCmatStencilFun : public gsFunction<T>
{
public:
    TFTCmatStencilFun(const gsMatrix<T> & gori, gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    : m_gori(gori), m_mm(mm), m_patch(patch), m_u(u), m_z(z)
    {
    }

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
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

SUITE(gsMaterialMatrixTFT_test)                 // The suite should have the same name as the file
{
    void MM_CHECK(const index_t material, const index_t impl, const bool Compressibility);
    void TFT_tangent_CHECK(const bool svk);

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

    void MM_CHECK(const index_t material, const index_t impl, const bool Compressibility)
    {
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
        mp_def.patch(0).coefs().col(1) *= 0.9; // over-contract the deformed geometry so that the point is wrinkled
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

        // Tolerances are relative to the largest entry of the reference quantity
        const bool compressible = Compressibility || material==0;
        const real_t relTol = compressible ? 1e-5 : 1e-6;

        /// STRAIN
        gsMatrix<> e = 1./2. * (Cmat - gori);
        e(2,0) *= 2;
        gsMatrix<> strain_MM = materialMatrix->eval3D_strain(0,pt,z);
        CHECK_MATRIX_CLOSE(e,strain_MM,TFT_relTol(e,relTol));

        /// WRINKLED STATE
        CHECK_EQUAL(0,(index_t)materialMatrixTFT->eval3D_tensionfield(0,pt,z,MaterialOutput::Generic)(0,0));

        /// MATRIX
        Sfun<real_t> Sfun(gori,materialMatrix.get(),0,pt,z);
        Sfun.deriv_into(e,resss);
        gsMatrix<> C_FD = resss.reshape(3,3); // Finite Differences
        gsMatrix<> C_MM = materialMatrix->eval3D_matrix(0,pt,z,MaterialOutput::Generic).reshape(3,3); // implemented
        CHECK_MATRIX_CLOSE(C_FD,C_MM,TFT_relTol(C_FD,relTol));

        /// dMATRIX
        Cfun<real_t> Cfun(gori,materialMatrix.get(),0,pt,z);
        Cfun.deriv_into(e,resss);
        gsMatrix<> dC_FD = resss.reshape(9,3);
        gsMatrix<> dC_MM = materialMatrix->eval3D_dmatrix(0,pt,z,MaterialOutput::Generic).reshape(9,3);
        dC_MM *= 2; // NOTE: dmatrix returns d(mm)/dC, and dC_FD is d(mm)/dE = 2*d(mm)/dC
        CHECK_MATRIX_CLOSE(dC_FD,dC_MM,TFT_relTol(dC_FD,relTol));

        /// STRESS TFT
        STFTfun<real_t> STFT(gori,materialMatrix.get(),0,pt,z);
        STFT.eval_into(e,resss);
        gsMatrix<> STFT_test = resss;
        gsMatrix<> STFT_MM   = vecTFT.piece(0).eval(pt).reshape(3,1);
        CHECK_MATRIX_CLOSE(STFT_test,STFT_MM,TFT_relTol(STFT_test,relTol));

        /// MATRIX TFT
        STFT.deriv_into(e,resss);
        // deriv_into returns d(S_c)/d(E_j) at index j + 3c, so reshaping column-wise gives the transposed tangent
        gsMatrix<> CTFT_FD = resss.reshape(3,3).transpose();
        gsMatrix<> CTFT_MM = matTFT.piece(0).eval(pt).reshape(3,3);
        CHECK_MATRIX_CLOSE(CTFT_FD,CTFT_MM,TFT_relTol(CTFT_MM,relTol));
    }


    /// Truncated-stress tangent on an affine patch pair with a prescribed Jacobian,
    /// compared with finite differences of the truncated stress.
    /// svk == false: incompressible Neo-Hookean (nonlinear _compute_C, rotated frame, shear in n1);
    /// svk == true : SvK with a skewed reference metric (linear _compute_C).
    void TFT_tangent_CHECK(const bool svk)
    {
        const real_t pi = EIGEN_PI;
        const real_t al = svk ? 0.9 : pi/6;
        gsMatrix<> R(2,2);
        R << math::cos(al), -math::sin(al), math::sin(al), math::cos(al);

        gsMatrix<> Jref(2,2), Jdef(2,2);
        if (svk)
        {
            Jref << 1.0, 0.3, 0.0, 0.8;
            gsMatrix<> D(2,2); D << 1.1, 0.0, 0.0, 0.75;
            gsMatrix<> F = R * D * R.transpose();
            Jdef = F * Jref;
        }
        else
        {
            const real_t lam1 = 1.2;
            const real_t lam2 = 0.8 * math::pow(lam1,-0.5);
            Jref.setIdentity();
            gsMatrix<> D(2,2); D << lam1, 0.0, 0.0, lam2;
            Jdef = R * D * R.transpose();
        }

        gsMultiPatch<> mp, mp_def;
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        mp.patch(0).coefs() = mp.patch(0).coefs() * Jref.transpose();
        mp.addAutoBoundaries();
        mp_def.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        mp_def.patch(0).coefs() = mp_def.patch(0).coefs() * Jdef.transpose();
        mp_def.addAutoBoundaries();

        const real_t mu = 4.225e5;
        const real_t PoissonRatio = svk ? 0.45 : 0.5;
        const real_t E_modulus = 2*mu*(1+PoissonRatio);
        gsFunctionExpr<> t(std::to_string(0.01),2);
        gsFunctionExpr<> E(std::to_string(E_modulus),2);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
        std::vector<gsFunctionSet<>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;

        gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
        gsMaterialMatrixBase<real_t>::uPtr materialMatrixTFT;
        gsOptionList options;
        if (svk)
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
            materialMatrixTFT = memory::make_unique(new gsMaterialMatrixTFT<2,real_t,true>(*materialMatrix));
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
            materialMatrixTFT = memory::make_unique(new gsMaterialMatrixTFT<2,real_t,false>(*materialMatrix));
        }
        materialMatrix->setDeformed(&mp_def);
        materialMatrixTFT->setDeformed(&mp_def);

        gsVector<> pt(2);
        pt.setConstant(0.25);
        gsMatrix<> z(1,1); z.setZero();

        // Reference metric (G11,G22,G12) and strain (engineering Voigt) at the point
        gsMatrix<> Gref = Jref.transpose() * Jref;
        gsMatrix<> Gdef = Jdef.transpose() * Jdef;
        gsMatrix<> gori(3,1);
        gori << Gref(0,0), Gref(1,1), Gref(0,1);
        gsMatrix<> Ebase = materialMatrix->eval3D_strain(0,pt,z);
        gsMatrix<> Eexp(3,1);
        Eexp << 0.5*(Gdef(0,0)-Gref(0,0)), 0.5*(Gdef(1,1)-Gref(1,1)), Gdef(0,1)-Gref(0,1);
        CHECK_MATRIX_CLOSE(Eexp,Ebase,1e-10);

        // 1. Production state is wrinkled
        CHECK_EQUAL(0,(index_t)materialMatrixTFT->eval3D_tensionfield(0,pt,z,MaterialOutput::Generic)(0,0));

        // 2. Root theta at the base point
        gsMaterialMatrixTFT<2,real_t,false> tftHelper(materialMatrix.get());
        gsMatrix<> C0 = materialMatrix->eval3D_matrix(0,pt,z,MaterialOutput::Generic);
        gsMatrix<> S0 = materialMatrix->eval3D_vector(0,pt,z,MaterialOutput::Generic);
        const real_t theta = tftHelper.eval_theta(C0,S0,Ebase)(0,0);
        if (!svk)
        {
            real_t d = std::fmod(theta - 2*pi/3, pi);
            if (d >  pi/2) d -= pi;
            if (d < -pi/2) d += pi;
            CHECK_CLOSE(0.0,d,1e-6);
        }

        // 3. gamma > 0 at the base point
        gsMatrix<> C33 = C0.reshape(3,3);
        gsVector<> n1_vec(3), n2_vec(3);
        n1_vec << math::pow(math::cos(theta),2), math::pow(math::sin(theta),2), 2*math::cos(theta)*math::sin(theta);
        n2_vec << -math::sin(theta)*math::cos(theta), math::sin(theta)*math::cos(theta), math::pow(math::cos(theta),2)-math::pow(math::sin(theta),2);
        const real_t n1Cn1 = (n1_vec.transpose() * C33 * n1_vec).value();
        const real_t gamma = - (n1_vec.transpose() * S0).value() / n1Cn1;
        CHECK(gamma > 0);

        // 4./5. FD tangent of the truncated stress, every stencil point wrinkled
        TFTStressStencilFun<real_t> stencil(gori,materialMatrix.get(),tftHelper,0,pt,z);
        gsMatrix<> res;
        stencil.deriv_into(Ebase,res);
        CHECK_EQUAL(12,stencil.columns());
        CHECK(stencil.allWrinkled());
        gsMatrix<> J_FD = res.reshape(3,3).transpose();
        gsMatrix<> J_MM = materialMatrixTFT->eval3D_matrix(0,pt,z,MaterialOutput::MatrixA).reshape(3,3);
        const real_t tol = 1e-6 * J_MM.cwiseAbs().maxCoeff();
        CHECK_MATRIX_CLOSE(J_FD,J_MM,tol);

        if (svk)
        {
            // The sign of the n2^T C n1 term in df/dE is only visible when n2^T C n1 != 0
            const real_t n2Cn1 = (n2_vec.transpose() * C33 * n1_vec).value();
            CHECK(math::abs(n2Cn1) / n1Cn1 >= 0.05);
        }
        else
        {
            // 6. The tangent is not symmetric, so the orientation matters
            CHECK(((J_FD - J_FD.transpose()).cwiseAbs().maxCoeff()) >= 100 * tol);

            // 7. gamma * H, H_id = sum_c dC_ic/dE_d n1_c, contributes to the tangent
            TFTCmatStencilFun<real_t> Cstencil(gori,materialMatrix.get(),0,pt,z);
            gsMatrix<> resC;
            Cstencil.deriv_into(Ebase,resC);
            gsMatrix<> H_FD(3,3);
            H_FD.setZero();
            for (index_t i=0; i!=3; i++)
                for (index_t d=0; d!=3; d++)
                    for (index_t c=0; c!=3; c++)
                        H_FD(i,d) += resC.data()[d + 3*(i + 3*c)] * n1_vec(c);
            CHECK(gamma * H_FD.cwiseAbs().maxCoeff() >= 100 * tol);
        }
    }

    TEST(TFT_tangent_NH_incomp_orientation)
    {
        TFT_tangent_CHECK(false);
    }

    TEST(TFT_tangent_SvK_linear_skewed)
    {
        TFT_tangent_CHECK(true);
    }

    /// Neo-Hookean compressible material on the undeformed geometry \a mp only; the material keeps
    /// non-owning references to the (static) thickness and parameter functions
    static gsMaterialMatrixBase<real_t>::uPtr HandOverUndeformed_material(const gsMultiPatch<> & mp)
    {
        const real_t mu = 4.225e5, nuv = 0.45;
        static gsFunctionExpr<> t("0.01",2);
        static gsFunctionExpr<> E(std::to_string(2*mu*(1+nuv)),2);
        static gsFunctionExpr<> nu(std::to_string(nuv),2);
        std::vector<gsFunctionSet<>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;

        gsOptionList options;
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",true);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        return getMaterialMatrix<2,real_t>(mp,t,parameters,options);
    }

    /// Wraps a material that has only the undeformed geometry in \a Decorator and checks that the
    /// decorator receives the undeformed geometry at construction and evaluates identically to a
    /// decorator whose geometries were set explicitly. A fresh material is built per call because
    /// the decorators forward setDeformed to the wrapped material.
    template<class Decorator>
    static void HandOverUndeformed_CHECK(const gsMultiPatch<> & mp, const gsMultiPatch<> & mp_def)
    {
        typedef gsMaterialMatrixBase<real_t>::Ptr material_ptr;

        gsMaterialMatrixBase<real_t>::uPtr material = HandOverUndeformed_material(mp);
        CHECK(material->hasUndeformed());
        CHECK(!material->hasDeformed());

        Decorator dec(material_ptr(memory::make_shared_not_owned(material.get())));
        CHECK(dec.getUndeformed() != nullptr);
        CHECK(!dec.hasDeformed());
        if (dec.getUndeformed() == nullptr)
            return;
        CHECK(dec.getUndeformed() == material->getUndeformed());

        gsMaterialMatrixBase<real_t>::uPtr materialRef = HandOverUndeformed_material(mp);
        Decorator ref(material_ptr(memory::make_shared_not_owned(materialRef.get())));
        ref.setUndeformed(materialRef->getUndeformed());
        ref.setDeformed(&mp_def);
        dec.setDeformed(&mp_def);

        gsVector<> pt(2);
        pt.setConstant(0.25);
        gsMatrix<> z(1,1);
        z.setZero();

        const gsMatrix<> Cdec = dec.eval3D_matrix(0,pt,z,MaterialOutput::MatrixA);
        const gsMatrix<> Cref = ref.eval3D_matrix(0,pt,z,MaterialOutput::MatrixA);
        const gsMatrix<> Sdec = dec.eval3D_vector(0,pt,z,MaterialOutput::VectorN);
        const gsMatrix<> Sref = ref.eval3D_vector(0,pt,z,MaterialOutput::VectorN);
        CHECK(Cref.cwiseAbs().maxCoeff() > 0);
        CHECK(Sref.cwiseAbs().maxCoeff() > 0);
        CHECK_MATRIX_CLOSE(Cref,Cdec,1e-12*Cref.cwiseAbs().maxCoeff());
        CHECK_MATRIX_CLOSE(Sref,Sdec,1e-12*Sref.cwiseAbs().maxCoeff());
    }

    TEST(TFT_wrapped_material_hands_over_undeformed)
    {
        gsMultiPatch<> mp_def;
        mp_def.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        mp_def.addAutoBoundaries();
        gsMultiPatch<> mp = mp_def;
        mp.patch(0).coefs()(2,0) += 0.1;
        mp_def.patch(0).coefs().col(1) *= 0.9;

        gsMaterialMatrixLinear<2,real_t> bare;
        CHECK(!bare.hasUndeformed());
        bare.setDeformed(&mp_def);
        CHECK(bare.hasDeformed());
        CHECK(!bare.hasUndeformed());

        HandOverUndeformed_CHECK< gsMaterialMatrixTFT<2,real_t,false> >(mp,mp_def);
        HandOverUndeformed_CHECK< gsMaterialMatrixTFTSplit<2,real_t,false> >(mp,mp_def);
    }

    /// Largest deviation of the columns of the through-thickness blocks 1..nz-1 from block 0.
    /// The columns are ordered [(u1,z1) .. (un,z1), (u1,z2) .. (un,z2), ..], n = number of in-plane points.
    static real_t TFT_zVariation(const gsMatrix<real_t> & M, const index_t n, const index_t nz)
    {
        real_t d = 0;
        for (index_t j=1; j<nz; ++j)
            d = math::max(d,(M.middleCols(j*n,n) - M.leftCols(n)).cwiseAbs().maxCoeff());
        return d;
    }

    /// Linear (SvK) material on a flat reference square deformed into a curved surface (stretches plus a quadratic
    /// out-of-plane bump of the deformed configuration) wrapped in the
    /// TFT decorator (linear = true) and in the Split decorator with Eta = 0. With a linear base law the
    /// wrinkle state is the mid-surface one at every through-thickness point, so the TFT output must
    /// not depend on z and must coincide with the Split decorator at Eta = 0.
    /// \a stretch1, \a stretch2 are the in-plane stretches of the deformed surface; \a expectedState
    /// the wrinkling state (1 taut, 0 wrinkled, -1 slack) they are chosen to produce.
    static void TFT_linear_curved_CHECK(const real_t stretch1, const real_t stretch2, const index_t expectedState)
    {
        typedef real_t T;
        const T bump = 0.15;

        gsMultiPatch<T> ref, def;
        ref.addPatch(gsNurbsCreator<T>::BSplineSquare(1));
        ref.patch(0).degreeElevate(1);
        ref.embed(3);
        def = ref;
        gsMatrix<T> & c = def.patch(0).coefs();
        for (index_t i=0; i!=c.rows(); ++i)
        {
            const T x = c(i,0), y = c(i,1);
            c(i,0) = stretch1 * x;
            c(i,1) = stretch2 * y;
            c(i,2) = bump * (x*(1-x) + y*(1-y) + 2*x*y);
        }
        ref.addAutoBoundaries();
        def.addAutoBoundaries();

        gsConstantFunction<T> t(0.01,3), E(1e5,3), nu(0.3,3);
        std::vector<gsFunctionSet<T>*> par(2);
        par[0] = &E;
        par[1] = &nu;
        gsOptionList o;
        o.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        o.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        gsMaterialMatrixBase<T>::uPtr wrapped = getMaterialMatrix<3,T>(ref,t,par,o);

        gsMaterialMatrixTFT<3,T,true> tft(wrapped.get());
        gsMaterialMatrixTFTSplit<3,T,true> split(wrapped.get());
        split.options().setReal("Eta",0);
        gsMaterialMatrixBase<T> * all[3] = { &tft, &split, wrapped.get() };
        for (index_t i=0; i!=3; ++i)
        {
            all[i]->setUndeformed(&ref);
            all[i]->setDeformed(&def);
        }

        // Three in-plane points; per point a different, non-uniform set of four through-thickness
        // heights spanning the thickness (z is the through-thickness coordinate in [-1/2,1/2], scaled by the thickness inside the material)
        gsMatrix<T> u(2,3);
        u << 0.3, 0.7, 0.5,
             0.4, 0.6, 0.2;
        gsMatrix<T> z(4,3);
        z << -0.5,  -0.4,  -0.3,
              0.0,   0.1,  -0.1,
              0.25,  0.3,   0.2,
              0.5,   0.45,  0.35;
        const index_t n = u.cols(), nz = z.rows();

        const T eps  = std::numeric_limits<T>::epsilon();
        const T rel  = math::max(T(1e-12), T(100)*eps);
        const T tiny = std::numeric_limits<T>::min();

        // The case produces the requested state at all points
        const gsMatrix<T> TF = tft.eval3D_tensionfield(0,u,z);
        CHECK_EQUAL(n*nz,TF.cols());
        for (index_t i=0; i!=TF.cols(); ++i)
            CHECK_EQUAL(expectedState,(index_t)TF(0,i));

        const gsMatrix<T> Et = tft.eval3D_strain(0,u,z);
        const gsMatrix<T> St = tft.eval3D_stress(0,u,z,MaterialOutput::VectorN);
        const gsMatrix<T> Pt = tft.eval3D_pstrain(0,u,z);
        const gsMatrix<T> Es = split.eval3D_strain(0,u,z);
        const gsMatrix<T> Ss = split.eval3D_stress(0,u,z,MaterialOutput::VectorN);

        // Control: the wrapped strain does depend on z on this surface, so the z-independence below is
        // not vacuous
        const gsMatrix<T> Ew = wrapped->eval3D_strain(0,u,z);
        CHECK(TFT_zVariation(Ew,n,nz) > 1e-4 * Ew.cwiseAbs().maxCoeff());

        CHECK_EQUAL(Es.rows(),Et.rows());
        CHECK_EQUAL(Es.cols(),Et.cols());
        CHECK_EQUAL(Ss.cols(),St.cols());

        // z-independence, exact up to rounding of the replicated mid-surface values
        CHECK_EQUAL(T(0),TFT_zVariation(Et,n,nz));
        CHECK_EQUAL(T(0),TFT_zVariation(St,n,nz));
        CHECK_EQUAL(T(0),TFT_zVariation(Pt,n,nz));

        // Equality with the Split decorator at Eta = 0 (strain: at wrinkled points the elastic strain
        // E + gamma n1 of Split is the TFT strain output)
        CHECK_MATRIX_CLOSE(Es,Et,rel*Es.cwiseAbs().maxCoeff() + tiny);
        CHECK_MATRIX_CLOSE(Ss,St,rel*Ss.cwiseAbs().maxCoeff() + tiny);

        // Principal strain, one in-plane point per call: the multi-point call evaluates all points with
        // the metric of the last one. At wrinkled points TFT returns the total and Split the elastic
        // principal strain, so the comparison is restricted to taut and slack points.
        if (expectedState != 0)
        {
            for (index_t k=0; k!=n; ++k)
            {
                const gsMatrix<T> uk = u.col(k), zk = z.col(k);
                const gsMatrix<T> P1 = tft.eval3D_pstrain(0,uk,zk);
                const gsMatrix<T> P2 = split.eval3D_pstrain(0,uk,zk);
                CHECK_MATRIX_CLOSE(P2,P1,rel*P2.cwiseAbs().maxCoeff() + tiny);
            }
        }

        // Single in-plane point and single through-thickness row give the same values as the block call
        for (index_t k=0; k!=n; ++k)
        {
            const gsMatrix<T> uk = u.col(k), zk = z.col(k);
            const gsMatrix<T> E1 = tft.eval3D_strain(0,uk,zk);
            const gsMatrix<T> S1 = tft.eval3D_stress(0,uk,zk,MaterialOutput::VectorN);
            for (index_t j=0; j!=nz; ++j)
            {
                CHECK_MATRIX_CLOSE(Et.col(j*n+k),E1.col(j),rel*Et.cwiseAbs().maxCoeff() + tiny);
                CHECK_MATRIX_CLOSE(St.col(j*n+k),S1.col(j),rel*St.cwiseAbs().maxCoeff() + tiny);
            }
        }
        const gsMatrix<T> z1 = gsMatrix<T>::Constant(1,n,0.37);
        CHECK_MATRIX_CLOSE(split.eval3D_strain(0,u,z1),tft.eval3D_strain(0,u,z1),rel*Es.cwiseAbs().maxCoeff() + tiny);
    }

    TEST(TFT_linear_curved_reference_wrinkled)
    {
        TFT_linear_curved_CHECK(1.2,0.8/math::sqrt(1.2),0);
    }

    TEST(TFT_linear_curved_reference_taut)
    {
        TFT_linear_curved_CHECK(1.2,0.97,1);
    }

    TEST(TFT_linear_curved_reference_slack)
    {
        TFT_linear_curved_CHECK(0.95,0.90,-1);
    }

    /// Curved d = 3 SvK state of TFT_linear_curved_CHECK: the wrapped law, a TFT decorator (linear = true) and a
    /// Split decorator around it, and the evaluation points (three in-plane points, four heights each).
    struct TFT_curved_setup
    {
        typedef real_t T;

        gsMultiPatch<T> ref, def;
        gsConstantFunction<T> t, E, nu;
        gsMaterialMatrixBase<T>::uPtr wrapped;
        memory::unique_ptr< gsMaterialMatrixTFT<3,T,true> >      tft;
        memory::unique_ptr< gsMaterialMatrixTFTSplit<3,T,true> > split;
        gsMatrix<T> u, z;

        TFT_curved_setup(const T stretch1, const T stretch2)
        :
        t(0.01,3), E(1e5,3), nu(0.3,3)
        {
            const T bump = 0.15;
            ref.addPatch(gsNurbsCreator<T>::BSplineSquare(1));
            ref.patch(0).degreeElevate(1);
            ref.embed(3);
            def = ref;
            gsMatrix<T> & c = def.patch(0).coefs();
            for (index_t i=0; i!=c.rows(); ++i)
            {
                const T x = c(i,0), y = c(i,1);
                c(i,0) = stretch1 * x;
                c(i,1) = stretch2 * y;
                c(i,2) = bump * (x*(1-x) + y*(1-y) + 2*x*y);
            }
            ref.addAutoBoundaries();
            def.addAutoBoundaries();

            std::vector<gsFunctionSet<T>*> par(2);
            par[0] = &E;
            par[1] = &nu;
            gsOptionList o;
            o.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            o.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            wrapped = getMaterialMatrix<3,T>(ref,t,par,o);

            tft   = memory::make_unique(new gsMaterialMatrixTFT<3,T,true>(wrapped.get()));
            split = memory::make_unique(new gsMaterialMatrixTFTSplit<3,T,true>(wrapped.get()));
            gsMaterialMatrixBase<T> * all[3] = { tft.get(), split.get(), wrapped.get() };
            for (index_t i=0; i!=3; ++i)
            {
                all[i]->setUndeformed(&ref);
                all[i]->setDeformed(&def);
            }

            u.resize(2,3);
            u << 0.3, 0.7, 0.5,
                 0.4, 0.6, 0.2;
            z.resize(4,3);
            z << -0.5,  -0.4,  -0.3,
                  0.0,   0.1,  -0.1,
                  0.25,  0.3,   0.2,
                  0.5,   0.45,  0.35;
        }
    };

    /// \f$\max|b-a| \le \mathrm{rel}\max|a| + \mathrm{tiny}\f$
    static void TFT_eta_CLOSE(const gsMatrix<real_t> & a, const gsMatrix<real_t> & b, const real_t rel)
    {
        CHECK_EQUAL(a.rows(),b.rows());
        CHECK_EQUAL(a.cols(),b.cols());
        if (a.rows()!=b.rows() || a.cols()!=b.cols())
            return;
        CHECK_MATRIX_CLOSE(a,b,rel*a.cwiseAbs().maxCoeff() + std::numeric_limits<real_t>::min());
    }

    /// Residual stiffness Eta of the TFT decorator against the Split decorator (same definition, Eqs. 35 and 40 of
    /// Zhang & Kiendl 2026): stress, tangent and strain at the wrinkled, taut and slack curved states.
    static void TFT_eta_split_CHECK(const real_t stretch1, const real_t stretch2, const index_t expectedState)
    {
        typedef real_t T;
        TFT_curved_setup s(stretch1,stretch2);
        const T rel = 1e-12;
        const MaterialOutput MA = MaterialOutput::MatrixA;
        const MaterialOutput VN = MaterialOutput::VectorN;

        const gsMatrix<T> TF = s.tft->eval3D_tensionfield(0,s.u,s.z);
        for (index_t i=0; i!=TF.cols(); ++i)
            CHECK_EQUAL(expectedState,(index_t)TF(0,i));

        s.tft->options().setReal("SlackMultiplier",0);
        s.split->options().setReal("SlackMultiplier",0);

        // Eta = 0: the tangents agree to the same tolerance (the root-searched theta of TFT enters to first order)
        const gsMatrix<T> Cs0 = s.split->eval3D_matrix(0,s.u,s.z,MA);
        const gsMatrix<T> Ct0 = s.tft->eval3D_matrix(0,s.u,s.z,MA);
        gsInfo<<"state "<<expectedState<<", Eta 0 tangent gap TFT vs Split (relative): "
              <<(Cs0-Ct0).cwiseAbs().maxCoeff()/(Cs0.cwiseAbs().maxCoeff() + std::numeric_limits<T>::min())<<"\n";
        TFT_eta_CLOSE(Cs0,Ct0,rel);

        const T eta = 1e-3;
        s.tft->options().setReal("Eta",eta);
        s.split->options().setReal("Eta",eta);
        const gsMatrix<T> Ct = s.tft->eval3D_matrix(0,s.u,s.z,MA);
        const gsMatrix<T> Cs = s.split->eval3D_matrix(0,s.u,s.z,MA);
        const gsMatrix<T> St = s.tft->eval3D_stress(0,s.u,s.z,VN);
        const gsMatrix<T> Ss = s.split->eval3D_stress(0,s.u,s.z,VN);
        const gsMatrix<T> Et = s.tft->eval3D_strain(0,s.u,s.z);
        const gsMatrix<T> Es = s.split->eval3D_strain(0,s.u,s.z);
        gsInfo<<"state "<<expectedState<<", Eta 1e-3 gaps TFT vs Split (relative): tangent "
              <<(Cs-Ct).cwiseAbs().maxCoeff()/(Cs.cwiseAbs().maxCoeff() + std::numeric_limits<T>::min())<<", stress "
              <<(Ss-St).cwiseAbs().maxCoeff()/(Ss.cwiseAbs().maxCoeff() + std::numeric_limits<T>::min())<<", strain "
              <<(Es-Et).cwiseAbs().maxCoeff()/(Es.cwiseAbs().maxCoeff() + std::numeric_limits<T>::min())<<"\n";
        TFT_eta_CLOSE(Cs,Ct,rel);
        TFT_eta_CLOSE(Ss,St,rel);
        TFT_eta_CLOSE(Es,Et,rel);

        // Control: Eta changes the tangent and, away from taut points, the stress
        s.tft->options().setReal("Eta",0);
        const gsMatrix<T> C0 = s.tft->eval3D_matrix(0,s.u,s.z,MA);
        const gsMatrix<T> S0 = s.tft->eval3D_stress(0,s.u,s.z,VN);
        s.tft->options().setReal("Eta",eta);
        const gsMatrix<T> S1 = s.tft->eval3D_stress(0,s.u,s.z,VN);
        if (expectedState != 1)
        {
            CHECK((Ct-C0).cwiseAbs().maxCoeff() > 1e-6*Ct.cwiseAbs().maxCoeff());
            CHECK((S1-S0).cwiseAbs().maxCoeff() > 1e-6*S0.cwiseAbs().maxCoeff());
        }
        else
        {
            CHECK_EQUAL(T(0),(Ct-C0).cwiseAbs().maxCoeff());
            CHECK_EQUAL(T(0),(S1-S0).cwiseAbs().maxCoeff());
        }
    }

    TEST(TFT_eta_vs_split_wrinkled)
    {
        TFT_eta_split_CHECK(1.2,0.8/math::sqrt(1.2),0);
    }

    TEST(TFT_eta_vs_split_taut)
    {
        TFT_eta_split_CHECK(1.2,0.97,1);
    }

    TEST(TFT_eta_vs_split_slack)
    {
        TFT_eta_split_CHECK(0.95,0.90,-1);
    }

    /// Tangent-only residual stiffness: stress and strain are those of Eta = 0 exactly, the tangent is blended.
    static void TFT_eta_tangent_only_CHECK(const real_t stretch1, const real_t stretch2, const index_t expectedState)
    {
        typedef real_t T;
        TFT_curved_setup s(stretch1,stretch2);
        const MaterialOutput MA = MaterialOutput::MatrixA;
        const MaterialOutput VN = MaterialOutput::VectorN;

        const gsMatrix<T> TF = s.tft->eval3D_tensionfield(0,s.u,s.z);
        for (index_t i=0; i!=TF.cols(); ++i)
            CHECK_EQUAL(expectedState,(index_t)TF(0,i));

        const gsMatrix<T> S0 = s.tft->eval3D_stress(0,s.u,s.z,VN);
        const gsMatrix<T> E0 = s.tft->eval3D_strain(0,s.u,s.z);
        const gsMatrix<T> M0 = s.tft->eval3D_matrix(0,s.u,s.z,MA);
        const gsMatrix<T> Cb = s.wrapped->eval3D_matrix(0,s.u,s.z,MA);

        const T eta = 1e-3;
        s.tft->options().setReal("Eta",eta);
        s.tft->options().setSwitch("EtaTangentOnly",true);
        CHECK_EQUAL(T(0),(s.tft->eval3D_stress(0,s.u,s.z,VN)-S0).cwiseAbs().maxCoeff());
        CHECK_EQUAL(T(0),(s.tft->eval3D_strain(0,s.u,s.z)-E0).cwiseAbs().maxCoeff());

        gsMatrix<T> Mexp;
        if (expectedState == 0)
            Mexp = (1-eta)*M0 + eta*Cb;
        else if (expectedState == -1)
            Mexp = eta*Cb;
        else
            Mexp = Cb;
        const gsMatrix<T> M = s.tft->eval3D_matrix(0,s.u,s.z,MA);
        TFT_eta_CLOSE(Mexp,M,1e-14);

        // The blend is not vacuous
        if (expectedState != 1)
            CHECK((M-M0).cwiseAbs().maxCoeff() > 1e-6*M.cwiseAbs().maxCoeff());
    }

    TEST(TFT_eta_tangent_only_wrinkled)
    {
        TFT_eta_tangent_only_CHECK(1.2,0.8/math::sqrt(1.2),0);
    }

    TEST(TFT_eta_tangent_only_taut)
    {
        TFT_eta_tangent_only_CHECK(1.2,0.97,1);
    }

    TEST(TFT_eta_tangent_only_slack)
    {
        TFT_eta_tangent_only_CHECK(0.95,0.90,-1);
    }

    /// \f$\max|a|\f$ with a floor, for relative comparisons
    static real_t TFT_eta_maxAbs(const gsMatrix<real_t> & a)
    {
        return a.cwiseAbs().maxCoeff() + std::numeric_limits<real_t>::min();
    }

    /// Principal strains under Eta > 0, default and tangent-only variants (SlackMultiplier = 0). Taut: untouched. Slack:
    /// scaled by Eta (default) or by SlackMultiplier = 0 (tangent-only, as at Eta = 0); equal to the Split decorator
    /// in the default variant. Wrinkled: neither variant blends the strain, TFT returns the principal strains of the
    /// wrapped law (total strain), whereas Split returns those of the elastic strain, so only TFT-internal values are
    /// compared.
    static void TFT_eta_pstrain_CHECK(const real_t stretch1, const real_t stretch2, const index_t expectedState)
    {
        typedef real_t T;
        TFT_curved_setup s(stretch1,stretch2);
        const T rel = 1e-12;
        const T eta = 1e-3;

        const gsMatrix<T> TF = s.tft->eval3D_tensionfield(0,s.u,s.z);
        for (index_t i=0; i!=TF.cols(); ++i)
            CHECK_EQUAL(expectedState,(index_t)TF(0,i));

        s.tft->options().setReal("SlackMultiplier",0);
        s.split->options().setReal("SlackMultiplier",0);

        const gsMatrix<T> P0 = s.tft->eval3D_pstrain(0,s.u,s.z);
        const gsMatrix<T> Pw = s.wrapped->eval3D_pstrain(0,s.u,gsMatrix<T>::Zero(1,s.u.cols())).replicate(1,s.z.rows());

        s.tft->options().setReal("Eta",eta);
        s.split->options().setReal("Eta",eta);
        const gsMatrix<T> Pd = s.tft->eval3D_pstrain(0,s.u,s.z);
        s.tft->options().setSwitch("EtaTangentOnly",true);
        const gsMatrix<T> Pt = s.tft->eval3D_pstrain(0,s.u,s.z);

        // Tangent-only does not touch the strain output
        CHECK_EQUAL(T(0),(Pt-P0).cwiseAbs().maxCoeff());

        if (expectedState == 1)
        {
            CHECK_EQUAL(T(0),(Pd-P0).cwiseAbs().maxCoeff());
            const gsMatrix<T> Ps = s.split->eval3D_pstrain(0,s.u,s.z);
            TFT_eta_CLOSE(Ps.topRows(2),Pd.topRows(2),rel);
        }
        else if (expectedState == -1)
        {
            // slack factor Eta with SlackMultiplier = 0; zero at Eta = 0
            CHECK_EQUAL(T(0),P0.cwiseAbs().maxCoeff());
            CHECK(Pd.cwiseAbs().maxCoeff() > 1e-6*TFT_eta_maxAbs(Pw));
            const gsMatrix<T> Ps = s.split->eval3D_pstrain(0,s.u,s.z);
            TFT_eta_CLOSE(Ps.topRows(2),Pd.topRows(2),rel);
            TFT_eta_CLOSE(eta*Pw.topRows(2),Pd.topRows(2),1e-14);
        }
        else
        {
            CHECK_EQUAL(T(0),(Pd-P0).cwiseAbs().maxCoeff());
            CHECK_EQUAL(T(0),(Pd-Pw).cwiseAbs().maxCoeff());
            CHECK(P0.cwiseAbs().maxCoeff() > 0);
        }
    }

    TEST(TFT_eta_pstrain_wrinkled)
    {
        TFT_eta_pstrain_CHECK(1.2,0.8/math::sqrt(1.2),0);
    }

    TEST(TFT_eta_pstrain_taut)
    {
        TFT_eta_pstrain_CHECK(1.2,0.97,1);
    }

    TEST(TFT_eta_pstrain_slack)
    {
        TFT_eta_pstrain_CHECK(0.95,0.90,-1);
    }

    /// SlackMultiplier > 0 together with Eta > 0: the slack factor is SlackMultiplier (the rule of the Split
    /// decorator), not Eta. Stress, strain, principal strains and tangent agree with Split; with SlackMultiplier = 0
    /// the factor is Eta and the result differs.
    TEST(TFT_eta_slackmultiplier_overrides_eta)
    {
        typedef real_t T;
        TFT_curved_setup s(0.95,0.90);
        const T rel = 1e-12;
        const MaterialOutput MA = MaterialOutput::MatrixA;
        const MaterialOutput VN = MaterialOutput::VectorN;

        const gsMatrix<T> TF = s.tft->eval3D_tensionfield(0,s.u,s.z);
        for (index_t i=0; i!=TF.cols(); ++i)
            CHECK_EQUAL(-1,(index_t)TF(0,i));

        const T eta = 0.1, slack = 1e-3;
        s.tft->options().setReal("Eta",eta);
        s.split->options().setReal("Eta",eta);

        s.tft->options().setReal("SlackMultiplier",0);
        const gsMatrix<T> S_eta = s.tft->eval3D_stress(0,s.u,s.z,VN);
        const gsMatrix<T> C_eta = s.tft->eval3D_matrix(0,s.u,s.z,MA);
        const gsMatrix<T> E_eta = s.tft->eval3D_strain(0,s.u,s.z);
        const gsMatrix<T> P_eta = s.tft->eval3D_pstrain(0,s.u,s.z);

        s.tft->options().setReal("SlackMultiplier",slack);
        s.split->options().setReal("SlackMultiplier",slack);
        const gsMatrix<T> S = s.tft->eval3D_stress(0,s.u,s.z,VN);
        const gsMatrix<T> C = s.tft->eval3D_matrix(0,s.u,s.z,MA);
        const gsMatrix<T> E = s.tft->eval3D_strain(0,s.u,s.z);
        const gsMatrix<T> P = s.tft->eval3D_pstrain(0,s.u,s.z);

        TFT_eta_CLOSE(s.split->eval3D_stress(0,s.u,s.z,VN),S,rel);
        TFT_eta_CLOSE(s.split->eval3D_matrix(0,s.u,s.z,MA),C,rel);
        TFT_eta_CLOSE(s.split->eval3D_strain(0,s.u,s.z),E,rel);
        TFT_eta_CLOSE(s.split->eval3D_pstrain(0,s.u,s.z).topRows(2),P.topRows(2),rel);

        // Controls: Eta and SlackMultiplier are two orders apart, so the two factors give clearly different results
        CHECK((S-S_eta).cwiseAbs().maxCoeff() > 1e-6*TFT_eta_maxAbs(S_eta));
        CHECK((C-C_eta).cwiseAbs().maxCoeff() > 1e-6*TFT_eta_maxAbs(C_eta));
        CHECK((E-E_eta).cwiseAbs().maxCoeff() > 1e-6*TFT_eta_maxAbs(E_eta));
        CHECK((P-P_eta).cwiseAbs().maxCoeff() > 1e-6*TFT_eta_maxAbs(P_eta));
        TFT_eta_CLOSE(S_eta*(slack/eta),S,1e-12);
    }

    /// A negative SlackMultiplier is rejected at the first evaluation that reads it, with an admissible Eta and
    /// also on the EtaTangentOnly path, where the slack factor is the SlackMultiplier itself.
    TEST(TFT_slackmultiplier_negative_throws)
    {
        typedef real_t T;
        TFT_curved_setup s(0.95,0.90);
        const MaterialOutput MA = MaterialOutput::MatrixA;
        const MaterialOutput VN = MaterialOutput::VectorN;

        const gsMatrix<T> S0 = s.tft->eval3D_stress(0,s.u,s.z,VN);
        CHECK(S0.allFinite());

        s.tft->options().setReal("SlackMultiplier",-1e-3);
        CHECK_THROW(s.tft->eval3D_stress(0,s.u,s.z,VN), std::runtime_error);
        CHECK_THROW(s.tft->eval3D_matrix(0,s.u,s.z,MA), std::runtime_error);
        CHECK_THROW(s.tft->eval3D_strain(0,s.u,s.z), std::runtime_error);
        CHECK_THROW(s.tft->eval3D_energy(0,s.u,s.z), std::runtime_error);

        s.tft->options().setReal("Eta",0.1);
        CHECK_THROW(s.tft->eval3D_stress(0,s.u,s.z,VN), std::runtime_error);

        s.tft->options().setSwitch("EtaTangentOnly",true);
        CHECK_THROW(s.tft->eval3D_stress(0,s.u,s.z,VN), std::runtime_error);
    }

    /// Affine d = 2 setup of TFT_tangent_CHECK (svk = linear): with \a zeroOnly, only the check that Eta = 0 with
    /// EtaTangentOnly set gives exactly the outputs of a default decorator; otherwise the tangent of the default
    /// variant at Eta = 0.1 against finite differences of the decorator's own stress with respect to the deformed
    /// metric.
    template<bool linear>
    static void TFT_eta_affine_CHECK(const bool zeroOnly)
    {
        typedef real_t T;
        const real_t pi = EIGEN_PI;
        const bool svk = linear;
        const real_t al = svk ? 0.9 : pi/6;
        gsMatrix<> R(2,2);
        R << math::cos(al), -math::sin(al), math::sin(al), math::cos(al);

        gsMatrix<> Jref(2,2), Jdef(2,2);
        if (svk)
        {
            Jref << 1.0, 0.3, 0.0, 0.8;
            gsMatrix<> D(2,2); D << 1.1, 0.0, 0.0, 0.75;
            Jdef = (R * D * R.transpose()) * Jref;
        }
        else
        {
            const real_t lam1 = 1.2;
            const real_t lam2 = 0.8 * math::pow(lam1,-0.5);
            Jref.setIdentity();
            gsMatrix<> D(2,2); D << lam1, 0.0, 0.0, lam2;
            Jdef = R * D * R.transpose();
        }

        gsMultiPatch<> mp, mp_def;
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        mp.patch(0).coefs() = mp.patch(0).coefs() * Jref.transpose();
        mp.addAutoBoundaries();
        mp_def.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        mp_def.patch(0).coefs() = mp_def.patch(0).coefs() * Jdef.transpose();
        mp_def.addAutoBoundaries();

        const real_t mu = 4.225e5;
        const real_t PoissonRatio = svk ? 0.45 : 0.5;
        const real_t E_modulus = 2*mu*(1+PoissonRatio);
        gsFunctionExpr<> t(std::to_string(0.01),2);
        gsFunctionExpr<> E(std::to_string(E_modulus),2);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
        std::vector<gsFunctionSet<>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;

        gsOptionList options;
        if (svk)
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        }
        gsMaterialMatrixBase<real_t>::uPtr materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
        materialMatrix->setDeformed(&mp_def);
        gsMaterialMatrixTFT<2,real_t,linear> tft(*materialMatrix);
        tft.setDeformed(&mp_def);

        gsVector<> pt(2);
        pt.setConstant(0.25);
        gsMatrix<> z(1,1); z.setZero();
        const MaterialOutput MA = MaterialOutput::MatrixA;
        const MaterialOutput VN = MaterialOutput::VectorN;
        CHECK_EQUAL(0,(index_t)tft.eval3D_tensionfield(0,pt,z,MaterialOutput::Generic)(0,0));

        if (zeroOnly)
        {
            gsMaterialMatrixTFT<2,real_t,linear> tft0(*materialMatrix);
            tft0.setDeformed(&mp_def);
            tft0.options().setReal("Eta",0);
            tft0.options().setSwitch("EtaTangentOnly",true);
            const gsMatrix<> S = tft.eval3D_stress(0,pt,z,VN);
            CHECK(S.cwiseAbs().maxCoeff() > 0);
            CHECK_EQUAL(T(0),(tft0.eval3D_stress(0,pt,z,VN)-S).cwiseAbs().maxCoeff());
            CHECK_EQUAL(T(0),(tft0.eval3D_matrix(0,pt,z,MA)-tft.eval3D_matrix(0,pt,z,MA)).cwiseAbs().maxCoeff());
            CHECK_EQUAL(T(0),(tft0.eval3D_strain(0,pt,z)-tft.eval3D_strain(0,pt,z)).cwiseAbs().maxCoeff());
            return;
        }

        // Tangent (i,d) = dS_i/dE_d with the engineering strain E = (E11, E22, 2 E12)
        const gsMatrix<> C0 = tft.eval3D_matrix(0,pt,z,MA).reshape(3,3);
        const gsMatrix<> S0 = tft.eval3D_stress(0,pt,z,VN);
        const gsMatrix<> Cb = materialMatrix->eval3D_matrix(0,pt,z,MA).reshape(3,3);

        const T eta = 0.1;
        tft.options().setReal("Eta",eta);
        const gsMatrix<> J_MM = tft.eval3D_matrix(0,pt,z,MA).reshape(3,3);
        const gsMatrix<> S1 = tft.eval3D_stress(0,pt,z,VN);
        const real_t tol = 1e-6 * J_MM.cwiseAbs().maxCoeff();

        // The blend is visible at this tolerance, and the stress differs from the Eta = 0 one
        CHECK(eta * (Cb - C0).cwiseAbs().maxCoeff() >= 100 * tol);
        CHECK((S1 - S0).cwiseAbs().maxCoeff() > 1e-3 * S0.cwiseAbs().maxCoeff());

        // Deformed metric G' = G + dG as the upper Cholesky factor U (U^T U = G'), patch coefs * U^T; the decorator
        // depends on the deformed geometry through the metric only. Stress of the decorator itself.
        gsMatrix<> Gdef = Jdef.transpose() * Jdef;
        // Perturbed geometry, kept alive for as long as the decorator references it
        gsMultiPatch<> m;
        auto stressAt = [&](const gsMatrix<> & G, bool & wrinkled) -> gsMatrix<>
        {
            gsMatrix<> U(2,2);
            U << math::sqrt(G(0,0)), G(0,1)/math::sqrt(G(0,0)),
                 0.0,                math::sqrt(G(1,1) - G(0,1)*G(0,1)/G(0,0));
            m.clear();
            m.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
            m.patch(0).coefs() = m.patch(0).coefs() * U.transpose();
            m.addAutoBoundaries();
            tft.setDeformed(&m);
            wrinkled = (tft.eval3D_tensionfield(0,pt,z,MaterialOutput::Generic)(0,0) == 0);
            return tft.eval3D_stress(0,pt,z,VN);
        };

        bool w;
        const gsMatrix<> Schol = stressAt(Gdef,w);
        CHECK(w);
        CHECK_MATRIX_CLOSE(S1,Schol,1e-10*S1.cwiseAbs().maxCoeff());

        // Fourth-order central stencil (-S(2h) + 8 S(h) - 8 S(-h) + S(-2h)) / (12 h), truncation O(h^4)
        const real_t h = 1e-5;
        gsMatrix<> J_FD(3,3);
        for (index_t d=0; d!=3; ++d)
        {
            const real_t scale = (d==2) ? 1.0 : 2.0;   // dG11 = 2 dE_0, dG22 = 2 dE_1, dG12 = dE_2
            gsMatrix<> dG(2,2);
            dG.setZero();
            gsMatrix<> Sst[4];
            const real_t m[4] = {-2.0, -1.0, 1.0, 2.0};
            for (index_t s=0; s!=4; ++s)
            {
                if (d==0)      dG(0,0) = scale*m[s]*h;
                else if (d==1) dG(1,1) = scale*m[s]*h;
                else           dG(0,1) = dG(1,0) = scale*m[s]*h;
                Sst[s] = stressAt(Gdef + dG,w);
                CHECK(w);
            }
            J_FD.col(d) = (-Sst[3] + 8*Sst[2] - 8*Sst[1] + Sst[0]) / (12*h);
        }
        tft.setDeformed(&mp_def);
        CHECK_MATRIX_CLOSE(J_FD,J_MM,tol);
    }

    TEST(TFT_eta_tangent_FD_SvK_linear)
    {
        TFT_eta_affine_CHECK<true>(false);
    }

    TEST(TFT_eta_tangent_FD_NH_incomp)
    {
        TFT_eta_affine_CHECK<false>(false);
    }

    TEST(TFT_eta_zero_is_default)
    {
        const MaterialOutput MA = MaterialOutput::MatrixA;
        const MaterialOutput VN = MaterialOutput::VectorN;
        const real_t st[3][3] = { {1.2, 0.8/math::sqrt(1.2), 0}, {1.2, 0.97, 1}, {0.95, 0.90, -1} };
        for (index_t i=0; i!=3; ++i)
        {
            TFT_curved_setup s(st[i][0],st[i][1]);
            gsMaterialMatrixTFT<3,real_t,true> tft0(s.wrapped.get());
            tft0.setUndeformed(&s.ref);
            tft0.setDeformed(&s.def);
            tft0.options().setReal("Eta",0);
            tft0.options().setSwitch("EtaTangentOnly",true);
            CHECK_EQUAL((index_t)st[i][2],(index_t)s.tft->eval3D_tensionfield(0,s.u,s.z)(0,0));
            CHECK_EQUAL(0.0,(tft0.eval3D_stress(0,s.u,s.z,VN)-s.tft->eval3D_stress(0,s.u,s.z,VN)).cwiseAbs().maxCoeff());
            CHECK_EQUAL(0.0,(tft0.eval3D_matrix(0,s.u,s.z,MA)-s.tft->eval3D_matrix(0,s.u,s.z,MA)).cwiseAbs().maxCoeff());
            CHECK_EQUAL(0.0,(tft0.eval3D_strain(0,s.u,s.z)-s.tft->eval3D_strain(0,s.u,s.z)).cwiseAbs().maxCoeff());
            CHECK_EQUAL(0.0,(tft0.eval3D_pstrain(0,s.u,s.z)-s.tft->eval3D_pstrain(0,s.u,s.z)).cwiseAbs().maxCoeff());
        }
        TFT_eta_affine_CHECK<false>(true);
        TFT_eta_affine_CHECK<true>(true);
    }

    /// Curved d = 3 setup of TFT_curved_setup with an incompressible Neo-Hookean wrapped law and the nonlinear
    /// decorator (linear = false), whose wrapped energy and wrinkle state depend on the through-thickness position.
    struct TFT_energy_curved_NH_setup
    {
        typedef real_t T;

        gsMultiPatch<T> ref, def;
        gsConstantFunction<T> t, E, nu;
        gsMaterialMatrixBase<T>::uPtr wrapped;
        memory::unique_ptr< gsMaterialMatrixTFT<3,T,false> > tft;
        gsMatrix<T> u, z;

        TFT_energy_curved_NH_setup(const T stretch1, const T stretch2)
        :
        t(0.01,3), E(1e5,3), nu(0.5,3)
        {
            const T bump = 0.15;
            ref.addPatch(gsNurbsCreator<T>::BSplineSquare(1));
            ref.patch(0).degreeElevate(1);
            ref.embed(3);
            def = ref;
            gsMatrix<T> & c = def.patch(0).coefs();
            for (index_t i=0; i!=c.rows(); ++i)
            {
                const T x = c(i,0), y = c(i,1);
                c(i,0) = stretch1 * x;
                c(i,1) = stretch2 * y;
                c(i,2) = bump * (x*(1-x) + y*(1-y) + 2*x*y);
            }
            ref.addAutoBoundaries();
            def.addAutoBoundaries();

            std::vector<gsFunctionSet<T>*> par(2);
            par[0] = &E;
            par[1] = &nu;
            gsOptionList o;
            o.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
            o.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
            o.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            wrapped = getMaterialMatrix<3,T>(ref,t,par,o);

            tft = memory::make_unique(new gsMaterialMatrixTFT<3,T,false>(wrapped.get()));
            gsMaterialMatrixBase<T> * all[2] = { tft.get(), wrapped.get() };
            for (index_t i=0; i!=2; ++i)
            {
                all[i]->setUndeformed(&ref);
                all[i]->setDeformed(&def);
            }

            u.resize(2,3);
            u << 0.3, 0.7, 0.5,
                 0.4, 0.6, 0.2;
            z.resize(4,3);
            z << -0.5,  -0.4,  -0.3,
                  0.0,   0.1,  -0.1,
                  0.25,  0.3,   0.2,
                  0.5,   0.45,  0.35;
        }
    };

    /// Flat d = 3 SvK patch whose in-plane stretch varies along the first parametric direction, so that three
    /// in-plane points of one patch are slack (x = 0.1: stretches 0.95, 0.90), taut (x = 0.5: 1.2, 0.97) and
    /// wrinkled (x = 0.9: 1.45, 0.73). The deformation is x' = a(x), y' = b(x) y with quadratic a, b (exactly
    /// representable by the degree-2 patch), so the stretches are the derivatives a'(x) and b(x) up to a small
    /// shear y b'(x).
    struct TFT_energy_mixed_setup
    {
        typedef real_t T;

        gsMultiPatch<T> ref, def;
        gsConstantFunction<T> t, E, nu;
        gsMaterialMatrixBase<T>::uPtr wrapped;
        memory::unique_ptr< gsMaterialMatrixTFT<3,T,true> > tft;
        gsMatrix<T> u, z;

        /// Control point of the quadratic q(x) = q0 + q1 x + q2 x^2 on [0,1] at the Greville abscissa \a xg
        static T quadCtrl(const T q0, const T q1, const T q2, const T xg)
        {
            const T P0 = q0, P2 = q0 + q1 + q2;
            const T qh = q0 + T(0.5)*q1 + T(0.25)*q2;
            if (xg == 0) return P0;
            if (xg == 1) return P2;
            return 2*qh - T(0.5)*(P0 + P2);
        }

        TFT_energy_mixed_setup()
        :
        t(0.01,3), E(1e5,3), nu(0.3,3)
        {
            ref.addPatch(gsNurbsCreator<T>::BSplineSquare(1));
            ref.patch(0).degreeElevate(1);
            ref.embed(3);
            def = ref;
            gsMatrix<T> & c = def.patch(0).coefs();
            for (index_t i=0; i!=c.rows(); ++i)
            {
                const T x = c(i,0), y = c(i,1);
                c(i,0) = quadCtrl(0,0.8875,0.3125,x);
                c(i,1) = quadCtrl(0.8340625,0.75625,-0.96875,x) * y;
                c(i,2) = 0;
            }
            ref.addAutoBoundaries();
            def.addAutoBoundaries();

            std::vector<gsFunctionSet<T>*> par(2);
            par[0] = &E;
            par[1] = &nu;
            gsOptionList o;
            o.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            o.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            wrapped = getMaterialMatrix<3,T>(ref,t,par,o);

            tft = memory::make_unique(new gsMaterialMatrixTFT<3,T,true>(wrapped.get()));
            gsMaterialMatrixBase<T> * all[2] = { tft.get(), wrapped.get() };
            for (index_t i=0; i!=2; ++i)
            {
                all[i]->setUndeformed(&ref);
                all[i]->setDeformed(&def);
            }

            u.resize(2,3);
            u << 0.1, 0.5, 0.9,
                 0.05, 0.05, 0.05;
            z.resize(4,3);
            z << -0.5,  -0.4,  -0.3,
                  0.0,   0.1,  -0.1,
                  0.25,  0.3,   0.2,
                  0.5,   0.45,  0.35;
        }
    };

    /// Sets the three options that enter the slack factor of the decorator.
    template<class S>
    static void TFT_energy_options(S & s, const real_t eta, const real_t slackMultiplier, const bool tangentOnly)
    {
        s.tft->options().setReal("Eta",eta);
        s.tft->options().setReal("SlackMultiplier",slackMultiplier);
        s.tft->options().setSwitch("EtaTangentOnly",tangentOnly);
    }

    /// Asserts that every column of the evaluation grid of \a s is in the wrinkle state \a expectedState
    template<class S>
    static void TFT_energy_state_CHECK(const S & s, const index_t expectedState)
    {
        const gsMatrix<real_t> TF = s.tft->eval3D_tensionfield(0,s.u,s.z);
        CHECK_EQUAL(s.u.cols()*s.z.rows(),TF.cols());
        for (index_t i=0; i!=TF.cols(); ++i)
            CHECK_EQUAL(expectedState,(index_t)TF(0,i));
    }

    /// Taut state: the decorator energy is the wrapped membrane energy (VectorN output), column for column and
    /// independent of the residual stiffness options.
    template<class S>
    static void TFT_energy_taut_CHECK(S & s)
    {
        typedef real_t T;
        TFT_energy_state_CHECK(s,1);

        const gsMatrix<T> Ww = s.wrapped->eval3D_energy(0,s.u,s.z,MaterialOutput::VectorN);
        CHECK_EQUAL(1,Ww.rows());
        CHECK_EQUAL(s.u.cols()*s.z.rows(),Ww.cols());
        CHECK(Ww.cwiseAbs().maxCoeff() > 0);

        const gsMatrix<T> W0 = s.tft->eval3D_energy(0,s.u,s.z);
        CHECK_EQUAL(1,W0.rows());
        CHECK_EQUAL(Ww.cols(),W0.cols());
        CHECK_EQUAL(T(0),(W0-Ww).cwiseAbs().maxCoeff());

        TFT_energy_options(s,0.1,1e-3,false);
        CHECK_EQUAL(T(0),(s.tft->eval3D_energy(0,s.u,s.z)-Ww).cwiseAbs().maxCoeff());
        TFT_energy_options(s,0.1,0,true);
        CHECK_EQUAL(T(0),(s.tft->eval3D_energy(0,s.u,s.z,MaterialOutput::Generic)-Ww).cwiseAbs().maxCoeff());
    }

    /// Slack state: W = sigma W_wrapped, and the stress carries the same factor sigma, which ties the energy to the
    /// stress path (the tangent uses the unswitched Eta). The slack factor of the decorator is
    /// SlackMultiplier if Eta = 0 or SlackMultiplier > 0 and Eta otherwise, where EtaTangentOnly sets the Eta of the
    /// stress path to 0.
    template<class S>
    static void TFT_energy_slack_CHECK(S & s)
    {
        typedef real_t T;
        const MaterialOutput VN = MaterialOutput::VectorN;
        TFT_energy_state_CHECK(s,-1);

        const gsMatrix<T> Ww = s.wrapped->eval3D_energy(0,s.u,s.z,VN);
        const gsMatrix<T> Sw = s.wrapped->eval3D_stress(0,s.u,s.z,VN);
        CHECK(Ww.cwiseAbs().maxCoeff() > 0);
        CHECK(Sw.cwiseAbs().maxCoeff() > 0);

        struct Case { T eta, slackMultiplier; bool tangentOnly; T sigma; };
        const Case cases[6] = {
            { 0.0, 0.0,  false, 0.0  },   // no slack factor: zero energy
            { 0.1, 0.0,  false, 0.1  },   // sigma = Eta
            { 0.0, 1e-3, false, 1e-3 },   // sigma = SlackMultiplier
            { 0.1, 1e-3, false, 1e-3 },   // SlackMultiplier overrides Eta
            { 0.1, 0.0,  true,  0.0  },   // EtaTangentOnly: Eta taken as 0
            { 0.1, 1e-3, true,  1e-3 } }; // EtaTangentOnly, SlackMultiplier > 0
        for (index_t i=0; i!=6; ++i)
        {
            TFT_energy_options(s,cases[i].eta,cases[i].slackMultiplier,cases[i].tangentOnly);
            const gsMatrix<T> W = s.tft->eval3D_energy(0,s.u,s.z);
            const gsMatrix<T> St = s.tft->eval3D_stress(0,s.u,s.z,VN);
            CHECK_EQUAL(1,W.rows());
            CHECK_EQUAL(Ww.cols(),W.cols());
            if (cases[i].sigma == 0)
            {
                CHECK_EQUAL(T(0),W.cwiseAbs().maxCoeff());
                CHECK_EQUAL(T(0),St.cwiseAbs().maxCoeff());
            }
            else
            {
                TFT_eta_CLOSE(cases[i].sigma*Ww,W,1e-14);
                TFT_eta_CLOSE(cases[i].sigma*Sw,St,1e-14);
            }
        }
    }

    /// Wrinkled state: the energy is NaN in every column, with and without residual stiffness, while the wrapped
    /// energy is finite there.
    template<class S>
    static void TFT_energy_wrinkled_CHECK(S & s)
    {
        typedef real_t T;
        TFT_energy_state_CHECK(s,0);

        const gsMatrix<T> Ww = s.wrapped->eval3D_energy(0,s.u,s.z,MaterialOutput::VectorN);
        for (index_t i=0; i!=Ww.cols(); ++i)
            CHECK(!(std::isnan)(Ww(0,i)));

        const T etas[2] = { 0.0, 0.1 };
        for (index_t e=0; e!=2; ++e)
        {
            TFT_energy_options(s,etas[e],0,false);
            const gsMatrix<T> W = s.tft->eval3D_energy(0,s.u,s.z);
            CHECK_EQUAL(1,W.rows());
            CHECK_EQUAL(Ww.cols(),W.cols());
            for (index_t i=0; i!=W.cols(); ++i)
                CHECK((std::isnan)(W(0,i)));
        }
    }

    /// Output layout 1 x (n nz), column j*n+k = point k at row j of z: the block call equals the calls with a single
    /// in-plane point and the matching column of z.
    template<class S>
    static void TFT_energy_shape_CHECK(S & s)
    {
        typedef real_t T;
        const index_t n = s.u.cols(), nz = s.z.rows();
        const gsMatrix<T> W = s.tft->eval3D_energy(0,s.u,s.z);
        CHECK_EQUAL(1,W.rows());
        CHECK_EQUAL(n*nz,W.cols());
        const T tol = 1e-14 * W.cwiseAbs().maxCoeff();
        CHECK(W.cwiseAbs().maxCoeff() > 0);
        for (index_t k=0; k!=n; ++k)
        {
            const gsMatrix<T> uk = s.u.col(k), zk = s.z.col(k);
            const gsMatrix<T> W1 = s.tft->eval3D_energy(0,uk,zk);
            CHECK_EQUAL(1,W1.rows());
            CHECK_EQUAL(nz,W1.cols());
            for (index_t j=0; j!=nz; ++j)
                CHECK_CLOSE(W1(0,j),W(0,j*n+k),tol);
        }
    }

    TEST(TFT_energy_taut_linear)
    {
        TFT_curved_setup s(1.2,0.97);
        TFT_energy_taut_CHECK(s);
    }

    TEST(TFT_energy_slack_linear)
    {
        TFT_curved_setup s(0.95,0.90);
        TFT_energy_slack_CHECK(s);
    }

    TEST(TFT_energy_wrinkled_linear)
    {
        TFT_curved_setup s(1.2,0.8/math::sqrt(1.2));
        TFT_energy_wrinkled_CHECK(s);
    }

    TEST(TFT_energy_taut_nonlinear)
    {
        TFT_energy_curved_NH_setup s(1.2,0.97);
        TFT_energy_taut_CHECK(s);
    }

    TEST(TFT_energy_slack_nonlinear)
    {
        TFT_energy_curved_NH_setup s(0.95,0.90);
        TFT_energy_slack_CHECK(s);
    }

    TEST(TFT_energy_wrinkled_nonlinear)
    {
        TFT_energy_curved_NH_setup s(1.2,0.8/math::sqrt(1.2));
        TFT_energy_wrinkled_CHECK(s);
    }

    /// One call with three in-plane points in the slack, taut and wrinkled state: each column is checked against
    /// its own expectation (sigma W_wrapped, W_wrapped, NaN).
    TEST(TFT_energy_mixed_states_linear)
    {
        typedef real_t T;
        TFT_energy_mixed_setup s;
        const index_t n = s.u.cols(), nz = s.z.rows();
        const T sigma = 0.1;
        s.tft->options().setReal("Eta",sigma);

        const gsMatrix<T> TF = s.tft->eval3D_tensionfield(0,s.u,s.z);
        const gsMatrix<T> Ww = s.wrapped->eval3D_energy(0,s.u,s.z,MaterialOutput::VectorN);
        const gsMatrix<T> W  = s.tft->eval3D_energy(0,s.u,s.z);
        CHECK_EQUAL(1,W.rows());
        CHECK_EQUAL(n*nz,W.cols());
        const index_t state[3] = { -1, 1, 0 };
        for (index_t j=0; j!=nz; ++j)
        {
            for (index_t k=0; k!=n; ++k)
            {
                const index_t c = j*n+k;
                CHECK_EQUAL(state[k],(index_t)TF(0,c));
                CHECK(Ww(0,c) > 0);
                if (state[k] == 1)
                    CHECK_EQUAL(T(0),math::abs(W(0,c)-Ww(0,c)));
                else if (state[k] == -1)
                    CHECK_CLOSE(sigma*Ww(0,c),W(0,c),1e-14*Ww(0,c));
                else
                    CHECK((std::isnan)(W(0,c)));
            }
        }
    }

    /// Output shape and per-point consistency; the linear wrapped energy is the membrane one and so does not vary
    /// through the thickness.
    TEST(TFT_energy_shape_linear)
    {
        TFT_curved_setup s(1.2,0.97);
        TFT_energy_state_CHECK(s,1);
        TFT_energy_shape_CHECK(s);
        const gsMatrix<real_t> W = s.tft->eval3D_energy(0,s.u,s.z);
        CHECK_EQUAL(real_t(0),TFT_zVariation(W,s.u.cols(),s.z.rows()));
    }

    /// As TFT_energy_shape_linear; for the Neo-Hookean wrapped law the energy depends on z, which makes the
    /// z-ordering of the columns observable. In the taut state the decorator output is the wrapped one column
    /// for column.
    TEST(TFT_energy_shape_nonlinear)
    {
        typedef real_t T;
        TFT_energy_curved_NH_setup s(1.2,0.97);
        TFT_energy_state_CHECK(s,1);
        TFT_energy_shape_CHECK(s);
        const gsMatrix<T> W  = s.tft->eval3D_energy(0,s.u,s.z);
        const gsMatrix<T> Ww = s.wrapped->eval3D_energy(0,s.u,s.z,MaterialOutput::VectorN);
        CHECK_EQUAL(T(0),(W-Ww).cwiseAbs().maxCoeff());
        CHECK(TFT_zVariation(W,s.u.cols(),s.z.rows()) > 0);
    }

    /// Directional derivative of the energy along a metric path G(t) = G0 + t dG at fixed reference metric I:
    /// with the engineering Green-Lagrange strain E = (E11, E22, 2 E12) = ((G11-1)/2, (G22-1)/2, G12) one has
    /// dW/dt = S . dE/dt = S . (dG11/2, dG22/2, dG12), S the VectorN stress of the decorator. Fourth-order central
    /// differences of the decorator energy at z = 0 are compared with the contraction, in the taut state and in
    /// the slack state (with a nonzero slack factor, otherwise both sides vanish). The state is asserted at the
    /// base point and at every stencil point.
    template<bool linear>
    static void TFT_energy_FD_CHECK(const bool slack)
    {
        typedef real_t T;
        const MaterialOutput VN = MaterialOutput::VectorN;
        const index_t expectedState = slack ? -1 : 1;

        const T al = 0.9;
        gsMatrix<T> R(2,2);
        R << math::cos(al), -math::sin(al), math::sin(al), math::cos(al);
        gsMatrix<T> D(2,2);
        D << (slack ? 0.95 : 1.2), 0.0, 0.0, (slack ? 0.90 : 0.97);
        const gsMatrix<T> Jdef = R * D * R.transpose();
        const gsMatrix<T> G0 = Jdef.transpose() * Jdef;

        gsMultiPatch<T> mp, mp_def, m;
        mp.addPatch( gsNurbsCreator<T>::BSplineSquare(1) );
        mp.addAutoBoundaries();
        mp_def.addPatch( gsNurbsCreator<T>::BSplineSquare(1) );
        mp_def.patch(0).coefs() = mp_def.patch(0).coefs() * Jdef.transpose();
        mp_def.addAutoBoundaries();

        const T mu = 4.225e5;
        const T PoissonRatio = linear ? 0.45 : 0.5;
        gsConstantFunction<T> t(0.01,2), E(2*mu*(1+PoissonRatio),2), nu(PoissonRatio,2);
        std::vector<gsFunctionSet<T>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;
        gsOptionList options;
        if (linear)
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        }
        gsMaterialMatrixBase<T>::uPtr materialMatrix = getMaterialMatrix<2,T>(mp,t,parameters,options);
        materialMatrix->setDeformed(&mp_def);
        gsMaterialMatrixTFT<2,T,linear> tft(*materialMatrix);
        tft.setDeformed(&mp_def);

        gsVector<T> pt(2);
        pt.setConstant(0.25);
        gsMatrix<T> z(1,1);
        z.setZero();

        // Deformed metric G as the upper Cholesky factor U (U^T U = G): patch coefs * U^T has metric G exactly
        auto setMetric = [&](const gsMatrix<T> & G)
        {
            gsMatrix<T> U(2,2);
            U << math::sqrt(G(0,0)), G(0,1)/math::sqrt(G(0,0)),
                 0.0,                math::sqrt(G(1,1) - G(0,1)*G(0,1)/G(0,0));
            m.clear();
            m.addPatch( gsNurbsCreator<T>::BSplineSquare(1) );
            m.patch(0).coefs() = m.patch(0).coefs() * U.transpose();
            m.addAutoBoundaries();
            tft.setDeformed(&m);
        };
        auto energyAt = [&](const gsMatrix<T> & G) -> T
        {
            setMetric(G);
            CHECK_EQUAL(expectedState,(index_t)tft.eval3D_tensionfield(0,pt,z)(0,0));
            return tft.eval3D_energy(0,pt,z)(0,0);
        };

        gsMatrix<T> dG(2,2);
        dG << 0.3, 0.15, 0.15, -0.2;
        gsVector<T> dEhat(3);
        dEhat << 0.5*dG(0,0), 0.5*dG(1,1), dG(0,1);

        // Slack factors sigma of the two ways to make the slack energy nonzero
        const index_t nsets = slack ? 2 : 1;
        const T etas[2] = { 0.1, 0.0 }, slackMults[2] = { 0.0, 1e-3 };
        const T h = 1e-5;
        for (index_t i=0; i!=nsets; ++i)
        {
            tft.options().setReal("Eta",slack ? etas[i] : 0.0);
            tft.options().setReal("SlackMultiplier",slack ? slackMults[i] : 0.0);

            const T W0 = energyAt(G0);
            const gsMatrix<T> S = tft.eval3D_stress(0,pt,z,VN);
            CHECK_EQUAL(3,S.rows());
            CHECK(math::abs(W0) > 0);

            const T dWdt = S.col(0).dot(dEhat);
            const T scale = S.cwiseAbs().maxCoeff() * dEhat.cwiseAbs().maxCoeff();
            CHECK(math::abs(dWdt) > 1e-3 * scale);

            // Fourth-order central stencil (-W(2h) + 8 W(h) - 8 W(-h) + W(-2h)) / (12 h), truncation O(h^4)
            const T mult[4] = {-2.0, -1.0, 1.0, 2.0};
            T Wst[4];
            for (index_t s=0; s!=4; ++s)
                Wst[s] = energyAt(G0 + mult[s]*h*dG);
            const T dWdt_FD = (-Wst[3] + 8*Wst[2] - 8*Wst[1] + Wst[0]) / (12*h);
            CHECK_CLOSE(dWdt,dWdt_FD,1e-6*scale);
        }
        tft.setDeformed(&mp_def);
    }

    TEST(TFT_energy_FD_linear)
    {
        TFT_energy_FD_CHECK<true>(false);
        TFT_energy_FD_CHECK<true>(true);
    }

    TEST(TFT_energy_FD_nonlinear)
    {
        TFT_energy_FD_CHECK<false>(false);
        TFT_energy_FD_CHECK<false>(true);
    }

}
