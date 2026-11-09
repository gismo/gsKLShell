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

}
