/** @file example_shell2D.cpp

    @brief Simple 2D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.h>

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

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

SUITE(gsMaterialMatrix_test)                 // The suite should have the same name as the file
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
     index_t mat = 1, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Comp_Generic)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Incomp_Spectral)
    {
     index_t mat = 1, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_NH_Comp_Spectral)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }

    TEST(MM_MR_Incomp_Analytical)
    {
     index_t mat = 1, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Comp_Analytical)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Incomp_Generic)
    {
     index_t mat = 1, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Comp_Generic)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Incomp_Spectral)
    {
     index_t mat = 1, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_MR_Comp_Spectral)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }

    TEST(MM_OG_Incomp_Spectral)
    {
     index_t mat = 1, impl = 1;
     bool comp = false;
     MM_CHECK(mat, impl, comp);
    }
    TEST(MM_OG_Comp_Spectral)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     MM_CHECK(mat, impl, comp);
    }

    void MM_CHECK(const index_t material, const index_t impl, const bool Compressibility)
    {
        if (material==4 && impl!=3)
            CHECK(true);

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

        gsMaterialMatrixBase<real_t> * materialMatrix;
        gsMaterialMatrixBase<real_t> * materialMatrixTFT;
        gsOptionList options;
        if      (material==0)
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
            materialMatrixTFT = new gsMaterialMatrixTFT<2,real_t,true>(materialMatrix);
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
            materialMatrixTFT = new gsMaterialMatrixTFT<2,real_t,false>(materialMatrix);
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
        Sfun<real_t> Sfun(gori,materialMatrix,0,pt,z);
        Sfun.deriv_into(e,resss);
        gsMatrix<> C_FD = resss.reshape(3,3); // Finite Differences
        gsMatrix<> C_MM = materialMatrix->eval3D_matrix(0,pt,z,MaterialOutput::Generic).reshape(3,3); // implemented
        CHECK_MATRIX_CLOSE(C_FD,C_MM,1e-3);

        /// dMATRIX
        Cfun<real_t> Cfun(gori,materialMatrix,0,pt,z);
        Cfun.deriv_into(e,resss);
        gsMatrix<> dC_FD = resss.reshape(9,3);
        gsMatrix<> dC_MM = materialMatrix->eval3D_dmatrix(0,pt,z,MaterialOutput::Generic).reshape(9,3);
        dC_MM *= 2; // NOTE: dmatrix returns d(mm)/dC, and dC_FD is d(mm)/dE = 2*d(mm)/dC
        CHECK_MATRIX_CLOSE(dC_FD,dC_MM,1e-3);

        /// STRESS TFT
        STFTfun<real_t> STFT(gori,materialMatrix,0,pt,z);
        STFT.eval_into(e,resss);
        gsMatrix<> STFT_test = resss;
        gsMatrix<> STFT_MM   = vecTFT.piece(0).eval(pt).reshape(3,1);
        CHECK_MATRIX_CLOSE(STFT_test,STFT_test,1e-3);

        /// MATRIX TFT
        STFT.deriv_into(e,resss);
        gsMatrix<> CTFT_FD = resss.reshape(3,3);
        gsMatrix<> CTFT_MM = matTFT.piece(0).eval(pt).reshape(3,3);
        CHECK_MATRIX_CLOSE(CTFT_FD,CTFT_MM,1e-3);

        delete materialMatrix;
        delete materialMatrixTFT;
    }


}