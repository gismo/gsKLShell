/** @file gsThinShellAssemblerWrinkle_test.cpp

    @brief Unit tests of gsThinShellAssembler<2,real_t,false> with the exact-wrinkling law gsMaterialMatrixTFTSplit
    (Exact on) on a homogeneous, fully wrinkled membrane patch.

    The unit square carries an incompressible Neo-Hookean membrane (mu = 1.5e6, thickness 1e-3). The boundary
    displacement is the affine field u = ((l1-1) X, (l2-1) Y) with l1 = 1.5 and l2 = 0.6. The free lateral contraction
    of uniaxial incompressible tension is l1^(-1/2) = 0.816 > l2, so the membrane is wrinkled at every point and the
    affine field is the exact equilibrium (zero load). The oracles are closed forms, with G = I and the wrinkle
    direction n = (0,1):
    - wrinkle intensity    gamma = (1/l1 - l2^2)/2,
    - relaxed stress       S' = (mu (1 - l1^-3), 0, 0),
    - unrelaxed stress     S(E) = mu (1 - 1/(l1^4 l2^2), 1 - 1/(l1^2 l2^4), 0),
    - returned stress      S_eff = (1 - Eta) S' + Eta S(E),
    - reaction on the east side (-t l1 S_eff^11, 0).

    Checks: the wrinkled state, the material point, the symmetry and the finite-difference consistency of the tangent,
    Newton convergence with a direct (LU) and an iterative (CGDiagonal) linear solver, the reaction force, and the
    fallback counters of the material.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <random>
#include <vector>

SUITE(gsThinShellAssemblerWrinkle_test)                 // The suite should have the same name as the file
{
    static const real_t WPmu  = 1.5e6;
    static const real_t WPnu  = 0.5;
    static const real_t WPt   = 1e-3;
    static const real_t WPeta = 1e-6;
    static const real_t WPl1  = 1.5;
    static const real_t WPl2  = 0.6;

    /// Wrinkle intensity of the exact solution
    static real_t WPgammaRef()
    {
        return 0.5 * (1.0/WPl1 - WPl2*WPl2);
    }

    /// Contravariant (S11,S22,S12) returned by the law: (1-Eta) S' + Eta S(E)
    static gsVector<real_t,3> WPstressRef()
    {
        gsVector<real_t,3> relaxed, unrelaxed;
        relaxed   << WPmu * (1.0 - 1.0/(WPl1*WPl1*WPl1)), 0.0, 0.0;
        unrelaxed << WPmu * (1.0 - 1.0/(math::pow(WPl1,4)*WPl2*WPl2)),
                     WPmu * (1.0 - 1.0/(WPl1*WPl1*math::pow(WPl2,4))),
                     0.0;
        return (1.0 - WPeta) * relaxed + WPeta * unrelaxed;
    }

    /// 1x1 through-thickness coordinate matrix
    static gsMatrix<real_t> WPZ(const real_t z)
    {
        gsMatrix<real_t> zmat(1,1);
        zmat(0,0) = z;
        return zmat;
    }

    /// Wrinkling law options of the benchmark: exact solve on, tension-field criterion 2
    static void WPconfigure(gsMaterialMatrixTFTSplit<2,real_t,false> * split)
    {
        gsOptionList o = split->options();
        o.setReal("Eta",WPeta);
        o.setSwitch("Exact",true);
        o.setInt("TensionField",2);
        split->setOptions(o);
    }

    /// Incompressible Neo-Hookean base material of the patch
    static gsMaterialMatrixBase<real_t>::uPtr WPbase(const gsMultiPatch<real_t> & mp,
                                                      const gsFunctionSet<real_t> & t,
                                                      std::vector<gsFunctionSet<real_t>*> & parameters)
    {
        gsOptionList options;
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        return getMaterialMatrix<2,real_t>(mp,t,parameters,options);
    }

    static gsMultiPatch<real_t> WPgeometry()
    {
        gsMultiPatch<real_t> mp;
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        mp.degreeElevate(1);
        mp.uniformRefine();
        mp.uniformRefine();
        return mp;
    }

    /// Residual history of a Newton run
    struct WPNewton
    {
        std::vector<real_t> rnorm;   ///< |r_k|, k = 0..iterations
        gsVector<real_t>    x;       ///< last iterate
        index_t             iterations;
        bool                converged;
    };

    /**
     * Wrinkled patch: geometry, bases, boundary conditions, the base material and the wrinkling law integrated by the
     * assembler, and a second law on a separate base material for pointwise evaluations. The parameter functions and
     * the boundary functions are members because the material and the boundary conditions keep references to them.
     */
    struct WrinklePatch
    {
        gsMultiPatch<real_t> mp, mp_def, dispAffine, dispPerturbed;
        gsMultiBasis<real_t> dbasis;
        gsFunctionExpr<real_t> ux, uy;
        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> t, E, nu, force;
        std::vector<gsFunctionSet<real_t>*> parameters;
        gsMaterialMatrixBase<real_t>::uPtr base, base2;
        memory::unique_ptr<gsMaterialMatrixTFTSplit<2,real_t,false> > split, split2;
        memory::unique_ptr<gsThinShellAssembler<2,real_t,false> > assembler;
        gsMatrix<real_t> nodes;

        WrinklePatch()
        :
        mp(WPgeometry()), dbasis(mp),
        ux("0.5*x",2), uy("-0.4*y",2),
        t(WPt,2), E(2*WPmu*(1+WPnu),2), nu(WPnu,2),
        force(gsVector<real_t>::Zero(2),2)
        {
            mp_def = mp;

            // Fixed affine displacement and its fixed perturbation, as multipatch displacements
            dispAffine = mp;
            dispAffine.patch(0).coefs().col(0) *= (WPl1 - 1.0);
            dispAffine.patch(0).coefs().col(1) *= (WPl2 - 1.0);

            dispPerturbed = dispAffine;
            const real_t pi = static_cast<real_t>(EIGEN_PI);
            for (index_t i = 0; i != mp.patch(0).coefs().rows(); ++i)
            {
                const real_t X = mp.patch(0).coefs()(i,0), Y = mp.patch(0).coefs()(i,1);
                const real_t bump = 1e-3 * math::sin(pi*X) * math::sin(pi*Y);
                dispPerturbed.patch(0).coefs()(i,0) += bump;
                dispPerturbed.patch(0).coefs()(i,1) += bump;
            }

            bc.setGeoMap(mp);
            const boundary::side sides[4] = {boundary::west, boundary::east, boundary::south, boundary::north};
            for (index_t s = 0; s != 4; ++s)
            {
                bc.addCondition(sides[s], condition_type::dirichlet, &ux, 0, false, 0);
                bc.addCondition(sides[s], condition_type::dirichlet, &uy, 0, false, 1);
            }

            parameters.push_back(&E);
            parameters.push_back(&nu);
            base  = WPbase(mp,t,parameters);
            base2 = WPbase(mp,t,parameters);

            split  = memory::make_unique(new gsMaterialMatrixTFTSplit<2,real_t,false>(base.get()));
            split2 = memory::make_unique(new gsMaterialMatrixTFTSplit<2,real_t,false>(base2.get()));
            WPconfigure(split.get());
            WPconfigure(split2.get());
            split2->setUndeformed(&mp);

            // The assembler integrates the given law without cloning it, so its counters can be read afterwards
            assembler = memory::make_unique(new gsThinShellAssembler<2,real_t,false>(mp,dbasis,bc,force,split.get()));

            gsOptionList quOpt;
            quOpt.addReal("quA","",1.0);
            quOpt.addInt("quB","",1);
            nodes = gsQuadrature::getAllNodes(dbasis.basis(0),quOpt);
        }

        /// Deformed geometry x = (l1 X, l2 Y)
        gsMultiPatch<real_t> affineDeformed() const
        {
            gsMultiPatch<real_t> d = mp;
            d.patch(0).coefs().col(0) *= WPl1;
            d.patch(0).coefs().col(1) *= WPl2;
            return d;
        }

        /// Stiffness matrix at the free-dof vector \a x
        gsSparseMatrix<real_t> jacobian(const gsVector<real_t> & x)
        {
            assembler->constructSolution(x,mp_def);
            ThinShellAssemblerStatus status = assembler->assembleMatrix(mp_def);
            GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly of the matrix failed");
            gsSparseMatrix<real_t> m = assembler->matrix();
            return m;
        }

        /// Residual F_ext - F_int at the free-dof vector \a x
        gsVector<real_t> residual(const gsVector<real_t> & x)
        {
            assembler->constructSolution(x,mp_def);
            ThinShellAssemblerStatus status = assembler->assembleVector(mp_def);
            GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly of the vector failed");
            gsVector<real_t> r = assembler->rhs();
            return r;
        }

        /// Newton iteration x += K^{-1} r with the given linear solve; \a solve(K,r,k) returns the update of step k
        WPNewton newton(const gsVector<real_t> & x0,
                        const std::function<gsVector<real_t>(const gsSparseMatrix<real_t> &, const gsVector<real_t> &, index_t)> & solve,
                        const std::string & label)
        {
            const index_t maxIt = 20;
            WPNewton res;
            res.x = x0;
            res.converged = false;
            res.iterations = 0;

            gsVector<real_t> r = residual(res.x);
            res.rnorm.push_back(r.norm());
            const real_t r0 = res.rnorm[0];
            gsInfo << std::setprecision(16) << std::scientific;

            for (index_t k = 0; ; ++k)
            {
                gsInfo << label << " k=" << k << " |r|=" << res.rnorm[k] << " rel=" << res.rnorm[k]/r0 << "\n";
                if (res.rnorm[k] <= 1e-10 * r0)
                {
                    res.converged = true;
                    break;
                }
                if (k == maxIt)
                    break;

                gsSparseMatrix<real_t> K = jacobian(res.x);
                res.x += solve(K,r,k);
                ++res.iterations;
                r = residual(res.x);
                res.rnorm.push_back(r.norm());
            }
            return res;
        }
    };

    TEST(Wrinkled_state_and_dof_count)
    {
        WrinklePatch wp;
        CHECK_EQUAL(32,wp.assembler->numDofs());

        // Perturbed state: the minor stress is strongly negative, so every Gauss point is wrinkled
        const gsVector<real_t> x0 = wp.assembler->constructSolutionVector(wp.dispPerturbed);
        const gsMultiPatch<real_t> pert = wp.assembler->constructSolution(x0);
        wp.split2->setDeformed(&pert);

        index_t nBadState = 0, nBadGamma = 0;
        for (index_t k = 0; k != wp.nodes.cols(); ++k)
        {
            const gsMatrix<real_t> u = wp.nodes.col(k);
            const real_t state = wp.split2->eval3D_tensionfield(0,u,WPZ(0),MaterialOutput::Generic)(0,0);
            const real_t gamma = wp.split2->eval3D_gamma(0,u,WPZ(0),MaterialOutput::Generic)(0,0);
            if (state != 0)
                ++nBadState;
            if (!(math::isfinite(gamma) && gamma > 0))
                ++nBadGamma;
        }
        gsInfo << "wrinkled-state: nodes=" << wp.nodes.cols() << " non-wrinkled=" << nBadState
               << " bad-gamma=" << nBadGamma << "\n";
        CHECK_EQUAL(0,nBadState);
        CHECK_EQUAL(0,nBadGamma);
    }

    // (i) material point on the affine, exactly wrinkled state
    TEST(Material_point_affine_state)
    {
        WrinklePatch wp;
        const gsMultiPatch<real_t> aff = wp.affineDeformed();
        wp.split2->setDeformed(&aff);

        const real_t gammaRef = WPgammaRef();
        const gsVector<real_t,3> Sref = WPstressRef();

        real_t maxGammaErr = 0, maxStressErr = 0;
        for (index_t k = 0; k != wp.nodes.cols(); ++k)
        {
            const gsMatrix<real_t> u = wp.nodes.col(k);
            const real_t gamma = wp.split2->eval3D_gamma(0,u,WPZ(0),MaterialOutput::Generic)(0,0);
            const gsMatrix<real_t> S = wp.split2->eval3D_vector(0,u,WPZ(0));
            GISMO_ENSURE(S.rows()==3 && S.cols()==1,"Unexpected stress size");
            maxGammaErr  = math::max(maxGammaErr,math::abs(gamma - gammaRef) / math::abs(gammaRef));
            maxStressErr = math::max(maxStressErr,(S.col(0) - Sref).norm() / Sref.norm());
        }
        gsInfo << std::setprecision(16) << std::scientific
               << "(i) gamma_ref=" << gammaRef << " S_ref=" << Sref.transpose() << "\n"
               << "(i) max rel gamma error=" << maxGammaErr << " max rel stress error=" << maxStressErr << "\n";
        CHECK(maxGammaErr  <= 1e-9);
        CHECK(maxStressErr <= 1e-9);
    }

    // (ii) symmetry and finite-difference consistency of the tangent at the perturbed state
    TEST(Tangent_perturbed_state)
    {
        WrinklePatch wp;
        const index_t n = wp.assembler->numDofs();
        const gsVector<real_t> x0 = wp.assembler->constructSolutionVector(wp.dispPerturbed);

        const gsSparseMatrix<real_t> K = wp.jacobian(x0);
        const gsSparseMatrix<real_t> Kt = K.transpose();
        const real_t symmetry = (K - Kt).norm() / K.norm();
        gsInfo << std::setprecision(16) << std::scientific
               << "(ii) symmetry ||K-K^T||_F/||K||_F=" << symmetry << "\n";
        CHECK(symmetry <= 1e-12);

        // Fixed directions: portable map of a seeded Mersenne twister to [-1,1]
        std::mt19937 gen(1234);
        const real_t h = 1e-6;
        for (index_t d = 0; d != 3; ++d)
        {
            gsVector<real_t> v(n);
            for (index_t i = 0; i != n; ++i)
                v[i] = 2.0 * static_cast<real_t>(gen()) / 4294967295.0 - 1.0;

            // rhs = F_ext - F_int, so central differences of rhs approximate -K v
            const gsVector<real_t> fd = -(wp.residual(x0 + h*v) - wp.residual(x0 - h*v)) / (2*h);
            const gsVector<real_t> Kv = K * v;
            const real_t err = (Kv - fd).norm() / Kv.norm();
            gsInfo << "(ii) FD direction " << d << " ||Kv-fd||/||Kv||=" << err << "\n";
            CHECK(err <= 5e-6);
        }
    }

    // (iii) Newton with LU, (iv) Newton with CGDiagonal, (v) reaction force, and the material counters
    TEST(Newton_LU_CG_reaction_counters)
    {
        WrinklePatch wp;
        const index_t n = wp.assembler->numDofs();
        const gsVector<real_t> x0   = wp.assembler->constructSolutionVector(wp.dispPerturbed);
        const gsVector<real_t> xAff = wp.assembler->constructSolutionVector(wp.dispAffine);
        wp.split->resetCounters();

        // (iii) LU
        gsSparseSolver<real_t>::LU lu;
        const WPNewton resLU = wp.newton(x0,
            [&lu](const gsSparseMatrix<real_t> & K, const gsVector<real_t> & r, index_t)
            {
                lu.compute(K);
                return gsVector<real_t>(lu.solve(r));
            },"LU");

        gsInfo << std::setprecision(16) << std::scientific
               << "(iii) LU converged=" << resLU.converged << " iterations=" << resLU.iterations << "\n";
        CHECK(resLU.converged);
        CHECK(resLU.iterations <= 8);
        for (size_t k = 0; k + 1 < resLU.rnorm.size(); ++k)
        {
            const real_t rel = resLU.rnorm[k] / resLU.rnorm[0];
            if (rel > 1e-9 && rel <= 1e-2)
            {
                const real_t ratio = resLU.rnorm[k+1] / resLU.rnorm[k];
                gsInfo << "(iii) LU contraction k=" << k << " |r_{k+1}|/|r_k|=" << ratio << "\n";
                CHECK(ratio <= 0.1);
            }
        }
        const real_t dispErrLU = (resLU.x - xAff).norm() / xAff.norm();
        gsInfo << "(iii) LU relative displacement error=" << dispErrLU << "\n";
        CHECK(dispErrLU <= 1e-8);

        // (iv) CGDiagonal with an explicit tolerance and iteration cap
        gsSparseSolver<real_t>::CGDiagonal cg;
        cg.setTolerance(1e-10);
        cg.setMaxIterations(20*n);
        const WPNewton resCG = wp.newton(x0,
            [&cg](const gsSparseMatrix<real_t> & K, const gsVector<real_t> & r, index_t k)
            {
                cg.compute(K);
                gsVector<real_t> dx = cg.solve(r);
                gsInfo << "CG k=" << k << " solve iterations=" << cg.iterations() << " error=" << cg.error() << "\n";
                return dx;
            },"CG");

        gsInfo << "(iv) CG converged=" << resCG.converged << " iterations=" << resCG.iterations
               << " (LU iterations=" << resLU.iterations << ")\n";
        CHECK(resCG.converged);
        CHECK(resCG.iterations <= resLU.iterations + 2);

        // Counters of the law inside the assembler, accumulated over both runs
        gsInfo << "counters: exactFallbacks=" << wp.split->numExactFallbacks()
               << " exactReclassified=" << wp.split->numExactReclassified()
               << " thetaFallbacks=" << wp.split->numThetaFallbacks() << "\n";
        CHECK_EQUAL(0,wp.split->numExactFallbacks());
        CHECK_EQUAL(0,wp.split->numExactReclassified());
        CHECK_EQUAL(0,wp.split->numThetaFallbacks());

        // (v) reaction on the east side at the LU-converged state
        wp.assembler->constructSolution(resLU.x,wp.mp_def);
        const patchSide ps(0,boundary::east);
        const gsMatrix<real_t> F = wp.assembler->boundaryForce(wp.mp_def,ps);
        const real_t Fref = -WPt * WPl1 * WPstressRef()[0];
        gsInfo << "(v) F=(" << F(0,0) << ", " << F(1,0) << ") F_ref=" << Fref << "\n";
        CHECK(math::abs(F(0,0) - Fref) <= 1e-6 * math::abs(Fref));
        CHECK(math::abs(F(1,0)) <= 1e-6 * math::abs(F(0,0)));
    }
}
