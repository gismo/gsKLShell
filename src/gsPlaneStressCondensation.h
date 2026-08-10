/** @file gsPlaneStressCondensation.h

    @brief Batched plane-stress static condensation of 3D PFF material laws.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

// The whole content is only meaningful when the gsPhaseFieldFracture module
// (which provides gsMaterialBase / gsMaterialData) is enabled. The macro is set
// in the generated gsCore/gsConfigExt.h, pulled in by any gismo header above.
#ifdef gsPhaseFieldFracture_ENABLED

#include <gsPhaseFieldFracture/materials/gsMaterialBase.h>

namespace gismo
{

/**
 * @brief   Material-agnostic, batched plane-stress condensation of a 3D
 *          gsPhaseFieldFracture (PFF) material law.
 *
 * Given a batch of in-plane strain states (one per point, expressed in a single
 * LOCAL CARTESIAN frame), this helper enforces the plane-stress condition
 * \f$ S_{33}=0 \f$ by solving, per point, a scalar Newton iteration on the
 * out-of-plane strain \f$ E_{33} \f$, and then statically condenses the 3D
 * material tangent onto the in-plane (2D) tangent. It is the constitutive core
 * consumed by the legacy adapter (task 14) and the Phase-4 shell provider.
 *
 * It works PURELY on hand-filled \ref gsMaterialData: no geometry, no
 * gsFunctionSet. The caller is responsible for ALL curvilinear-to-Cartesian
 * transforms; inputs and outputs live in one local Cartesian frame.
 *
 * ### I/O layouts and conventions
 *  - Input  \c e2D : \f$3\times N\f$, shell Voigt-3 order
 *    \f$[E_{11},\,E_{22},\,2E_{12}]\f$ (ENGINEERING shear on the off-diagonal).
 *  - Input  \c params : per-parameter \f$1\times N\f$ rows, already evaluated at
 *    the N points (e.g. via \ref gsMaterialBase::precomputeParameters).
 *  - Output \c S2D : \f$3\times N\f$, shell Voigt-3 \f$[S_{11},\,S_{22},\,S_{12}]\f$.
 *  - Output \c C2D : \f$9\times N\f$, per point a symmetric \f$3\times3\f$ block
 *    via \c reshapeCol(k,3,3), ENGINEERING convention, ordered \f$[11,22,12]\f$
 *    to match gsMaterialMatrixNonlinear's condensed C
 *    (see gsMaterialMatrixNonlinear.hpp:2097-2107).
 *  - Output \c E33 : \f$1\times N\f$, the converged out-of-plane strain.
 *  - Optional \c dE33dE : \f$3\times N\f$, \f$ \partial E_{33}/\partial E_{2D} =
 *    -C_{3,p}/C_{3333} \f$ (implicit-function theorem on \f$S_{33}=0\f$), the
 *    chain-rule input used later by the TFT strategy.
 *  - Optional \c itersOut : the maximum, over all points, of the number of
 *    Newton *updates applied* (not law evaluations). A linear / St.Venant-
 *    Kirchhoff law is affine in \f$E_{33}\f$, so exactly ONE update lands on the
 *    solution and \c itersOut == 1.
 *
 * ### Voigt index maps (TRAP)
 * PFF full 6-Voigt order (gsPhaseFieldFractureUtils.h:21-41) is
 * \f$0{=}11,\,1{=}22,\,2{=}33,\,3{=}12,\,4{=}23,\,5{=}13\f$. The in-plane set
 * used for condensation is \c I = {0,1,3} (=11,22,12) and the condensed index is
 * 2 (=33). Do NOT conflate this with the RETURNED \f$3\times3\f$ shell-Voigt
 * matrices, whose slot 2 is the 12 shear: PFF flat "2" is the 33 component,
 * shell-Voigt slot "2" is the 12 shear -- two different "2"s (this exact
 * confusion is called out at gsMaterialMatrixNonlinear.hpp:2097-2107).
 *
 * ### Algorithm and complexity
 * Per Newton sweep there is exactly ONE batched stress eval + ONE batched
 * tangent eval over the whole point set (a per-point active mask skips converged
 * points but the batch is never repacked -- masked refill is cheaper than
 * repacking). Total cost = (#sweeps) x (one batched stress + one batched tangent
 * eval), i.e. what the per-output legacy path (_eval3D_Compressible_C33) pays
 * 4-6 separate times.
 *
 * ### Legacy-constant provenance
 * Mirrors gsMaterialMatrixNonlinear::_eval3D_Compressible_C33
 * (gsMaterialMatrixNonlinear.hpp:1865-1935): Newton on \f$c_{33}=C(2,2)\f$ with
 * initial guess \f$c_{33}=1/J_0^2\f$, update \f$dc_{33}=-2S_{33}/C_{3333}\f$,
 * \c tol 1e-10 on \f$|dc_{33}|\f$, \c itmax 100. In the \f$E_{33}\f$ variable
 * (\f$C_{33}=2E_{33}+1\f$) this is exactly \f$dE_{33}=-S_{33}/C_{3333}\f$ with
 * initial guess \f$E_{33}=\tfrac12(1/J_0^2-1)\f$, where
 * \f$1/J_0^2 = 1/\det(2E_{2D}+I)\f$ when \f$E_{13}=E_{23}=0\f$ (Kirchhoff-Love).
 * NOTE the factor 2: this class' default \c tol=1e-10 is on \f$|dE_{33}|\f$,
 * whereas legacy's \c 1e-10 on \f$|dc_{33}|=|2\,dE_{33}|\f$ corresponds to
 * \f$|dE_{33}|<5\times10^{-11}\f$ -- a factor-2 looser criterion here, ample for
 * the ~1e-8 adapter-vs-legacy parity target (task 15).
 *
 * The law is queried on BOTH \c data.strain (read e.g. by gsLinearMaterial) and
 * \c data.deformationGradient (read e.g. by the hyperelastics), filled
 * consistently so that \f$E=\tfrac12(F^\top F-I)\f$.
 *
 * @tparam T Real type
 * @ingroup KLShell
 */
template <class T>
class gsPlaneStressCondensation
{
public:

    /**
     * @param[in] law    Non-owning pointer to the 3D PFF material law (the
     *                   caller owns the law).
     * @param[in] tol    Newton tolerance on \f$|dE_{33}|\f$ (default 1e-10).
     * @param[in] itmax  Maximum number of Newton updates per point (default 100).
     */
    gsPlaneStressCondensation(const gsMaterialBase<T> * law,
                              T tol = 1e-10, index_t itmax = 100)
    :
    m_law(law), m_tol(tol), m_itmax(itmax)
    {
        GISMO_ASSERT(m_law != nullptr, "gsPlaneStressCondensation: null material law.");
    }

    /**
     * @brief   Condenses the 3D law to a plane-stress 2D response over N points.
     *
     * @param[in]  e2D      \f$3\times N\f$ in-plane strain (shell Voigt-3,
     *                      engineering shear), local Cartesian frame.
     * @param[in]  params   Per-parameter \f$1\times N\f$ rows at the N points.
     * @param[out] S2D      \f$3\times N\f$ condensed stress \f$[S_{11},S_{22},S_{12}]\f$.
     * @param[out] C2D      \f$9\times N\f$ condensed tangent, \f$3\times3\f$ per
     *                      point (reshapeCol), ordered \f$[11,22,12]\f$.
     * @param[out] E33      \f$1\times N\f$ converged out-of-plane strain.
     * @param[out] dE33dE   Optional \f$3\times N\f$, \f$\partial E_{33}/\partial E_{2D}\f$.
     * @param[out] itersOut Optional; max Newton updates applied over the batch.
     * @param[in]  wantC2D  If false, the condensed tangent is NOT produced: the
     *                      FINAL-state \c compute_matrix_into and the per-point
     *                      Schur condensation are skipped and \a C2D is left
     *                      untouched (not even resized). Requires
     *                      <tt>dE33dE == nullptr</tt> (the IFT derivative needs
     *                      the final tangent). Default \c true, i.e. the historic
     *                      behavior -- the legacy adapter is unaffected.
     *
     * @warning \a wantC2D gates ONLY the final tangent. The Newton loop itself
     *          needs BOTH \c compute_stress_into (for \f$S_{33}\f$) and
     *          \c compute_matrix_into (for \f$C_{3333}\f$) on EVERY sweep; those
     *          two calls are structurally NOT gateable here. Do not "optimize"
     *          them away: without \f$C_{3333}\f$ there is no Newton update.
     */
    void condense(const gsMatrix<T> & e2D,
                  const std::vector<gsMatrix<T>> & params,
                        gsMatrix<T> & S2D,
                        gsMatrix<T> & C2D,
                        gsMatrix<T> & E33,
                        gsMatrix<T> * dE33dE = nullptr,
                        index_t     * itersOut = nullptr,
                        bool          wantC2D = true) const
    {
        const index_t N = e2D.cols();
        GISMO_ASSERT(e2D.rows()==3, "e2D must be 3 x N (shell Voigt-3), got "<<e2D.rows()<<" rows.");
        GISMO_ASSERT((index_t)params.size()==(index_t)m_law->numParameters(),
                     "params has "<<params.size()<<" entries, law needs "<<m_law->numParameters()<<".");
        GISMO_ASSERT(wantC2D || dE33dE==nullptr,
                     "gsPlaneStressCondensation: dE33dE = -C(33,p)/C(33,33) is built from the FINAL "
                     "condensed tangent, which wantC2D=false does not compute.");

        // Hand-filled material data: pure 3D constitutive query, no geometry.
        gsMaterialData<T> data;
        data.dim   = 3;
        data.size  = N;
        data.patch = 0;
        data.parameters = params;               // per-point rows, already evaluated
        data.strain.resize(9,N);                // (dim*dim) x N full tensor, col-major
        data.deformationGradient.resize(9,N);   // (dim*dim) x N full tensor, col-major

        E33.resize(1,N);

        // ---- Initial guess: E33 = 1/2 (1/det(2E_2D+I) - 1) and first fill ----
        for (index_t k=0; k!=N; ++k)
        {
            // In-plane 2x2 block of C = 2E+I: [[a, b],[b, d]], with 2E12 = e2D(2,k).
            const T a = 2*e2D(0,k) + 1;
            const T d = 2*e2D(1,k) + 1;
            const T b = e2D(2,k);
            const T det = a*d - b*b;
            GISMO_ENSURE(det > 0, "gsPlaneStressCondensation: degenerate in-plane state at point "
                         <<k<<": det(2E_2D+I) = "<<det<<" <= 0.");
            E33(0,k) = 0.5*(1.0/det - 1.0);
            _fillPoint(data, k, e2D, E33(0,k));
        }

        // ---- Batched Newton on S33=0 ----
        // One stress + one tangent eval per sweep over the FULL batch; a per-point
        // active mask skips converged points without repacking the batch.
        gsMatrix<T> Sfull, Cvoigt;
        std::vector<char> active(N, 1);
        gsVector<index_t> iters = gsVector<index_t>::Zero(N); // updates applied per point
        index_t nactive = N;
        T worst = 0;
        for (index_t sweep=0; nactive>0 && sweep<=m_itmax; ++sweep)
        {
            m_law->compute_stress_into(data, Sfull);   // 9  x N full stress tensor
            m_law->compute_matrix_into(data, Cvoigt);  // 36 x N Voigt tangent
            worst = 0;
            for (index_t k=0; k!=N; ++k)
            {
                if (!active[k]) continue;
                // full tensor (2,2) -> col-major flat 2*3+2 = 8
                const T S33   = Sfull (8, k);
                // Voigt (2,2) = C_3333 -> col-major flat 2*6+2 = 14
                const T C3333 = Cvoigt(14, k);
                const T dE33  = -S33 / C3333;
                const T adE   = math::abs(dE33);
                if (adE > worst) worst = adE;
                if (adE < m_tol)
                {
                    active[k] = 0;
                    --nactive;
                }
                else
                {
                    E33(0,k) += dE33;
                    _fillPoint(data, k, e2D, E33(0,k));
                    ++iters[k];
                }
            }
        }
        GISMO_ENSURE(nactive==0, "gsPlaneStressCondensation: Newton did not converge in "
                     <<m_itmax<<" iterations; worst |dE33| = "<<worst<<" (tol = "<<m_tol<<").");
        if (itersOut) *itersOut = (N>0) ? iters.maxCoeff() : 0;

        // ---- Final sweep at the converged state, then extract + condense ----
        // The stress is always needed (S2D); the tangent only when it is asked for.
        m_law->compute_stress_into(data, Sfull);
        if (wantC2D)
            m_law->compute_matrix_into(data, Cvoigt);

        S2D.resize(3,N);
        if (wantC2D) C2D.resize(9,N);
        if (dE33dE) dE33dE->resize(3,N);

        // PFF Voigt index set for in-plane [11,22,12] and the condensed 33 index.
        const index_t I[3] = {0,1,3};
        const index_t i33  = 2;
        for (index_t k=0; k!=N; ++k)
        {
            const gsAsMatrix<T,Dynamic,Dynamic> S = Sfull.reshapeCol(k,3,3);   // full tensor

            // S2D = [S11, S22, S12] (S33 ~ 0 by construction).
            S2D(0,k) = S(0,0);
            S2D(1,k) = S(1,1);
            S2D(2,k) = S(0,1);

            if (!wantC2D) continue;

            const gsAsMatrix<T,Dynamic,Dynamic> C = Cvoigt.reshapeCol(k,6,6);  // Voigt 6x6

            // Static condensation of the 33 row/column (scalar-C3333 form, same
            // structure as legacy _Cijkl, gsMaterialMatrixNonlinear.hpp:1454-1460).
            const T inv33 = 1.0 / C(i33,i33);
            gsAsMatrix<T,Dynamic,Dynamic> Ccond = C2D.reshapeCol(k,3,3);
            for (index_t p=0; p!=3; ++p)
                for (index_t q=0; q!=3; ++q)
                    Ccond(p,q) = C(I[p],I[q]) - C(I[p],i33) * C(i33,I[q]) * inv33;

            if (dE33dE)
            {
                // IFT on S33=0: dE33/dE_p = -C(33,p)/C(33,33), the TFT chain-rule input.
                (*dE33dE)(0,k) = -C(i33,I[0]) * inv33;
                (*dE33dE)(1,k) = -C(i33,I[1]) * inv33;
                (*dE33dE)(2,k) = -C(i33,I[2]) * inv33;
            }
        }
    }

    // Reserved slot: an analytic incompressible plane-stress strategy
    // (IncompressibleAnalytic) will complement this compressible Newton path.
    // No interface scaffolding is provided here yet.

protected:

    /**
     * @brief   Fills the kinematic fields of point \a k from the in-plane strain
     *          and the current \a E33, consistently on BOTH strain and F.
     *
     * Builds the symmetric 3x3 strain \f$E\f$ (engineering shear halved,
     * \f$E_{13}=E_{23}=0\f$) and the deformation gradient \f$F=L^\top\f$ from the
     * Cholesky factor of \f$C=2E+I=L L^\top\f$ (SPD 3x3, fixed-size, no heap), so
     * \f$F^\top F=C\f$ and \f$\det F=\prod\mathrm{diag}(L)>0\f$ -- all any PFF law
     * consumes (they never use F's rotation part).
     */
    void _fillPoint(gsMaterialData<T> & data, const index_t k,
                    const gsMatrix<T> & e2D, const T E33val) const
    {
        gsAsMatrix<T,Dynamic,Dynamic> E = data.strain.reshapeCol(k,3,3);
        E.setZero();
        E(0,0) = e2D(0,k);
        E(1,1) = e2D(1,k);
        E(2,2) = E33val;
        E(0,1) = E(1,0) = 0.5*e2D(2,k);   // engineering shear halved; E13=E23=0

        // C = 2E + I is SPD; F = L^T with C = L L^T (Cholesky, fixed 3x3).
        // Native gsEigen fixed-size type: the gsMatrix<T,3,3> wrapper has no
        // Eigen evaluator specialization and does not instantiate inside LLT.
        gsEigen::Matrix<T,3,3> Cmat = 2.0 * E;
        Cmat(0,0) += 1.0; Cmat(1,1) += 1.0; Cmat(2,2) += 1.0;
        gsEigen::LLT<gsEigen::Matrix<T,3,3> > llt(Cmat);
        GISMO_ENSURE(llt.info()==gsEigen::Success,
                     "gsPlaneStressCondensation: C = 2E+I not SPD at point "<<k
                     <<" (state outside the physical range).");
        const gsEigen::Matrix<T,3,3> L = llt.matrixL();
        data.deformationGradient.reshapeCol(k,3,3) = L.transpose();
    }

    const gsMaterialBase<T> * m_law; ///< Non-owning pointer to the 3D law.
    T       m_tol;                   ///< Newton tolerance on |dE33|.
    index_t m_itmax;                 ///< Maximum Newton updates per point.

}; // class gsPlaneStressCondensation

} // namespace gismo

#endif // gsPhaseFieldFracture_ENABLED
