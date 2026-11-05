/** @file gsMaterialMatrixTFTSplit.h

    @brief Tension-field material decorator based on the Zhang-Kiendl split (linearised or exact wrinkled stress)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
*/

#pragma once

#include <functional>
#include <vector>

#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.h>

namespace gismo
{

namespace detail_tftsplit
{

/// Options of the pointwise core (read once per evaluation call from the class's gsOptionList)
template <class T>
struct WrinkleOptions
{
    index_t tensionField = 2;   ///< 0: strain-based, 1: stress-based, 2: mixed (default)
    T       eta = 0;            ///< Zhang (2026) Eq. 40 regularisation
    T       slackMultiplier = 0;///< slack stress factor; if <= 0, eta is used instead
    bool    exact = false;      ///< exact (non-truncated) wrinkled state; only honoured for linear=false
};

/// Result of the pointwise core for one point/fibre.
/// Principal values are ordered ascending: index "min" is the minor value.
template <class T>
struct WrinkleResult
{
    index_t state;          ///< criterion state: 1 taut, 0 wrinkled, -1 slack (S2crit NOT applied)
    T Smin, Smax;           ///< principal values of the UNMODIFIED S w.r.t. the reference metric
    T Emin, Emax;           ///< principal values of E w.r.t. the reference metric
    T theta;                ///< n = (cos θ, sin θ) covector components, θ ∈ [0,π); NaN unless wrinkled
    T gamma;                ///< γ with |n|_Euclid = 1 (TFT normalisation); NaN unless wrinkled
    T gammaHat;             ///< γ̂ = γ · G^{αβ} n_α n_β (parametrisation-independent); NaN unless wrinkled
    T residual;             ///< |f(θ)| at the θ used (linearised), max(|n1.S'|,|n2.S'|) at the converged iterate (exact); NaN unless wrinkled
    bool guardFallback;     ///< closed-form θ rejected by the isotropy guard
    bool tftFailure;        ///< guard failed AND TFT's eval_theta also failed its residual (closed form kept if usable, else TFT's root), or n1ᵀCn1 <= 0
    bool exactFallback;     ///< exact solve attempted and not converged / not a minimum: linearised result returned
    bool exactReclassified; ///< exact solve converged to gamma <= 0 (taut) or zero/negative wrinkle tension (slack)
    index_t exactIterations;///< number of S/C evaluations at E' made by the exact solve (0 if not attempted)
    gsMatrix<T> S;          ///< 3x1 S_eff
    gsMatrix<T> C;          ///< 3x3 C_eff (A_ij = dS_eff_i/dE_j); size 0x0 when not requested
    gsMatrix<T> E;          ///< 3x1 elastic strain (engineering Voigt), see eval3D_strain
    gsMatrix<T> Gcov;       ///< 2x2 reference metric the point was evaluated with
    gsMatrix<T> gdef;       ///< 2x2 covariant deformed metric (geometry path only; 0x0 otherwise)
    gsMatrix<T> Gbasis;     ///< 3x2 reference covariant basis vectors (columns, Cartesian; geometry path only, 0x0 otherwise)
    gsMatrix<T> gbasis;     ///< 3x2 deformed covariant basis vectors (geometry path only, 0x0 otherwise)
};

} // namespace detail_tftsplit

/**
 * @brief      Tension-field-theory material decorator following the split of Zhang and Kiendl
 *
 * Wraps any \ref gsMaterialMatrixBase and replaces its stress and tangent by those of a wrinkled
 * membrane (linearised or exact mode).
 *
 * <b>References</b>
 * - Zhang & Kiendl, CMAME 459 (2026) 119083: stress split S = S^+ + S^- with the mixed criterion
 *   (Eq. 27 spectral decomposition, Eq. 28-30 wrinkling magnitude, Eq. 34 wrinkled stress, Eq. 35 split
 *   and criterion, Eq. 36-38 tangent, Eq. 40 regularisation).
 * - Zhang & Kiendl, CMAME 432 (2024) 117386: strain-spectral split, wrinkling intensity |E_2|.
 * - Zhang & Kiendl, J. Mech. Phys. Solids 205 (2025) 106331: strain and stress split, mixed criterion, SvK.
 * - Verhelst, Moller & Den Besten, "A wrinkling model for general hyperelastic materials based on
 *   tension field theory", CMAME 441 (2025) 117955: tension-field formulation of
 *   \ref gsMaterialMatrixTFT (Eq. 34 criteria, Eq. 72-73 linearised wrinkled stress).
 *
 * <b>Conventions</b>
 * - Stress Voigt S = (S^11, S^22, S^12) (contravariant), engineering strain E = (E_11, E_22, 2E_12)
 *   (covariant), dS = C dE, and the tangent is A_ij = dS_i/dE_j (stored column-major, 9 entries).
 * - Point state: 1 taut, 0 wrinkled, -1 slack. The slack boundary is strict (E_max < 0), as in
 *   \ref gsMaterialMatrixBaseDim; the literature uses E_max <= 0.
 * - Principal values are ordered ascending (min = minor). They are the eigenvalues of the pencil
 *   S^{ab} n_b = lambda G^{ab} n_b, G being the reference metric (see below). Index 0 (minor) is
 *   S_p,1 of Verhelst et al. (Eq. 34) and S_2 of Zhang & Kiendl, who order S_1 >= S_2 (Eq. 27).
 * - The wrinkle direction n = (cos theta, sin theta) holds covector components, theta in [0,pi).
 *   gamma is the wrinkle strain amplitude for |n|_Euclid = 1 (the \ref gsMaterialMatrixTFT
 *   normalisation); gammaHat = gamma G^{ab} n_a n_b is independent of the parametrisation and equals
 *   the wrinkling magnitude mu_w of Zhang & Kiendl (2026, Eq. 28-30): their M_2 = N_2 (x) N_2 has N_2
 *   a unit vector in G, so n = |n|_G N_2 and gamma n (x) n = gamma |n|_G^2 M_2.
 * - Direction vectors. With n = (c,s), the covector m = (-s,c), and engineering Voigt form,
 *   n_1 = n (x) n = (c^2, s^2, 2cs), n_2 = sym(n (x) m) = (-sc, sc, c^2-s^2),
 *   n_3 = dn_2/dtheta = (s^2-c^2, c^2-s^2, -4cs), n_4 = n_3 + n_1 = m (x) m = (s^2, c^2, -2cs).
 *   Hence n_1.S = S^{ab} n_a n_b, n_2.S = S^{ab} n_a m_b and n_4.S = S^{ab} m_a m_b. The direction
 *   along the wrinkles is the vector t^a = (-sin theta, cos theta), which differs from the covector m
 *   on a skewed metric. At a converged exact state n is a null eigen-covector of S', so
 *   S'^{ab} = sigma t^a t^b with sigma = n_4.S', which is positive iff the stress along the wrinkles
 *   is tensile.
 * - Reference metric: the fibre metric G(z) of the undeformed shell for linear=false, the
 *   mid-surface metric A for linear=true (the wrapped \ref gsMaterialMatrixLinear is built on A).
 *
 * <b>Wrinkled state</b>
 * With n_1 and n_2 as in the conventions above, the wrinkled stress is the linearised
 * (Taylor-truncated) one,
 * @f[ S' = S + \gamma\, C n_1,\qquad \gamma = -\frac{n_1\cdot S}{n_1^T C n_1}, @f]
 * (Verhelst et al. Eq. 73; Zhang & Kiendl 2026, Eq. 30 and 34)
 * with theta the root of n_2 . S' = 0, found in closed form as the minor solution of the metric
 * pencil (Eq. 27 of Zhang & Kiendl, generalised to the pencil), validated by an isotropy guard
 * (|f| tolerance) and replaced by the root of
 * \ref gsMaterialMatrixTFT::eval_theta if the guard rejects it (counted, see numThetaFallbacks()).
 * The tangent J = C + C D + gamma H is derived in gsMaterialMatrixTFTSplit.hpp.
 * The regularised result is S_eff = (1-eta) S' + eta S, C_eff = (1-eta) J + eta C, i.e.
 * S_eff = S^+ + eta S^- with S^+ = S' (the TFT-modified stress) and S^- = S - S^+
 * (Zhang & Kiendl 2026, Eq. 35 and 40); slack points return SlackMultiplier*(S,C)
 * (Eta if SlackMultiplier <= 0).
 *
 * <b>Accuracy of the linearised mode</b>
 * Let H be the 3x3 matrix whose column d is (dC/dE_d) n_1. Then
 * @f[ S(E + \gamma n_1) = S + \gamma C n_1 + \tfrac12 \gamma^2 H n_1 + O(\gamma^3). @f]
 * The truncation error of the linearised S' is therefore \f$ \tfrac12\gamma^2 H n_1 \f$: O(gamma^2)
 * absolute, O(gamma) relative to the wrinkling correction gamma C n_1. The linearised gamma carries
 * an O(gamma) relative error as well. Whether this matters depends on the size of gammaHat in the
 * application, not on the order. The linearised mode is exact for a wrapped law that is linear in E
 * (SvK).
 * - For a hyperelastic wrapped law the linearised S' need not stay tensile along the wrinkles:
 *   n_4.S' <= 0 can occur, for example for incompressible neo-Hookean at unit stretch along the
 *   wrinkles. The linearised mode has no n_4.S' check and keeps such points wrinkled, so the stress
 *   jumps at the wrinkled/slack boundary, where the slack value is SlackMultiplier*S or Eta*S.
 *   Energy continuity across state switches holds for SvK and for the exact mode, not for the
 *   linearised mode with a hyperelastic law.
 * - In every mode the tangent jumps at the taut/wrinkled boundary (the relaxed energy W_r is C^1,
 *   not C^2). This is intrinsic to tension-field theory.
 *
 * <b>Exact mode</b>
 * With E' = E + gamma n_1, the pair (gamma, theta) solves n_1 . S(E') = n_2 . S(E') = 0 by a local
 * Newton method started from the linearised (gamma, theta). The iteration is scalar in gamma when the
 * closed-form theta passes the isotropy guard (coaxiality) and 2x2 otherwise. With
 * W(gamma,theta) = W(E + gamma n_1(theta)) one has dW/dgamma = n_1.S(E') and
 * dW/dtheta = 2 gamma n_2.S(E') (dn_1/dtheta = 2 n_2), so the solution of R = 0 is a stationary point
 * of W(gamma,theta) for gamma != 0, i.e. of Pipkin's relaxed energy W_r(E) = min W(E + gamma n (x) n).
 * The stress is S' = S(E') and the tangent is
 * C' - [C'n_1, C'n_2] diag(1,2 gamma) Jac_0^{-1} [n_1^T C'; n_2^T C'], with C' = C(E') and
 * Jac_0 = [[n_1.C'n_1, 2 gamma n_1.C'n_2], [n_2.C'n_1, n_4.S' + 2 gamma n_2.C'n_2]],
 * the R = 0 form of the Jacobian of R. The tangent is symmetric whenever C(E') is, and is the
 * consistent tangent up to O(||R||) at the converged iterate.
 *
 * At R = 0 the Hessian of W(gamma,theta) is diag(1, 2 gamma) Jac_0. For gamma > 0 the stationary
 * point is a minimum iff n_1.C'n_1 > 0 and det Jac_0 > 0; this is guaranteed when n_4.S' > 0 and C' is
 * positive definite on span{n_1, n_2}. In the major principal direction stationarity gives gamma < 0
 * (to first order gamma = -S_max/(n_1.C n_1)), which is not an admissible wrinkle; the iteration
 * therefore starts from the minor direction. The linearised gamma is the first scalar Newton step in
 * gamma from gamma = 0 at fixed theta (not so for the 2x2 iteration at a non-coaxial theta).
 *
 * Outcomes of a wrinkled point in exact mode, with the counters that record them:
 * - Not attempted if the direction search failed or n_1.C n_1 <= 0: counted in numThetaFallbacks()
 *   and numTFTFailures() (the closed form was already rejected by the guard or unusable, e.g. for a
 *   degenerate pencil), not in the exact-mode counters.
 *   A failed search with n_1.C n_1 > 0 gets the linearised result with the kept direction; a point with
 *   n_1.C n_1 <= 0 gets no wrinkle correction (S_eff = S, C_eff = C, state 0).
 * - A missing S/C provider where the solve is attempted is an error (GISMO_ERROR).
 * - Converged with gamma <= 0: taut. S, C and E are returned unmodified (eta is not applied).
 *   Counted in numExactReclassified(); the state becomes 1 and theta, gamma, gammaHat are NaN.
 * - Converged with n_4.S' <= 0 or NaN (no tension along the wrinkles): slack, with result
 *   sigma*(S,C,E), sigma = SlackMultiplier if > 0 else Eta. Counted in numExactReclassified(); the
 *   state becomes -1 and theta, gamma, gammaHat are NaN.
 * - Converged at a stationary point that is not a minimum (n_1.C'n_1 <= 0 or det Jac_0 <= 0), or not
 *   converged: the linearised result is returned and counted in numExactFallbacks(). Non-convergence
 *   covers non-finite S' or C', non-finite gamma or theta after a step, n_1.C'n_1 <= 0 in the
 *   scalar iteration, a singular 2x2 Jacobian, a residual that does not decrease over two iterations
 *   after the third, and no convergence within 20 iterations per mode.
 * - Otherwise accepted: S_eff = (1-eta) S(E') + eta S, E = E + gamma n_1,
 *   C_eff = (1-eta) J_ex + eta C.
 *
 * <b>Options</b>
 * - TensionField (int, 2): 0 strain-based, 1 stress-based, 2 mixed criterion. The decorator evaluates
 *   the criterion itself, with its own option (not the wrapped material's). 2 is the mixed criterion
 *   of Zhang & Kiendl (2026, Eq. 35) and Verhelst et al. (Eq. 34).
 * - Eta (real, 0): regularisation (Zhang & Kiendl 2026, Eq. 40)
 * - SlackMultiplier (real, 0): slack stress SlackMultiplier*S if > 0, else Eta*S; it does not
 *   affect wrinkled points
 * - S2crit (real, 0): display-only (S_2,crit of Zhang & Kiendl 2026); a wrinkled point is reported
 *   wrinkled by eval3D_tensionfield iff its unmodified S_min < S2crit (S and C are not affected)
 * - Midsurface (switch, false): evaluate once on the mid-surface (MatIntegration::Constant).
 *   For membrane assemblies only: the bending moment output VectorM is then undefined and raises an error.
 * - Exact (switch, false): exact (non-truncated) wrinkled stress, see above; for linear=true the exact
 *   and the linearised stress coincide and the option has no effect
 *
 * Cauchy outputs are S_eff/J, with J the stretch from the wrapped eval3D_detF, i.e. the components
 * of sigma = F S_eff F^T / J in the deformed basis. This holds in exact mode too: J is that of the
 * actual deformation, not of E' (the thickness stretch at E' is not computed).
 *
 * <b>Counters</b>
 * numThetaFallbacks(), numTFTFailures(), numExactFallbacks() and numExactReclassified() are
 * cumulative from construction until resetCounters() and are incremented atomically. Every output
 * call (eval3D_matrix, eval3D_vector, eval3D_pstress, eval3D_theta, ...) re-runs the pointwise core,
 * as does every public wrinkle() call, so a point evaluated by both eval3D_matrix and eval3D_vector
 * is counted twice. Each matching gsWarn message is printed at most once per program run and
 * template instantiation (dim, T, linear) (a function-local static flag).
 *
 * <b>Post-processing outputs</b> (columns j*u.cols()+k as for every multi-point output). All principal
 * quantities are ordered ascending: index/column 0 is the minor (minimum) value, which is S_p,1 in
 * Verhelst et al. (2025) and sigma_2 (E_2) in Zhang & Kiendl, who number the major value 1.
 * - eval3D_pstressDir: 9 rows, the column-major 3x3 matrix [d_min | d_max | N] of Cartesian unit vectors,
 *   the principal frame of S_eff in the tangent plane of the reference basis (geometry path); any state
 * - eval3D_wrinkleDir: 3 rows, unit vector normalise(-sin(theta) g_1 + cos(theta) g_2) along the wrinkles
 *   on the deformed surface (zero unless the criterion state is wrinkled)
 * - eval3D_wrinkleIntensity: max(0,-E_min), E_min the minor principal value of the total strain w.r.t. the
 *   reference metric (zero unless the criterion state is wrinkled)
 * - eval3D_gammaHat: gamma G^{ab} n_a n_b (NaN unless wrinkled)
 * The state is the criterion state after any exact-mode reclassification; S2crit is not applied.
 *
 * @tparam     dim     The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T       Real type
 * @tparam     linear  The wrapped material is linear in E (S = C E): the tangent has no dC/dE term
 *
 * @ingroup    KLShell
 */
template < short_t dim, class T, bool linear=false >
class gsMaterialMatrixTFTSplit : public gsMaterialMatrixBaseDim<dim,T>
{
public:
    /// Shared pointer for gsMaterialMatrixTFTSplit
    typedef memory::shared_ptr< gsMaterialMatrixTFTSplit > Ptr;

    /// Unique pointer for gsMaterialMatrixTFTSplit
    typedef memory::unique_ptr< gsMaterialMatrixTFTSplit > uPtr;

    /// Shared pointer to a \ref gsMaterialMatrixBase
    typedef typename gsMaterialMatrixBase<T>::Ptr material_ptr;

    /// Base class
    using Base = gsMaterialMatrixBaseDim<dim,T>;

    /// Shared pointer to a \ref gsFunctionSet
    typedef typename gsFunctionSet<T>::Ptr function_ptr;

    /// Provider of dC/dE_d (d=0,1,2), each 3x3 in engineering Voigt; called only for wrinkled points
    typedef std::function<void(std::vector<gsMatrix<T>> &)> dC_provider;

    /// Provider of the wrapped material's stress and tangent at a prescribed elastic strain:
    /// Eel (3x1, engineering Voigt) -> S (3x1), C (3x3, or 9x1 column-major), A_ij = dS_i/dE_j.
    /// Used by the exact wrinkling solve only.
    typedef std::function<void(const gsMatrix<T> & Eel, gsMatrix<T> & S, gsMatrix<T> & C)> SC_provider;

    // Define clone functions
    GISMO_CLONE_FUNCTION(gsMaterialMatrixTFTSplit)

public:

    /**
     * @brief      Constructs the decorator around a (non-owned) material matrix
     *
     * @param      materialMatrix  The wrapped material matrix
     */
    gsMaterialMatrixTFTSplit(gsMaterialMatrixBase<T> * materialMatrix)
    :
    gsMaterialMatrixTFTSplit(memory::make_shared_not_owned(materialMatrix))
    {}

    /**
     * @brief      Constructs the decorator around a clone of a material matrix
     *
     * @param      materialMatrix  The wrapped material matrix
     */
    gsMaterialMatrixTFTSplit(gsMaterialMatrixBase<T> & materialMatrix)
    :
    gsMaterialMatrixTFTSplit(memory::make_shared(materialMatrix.clone().release()))
    {}

    /**
     * @brief      Constructs the decorator around a shared material matrix
     *
     * @param      materialMatrix  The wrapped material matrix
     */
    gsMaterialMatrixTFTSplit(const material_ptr & materialMatrix);

    /// Updates the options (Exact is ignored, with a message, for linear=true)
    void setOptions(gsOptionList opt) override;

    /// Reads TensionField, Eta, SlackMultiplier and Exact from the option list
    detail_tftsplit::WrinkleOptions<T> wrinkleOptions() const;

    /// Number of wrinkled evaluations whose closed-form theta was rejected by the isotropy guard
    index_t numThetaFallbacks() const { return m_nThetaFallbacks; }

    /// Number of those where the eval_theta of \ref gsMaterialMatrixTFT also failed (the closed form is
    /// kept when usable, otherwise TFT's root), plus the cases with n1^T C n1 <= 0
    index_t numTFTFailures() const { return m_nTFTFailures; }

    /// Number of exact solves that did not converge or reached a saddle (linearised result used)
    index_t numExactFallbacks() const { return m_nExactFallbacks; }

    /// Number of exact solves that converged to a taut (gamma <= 0) or slack (zero/negative wrinkle tension) state
    index_t numExactReclassified() const { return m_nExactReclassified; }

    /// Sets all counters to zero
    void resetCounters() { m_nThetaFallbacks = 0; m_nTFTFailures = 0; m_nExactFallbacks = 0; m_nExactReclassified = 0; }

    /**
     * @brief      Pointwise core: classifies the point and returns the effective stress and tangent
     *
     * No geometry is involved. wrinkle() itself never calls the wrapped material: it is reached only
     * through the providers \a dCdE and \a SCatE. The root-finder fallback of
     * \ref gsMaterialMatrixTFT::eval_theta on guard failure uses only the S, C and E passed in.
     *
     * @param[in]  S           3x1 unmodified stress
     * @param[in]  C           3x3 unmodified tangent, or its 9x1 column-major form
     * @param[in]  E           3x1 engineering Voigt strain
     * @param[in]  Gcov        2x2 covariant reference metric
     * @param[in]  opt         The options
     * @param[in]  needTangent If false, the tangent is not computed and the result C is empty
     * @param[in]  dCdE        Provider of dC/dE_d, called at most once per call, only for a wrinkled
     *                         point with n_1.C n_1 > 0 whose linearised result is returned (also after
     *                         an exact-mode fallback) when the tangent is needed. An empty function
     *                         means dC/dE = 0. The class drivers bind eval3D_dmatrix of the wrapped
     *                         material on the geometry path, and central differences (step 1e-5 on
     *                         (g_11, g_22, 2 g_12)) of eval3D_matrix_C on the prescribed-metric path;
     *                         for linear=true the provider is empty.
     * @param[in]  SCatE       S and C of the wrapped material at E' (exact mode only), called once per
     *                         local Newton iteration. Required where the exact solve is attempted
     *                         (opt.exact, linear=false, wrinkled point with an admissible direction);
     *                         GISMO_ERROR if empty there. Both class drivers bind eval3D_vector_C and
     *                         eval3D_matrix_C at Cmat' = (G11 + 2E'_0, G22 + 2E'_1, G12 + E'_2).
     */
    detail_tftsplit::WrinkleResult<T> wrinkle(const gsMatrix<T> & S,
                                              const gsMatrix<T> & C,
                                              const gsMatrix<T> & E,
                                              const gsMatrix<T> & Gcov,
                                              const detail_tftsplit::WrinkleOptions<T> & opt,
                                              bool needTangent,
                                              const dC_provider & dCdE = dC_provider(),
                                              const SC_provider & SCatE = SC_provider()) const;

    /// See \ref gsMaterialMatrixBase for details. Constant if Midsurface, else Linear
    inline enum MatIntegration isMatIntegrated() const override
    { return m_options.getSwitch("Midsurface") ? MatIntegration::Constant : MatIntegration::Linear; }

    /// See \ref gsMaterialMatrixBase for details. Constant if Midsurface, else Linear
    inline enum MatIntegration isVecIntegrated() const override
    { return m_options.getSwitch("Midsurface") ? MatIntegration::Constant : MatIntegration::Linear; }

    /// See \ref gsMaterialMatrixBase for details
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const override
    { m_materialMat->density_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const override
    { m_materialMat->thickness_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void setParameters(const std::vector<gsFunctionSet<T>*> &pars) override
    { m_materialMat->setParameters(pars); }
    using Base::setParameters;

    /// See \ref gsMaterialMatrixBase for details
    void info() const override
    { m_materialMat->info(); }

    /// See \ref gsMaterialMatrixBase for details
    std::ostream &print(std::ostream &os) const override;

    /// See \ref gsMaterialMatrixBase for details
    void setUndeformed(function_ptr undeformed) override
    {
        Base::setUndeformed(undeformed);
        m_materialMat->setUndeformed(m_patches);
    }

    /// See \ref gsMaterialMatrixBase for details
    void setDeformed(function_ptr deformed) override
    {
        Base::setDeformed(deformed);
        m_materialMat->setDeformed(m_defpatches);
    }

    /// See \ref gsMaterialMatrixBase for details
    void setUndeformed(const gsFunctionSet<T> * undeformed) override
    {
        Base::setUndeformed(undeformed);
        m_materialMat->setUndeformed(m_patches);
    }

    /// See \ref gsMaterialMatrixBase for details
    void setDeformed(const gsFunctionSet<T> * deformed) override
    {
        Base::setDeformed(deformed);
        m_materialMat->setDeformed(m_defpatches);
    }

    /// See \ref gsMaterialMatrixBase for details
    const gsMaterialMatrixBase<T> * material() const override { return m_materialMat.get(); }

    /// See \ref gsMaterialMatrixBase for details
    gsMaterialMatrixBase<T> * material() override { return m_materialMat.get(); }

    /// See \ref gsMaterialMatrixBase for details (forwarded to the wrapped material)
    gsMatrix<T> eval3D_pstretch(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override;

    /// See \ref gsMaterialMatrixBase for details (forwarded to the wrapped material)
    gsMatrix<T> eval3D_pstretchDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override;

    /// See \ref gsMaterialMatrixBase for details (forwarded to the wrapped material)
    gsMatrix<T> eval3D_detF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// Effective tangent C_eff, 9 x (u.cols()*z.rows()), column-major with C_ij = dS_i/dE_j; identical for every \a out
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_matrix;

    /// Effective stress S_eff, 3 x (u.cols()*z.rows()). Midsurface=true and VectorM is an error
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_vector;

    /// Effective tangent (9x1) at the prescribed metric \a Cmat = (g_11,g_22,g_12) with the tensor g_12
    gsMatrix<T> eval3D_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T>& u, const T z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// Effective stress (3x1) at the prescribed metric \a Cmat = (g_11,g_22,g_12) with the tensor g_12
    gsMatrix<T> eval3D_vector_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T>& u, const T z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// Effective stress S_eff (3 rows)
    gsMatrix<T> eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;
    using Base::eval3D_stress;

    /// Elastic strain: E (taut), E + gamma n1 (wrinkled), SlackMultiplier*E or Eta*E (slack)
    gsMatrix<T> eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const override;
    using Base::eval3D_strain;

    /// Ascending principal values (2 rows) of S_eff with respect to the reference metric
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// Rows (min, max, 0) of the principal values of the elastic strain with respect to the reference metric
    gsMatrix<T> eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const override;

    /// Cauchy stress components in the deformed basis, S_eff / J
    gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;
    using Base::eval3D_CauchyStress;

    /// Same as \ref eval3D_CauchyStress
    gsMatrix<T> eval3D_CauchyVector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// Ascending eigenvalues (2 rows) of (1/J) g_cov S_eff^con, g_cov the deformed metric
    gsMatrix<T> eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// Displayed state per point: 1 taut, 0 wrinkled, -1 slack; a wrinkled point is displayed taut unless S_min < S2crit
    gsMatrix<T> eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_tensionfield;

    /// Wrinkle direction theta in [0,pi); NaN if the point is not wrinkled
    gsMatrix<T> eval3D_theta(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_theta;

    /// Wrinkle amplitude gamma (|n|_Euclid = 1); NaN if the point is not wrinkled
    gsMatrix<T> eval3D_gamma(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_gamma;

    /// Unit vector along the wrinkles on the deformed surface (3 rows); (0,0,0) unless the point is wrinkled
    gsMatrix<T> eval3D_wrinkleDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_wrinkleDir;

    /// Wrinkling intensity max(0,-E_min) of the total strain (1 row); 0 unless the point is wrinkled
    gsMatrix<T> eval3D_wrinkleIntensity(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_wrinkleIntensity;

    /// Parametrisation-independent wrinkle strain gamma G^{ab} n_a n_b (1 row); NaN if the point is not wrinkled
    gsMatrix<T> eval3D_gammaHat(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_gammaHat;

    /// Principal stress frame [d_min | d_max | N] of S_eff w.r.t. the reference basis (9 rows, column-major 3x3, ascending)
    gsMatrix<T> eval3D_pstressDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_pstressDir;

protected:
    /// Geometry-path driver: one result per column colIdx = j*u.cols()+k (point k, fibre j)
    std::vector<detail_tftsplit::WrinkleResult<T>> _evalCore(const index_t patch, const gsMatrix<T> & u,
                                                             const gsMatrix<T> & z, bool needTangent) const;

    /// Prescribed-metric driver
    detail_tftsplit::WrinkleResult<T> _evalCore_C(const gsMatrix<T> & Cmat, const index_t patch,
                                                  const gsVector<T> & u, const T z, bool needTangent) const;

    /// Exact wrinkled state: Newton on R = (n1.S(E+gamma n1), n2.S(E+gamma n1)) = 0; returns true if converged
    template <bool _linear>
    typename std::enable_if<!_linear, bool>::type
    _exactSolve(const gsVector<T,3> & E3, const T tolR, const SC_provider & SCatE, const bool scalar,
                T & gamma, T & theta, gsVector<T,3> & Sp, gsMatrix<T,3,3> & Cp, T & resNorm, index_t & nEval) const;

    /// Not available for a wrapped material that is linear in E (exact and linearised coincide); never called
    template <bool _linear>
    typename std::enable_if< _linear, bool>::type
    _exactSolve(const gsVector<T,3> & E3, const T tolR, const SC_provider & SCatE, const bool scalar,
                T & gamma, T & theta, gsVector<T,3> & Sp, gsMatrix<T,3,3> & Cp, T & resNorm, index_t & nEval) const;

    /// Stress S_eff of every result as columns
    gsMatrix<T> _collectStress(const std::vector<detail_tftsplit::WrinkleResult<T>> & res) const;

    /// Cauchy stress S_eff/J as columns
    gsMatrix<T> _collectCauchy(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const;

protected:

    mutable material_ptr m_materialMat;
    memory::shared_ptr< gsMaterialMatrixTFT<dim,T,linear> > m_tft;
    mutable index_t m_nThetaFallbacks, m_nTFTFailures, m_nExactFallbacks, m_nExactReclassified;

    using Base::m_patches;
    using Base::m_defpatches;
    using Base::m_options;

    using Base::m_data;
};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixTFTSplit.hpp)
#endif
