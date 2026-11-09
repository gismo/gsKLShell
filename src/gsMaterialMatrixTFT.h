/** @file gsMaterialMatrixTFT.h

    @brief Tension-field theory decorator of a material matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>

namespace gismo
{

namespace detail_tftsplit
{

/// Entries of a 2x2 matrix [[a,b],[c,d]]
template <class T>
struct Mat2
{
    T a, b, c, d;
};

/// Matrix G_cov * S_con of the pencil S^{ab} n_b = lambda G^{ab} n_b, S_con from the Voigt vector (S0,S1,S2)
template <class T>
inline Mat2<T> stressPencil(const gsMatrix<T> & Gcov, const T S0, const T S1, const T S2)
{
    return Mat2<T>{ Gcov(0,0)*S0 + Gcov(0,1)*S2, Gcov(0,0)*S2 + Gcov(0,1)*S1,
                    Gcov(1,0)*S0 + Gcov(1,1)*S2, Gcov(1,0)*S2 + Gcov(1,1)*S1 };
}

/// Matrix G_con * E_t of the pencil, E_t from the engineering Voigt vector (E0,E1,E2=2E_12)
template <class T>
inline Mat2<T> strainPencil(const Mat2<T> & Gcon, const T E0, const T E1, const T E2)
{
    const T h = T(0.5) * E2;
    return Mat2<T>{ Gcon.a*E0 + Gcon.b*h, Gcon.a*h + Gcon.b*E1,
                    Gcon.c*E0 + Gcon.d*h, Gcon.c*h + Gcon.d*E1 };
}

/// Eigenvalues tr/2 -+ sqrt(max(0,((a-d)/2)^2 + b c)) of a real 2x2 matrix with real spectrum.
/// The closed form avoids the determinant, which cancels for nearly equal eigenvalues.
template <class T>
inline void principalValues(const Mat2<T> & M, T & lmin, T & lmax)
{
    const T mid = T(0.5) * (M.a + M.d);
    const T dif = T(0.5) * (M.a - M.d);
    const T rad = math::sqrt(math::max(T(0), dif*dif + M.b*M.c));
    lmin = mid - rad;
    lmax = mid + rad;
}

/// Contravariant 2x2 metric from the covariant one
template <class T>
inline Mat2<T> invert2x2(const gsMatrix<T> & G)
{
    const T det = G(0,0)*G(1,1) - G(0,1)*G(1,0);
    return Mat2<T>{ G(1,1)/det, -G(0,1)/det, -G(1,0)/det, G(0,0)/det };
}

/// Voigt direction vectors n1 = (c^2,s^2,2cs), n2 = (-sc,sc,c^2-s^2), n3 = dn2/dtheta for n = (c,s)
template <class T>
struct VoigtDirs
{
    gsVector<T,3> n1, n2, n3;

    VoigtDirs(const T c, const T s)
    {
        n1 << c*c, s*s, 2*c*s;
        n2 << -s*c, s*c, c*c - s*s;
        n3 << s*s - c*c, c*c - s*s, -4*c*s;
    }
};

/// Unit vector along the wrinkles, normalise(-sin(theta) a_1 + cos(theta) a_2), from the 3x2 matrix of
/// covariant basis vectors a_alpha (columns, Cartesian). t^alpha = (-sin, cos) is annihilated by the
/// covector n = (cos, sin): t is a material line, so it is tangent to the wrinkles on the surface spanned by
/// a_alpha. The sign is arbitrary (theta and theta+pi give opposite vectors).
template <class T>
inline gsVector<T,3> wrinkleDirection(const T theta, const gsMatrix<T> & basis)
{
    gsVector<T,3> v = -math::sin(theta) * basis.col(0) + math::cos(theta) * basis.col(1);
    return v / v.norm();
}

/// gammaHat = gamma * G^{ab} n_a n_b, n = (cos theta, sin theta) covector, G^{ab} the inverse of the 2x2
/// covariant reference metric Gcov: the wrinkle strain independent of the parametrisation.
/// Invariant under theta -> theta + pi.
template <class T>
inline T gammaHat(const T gamma, const T theta, const gsMatrix<T> & Gcov)
{
    const Mat2<T> Gcon = invert2x2(Gcov);
    const T c = math::cos(theta);
    const T s = math::sin(theta);
    return gamma * (c*(Gcon.a*c + Gcon.b*s) + s*(Gcon.c*c + Gcon.d*s));
}

/// Wrinkling intensity of Zhang & Kiendl (2024), |E_2| with E_2 their (minor) second principal strain:
/// max(0, -Emin), Emin the minor principal value of the total Green-Lagrange strain w.r.t. the reference metric.
template <class T>
inline T wrinkleIntensity(const T Emin)
{
    return math::max(T(0), -Emin);
}

/// Principal frame of the in-plane stress S = (S^11,S^22,S^12) (contravariant Voigt) w.r.t. the reference
/// basis Gbasis (3x2, columns G_1, G_2, Cartesian). Writes the 3x3 matrix [d_min | d_max | N] (Cartesian
/// unit columns): N = normalise(G_1 x G_2), d_min/d_max the eigenvectors of S^{ab} G_a (x) G_b in the
/// tangent plane, in ascending eigenvalue order (column 0 = minor), d_max = N x d_min.
/// Degenerate basis (|G_1 x G_2| <= 10 eps (|G_1|^2 + |G_2|^2), non-finite G) or non-finite S: if G_1 and
/// G_2 both have an exactly zero z-component (taken as a planar problem), N = e_3 and the frame is built from
/// the longer of G_1, G_2; otherwise all nine entries are quiet NaN.
/// Complexity O(1) per point.
template <class T>
inline void principalStressFrame(const gsMatrix<T> & S, const gsMatrix<T> & Gbasis, gsMatrix<T> & frame)
{
    const gsVector<T,3> G1 = Gbasis.col(0);
    const gsVector<T,3> G2 = Gbasis.col(1);
    gsVector<T,3> N = G1.cross(G2);

    // Degenerate frame (|G_1 x G_2| <= 10 eps (|G_1|^2 + |G_2|^2), or non-finite data), same test as
    // gsMaterialMatrixBaseDim::_evalPStress. The dimension is not known here: a basis with exactly zero z-components
    // is treated as planar (dim 2, where the third row of m_gcov_ori is 0) and gets N = e_3 with the frame from the
    // longer of G_1, G_2. A dim-3 surface lying exactly in z = 0 therefore also gets this finite frame, where
    // _evalPStress returns NaN; any other degenerate basis gets NaN.
    const T nN = N.norm();
    if (!(nN > 10 * std::numeric_limits<T>::epsilon() * (G1.squaredNorm() + G2.squaredNorm())) || !S.allFinite())
    {
        const T nan = std::numeric_limits<T>::quiet_NaN();
        frame.resize(3,3);
        const T n1 = G1.norm(), n2 = G2.norm();
        const T nmax = (n1 >= n2) ? n1 : n2;
        if (S.allFinite() && G1(2) == 0 && G2(2) == 0 && G1.allFinite() && G2.allFinite() && nmax > 0)
        {
            gsVector<T,3> e3;
            e3 << 0, 0, 1;
            const gsVector<T,3> a = ((n1 >= n2) ? G1 : G2) / nmax;
            const gsVector<T,3> b = e3.cross(a);
            const T p[2][2] = { { a.dot(G1), a.dot(G2) }, { b.dot(G1), b.dot(G2) } };
            gsMatrix<T> B(2,2);
            const T Sc[2][2] = { { S(0), S(2) }, { S(2), S(1) } };
            for (index_t i = 0; i < 2; ++i)
                for (index_t j = 0; j < 2; ++j)
                {
                    T v = 0;
                    for (index_t k = 0; k < 2; ++k)
                        for (index_t l = 0; l < 2; ++l)
                            v += Sc[k][l] * p[i][k] * p[j][l];
                    B(i,j) = v;
                }
            typename gsEigen::SelfAdjointEigenSolver< typename gsMatrix<T>::Base > eig(B);
            const gsVector<T,3> dmin = eig.eigenvectors()(0,0) * a + eig.eigenvectors()(1,0) * b;
            frame.col(0) = dmin;
            frame.col(1) = e3.cross(dmin);
            frame.col(2) = e3;
        }
        else
            frame.setConstant(nan);
        return;
    }

    const gsVector<T,3> e1 = G1 / G1.norm();
    N /= nN;
    const gsVector<T,3> e2 = N.cross(e1);

    // B(i,j) = S^{kl} (e_i.G_k)(e_j.G_l): the stress tensor in the orthonormal tangent frame (e1,e2)
    const T p[2][2] = { { e1.dot(G1), e1.dot(G2) }, { e2.dot(G1), e2.dot(G2) } };
    const T Sc[2][2] = { { S(0), S(2) }, { S(2), S(1) } };
    gsMatrix<T> B(2,2);
    for (index_t i = 0; i < 2; ++i)
        for (index_t j = 0; j < 2; ++j)
        {
            T v = 0;
            for (index_t k = 0; k < 2; ++k)
                for (index_t l = 0; l < 2; ++l)
                    v += Sc[k][l] * p[i][k] * p[j][l];
            B(i,j) = v;
        }

    typename gsEigen::SelfAdjointEigenSolver< typename gsMatrix<T>::Base > eig(B);
    const gsVector<T,3> dmin = eig.eigenvectors()(0,0) * e1 + eig.eigenvectors()(1,0) * e2;
    const gsVector<T,3> dmax = N.cross(dmin);
    frame.resize(3,3);
    frame.col(0) = dmin;
    frame.col(1) = dmax;
    frame.col(2) = N;
}

} // namespace detail_tftsplit

/**
 * @brief      Tension-field theory decorator of a material model
 *
 * Principal quantities (eval3D_pstress, eval3D_CauchyPStress, eval3D_pstressDir) are ordered ascending:
 * index/column 0 is the minor (minimum) value, S_p,1 in Verhelst et al. (2025) and sigma_2 in Zhang & Kiendl, who
 * number the major value 1. They, and the Cauchy stress, derive from the modified (tension-field) stress.
 * eval3D_wrinkleDir and eval3D_wrinkleIntensity are zero unless the point is wrinkled (state 0); eval3D_gammaHat
 * is NaN unless wrinkled. Definitions are shared with gsMaterialMatrixTFTSplit (namespace detail_tftsplit).
 *
 * The wrinkled stress is the linearised S' = S(E) + gamma C(E) n_1 (Verhelst et al., CMAME 441 (2025)
 * 117955, Eq. 72-73). This is exact for a wrapped law that is linear in E and a first-order Taylor
 * truncation for hyperelastic laws. See \ref gsMaterialMatrixTFTSplit (option Exact) for the
 * non-truncated state and the error order.
 *
 * Residual stiffness (option Eta, 0 <= Eta <= 1, default 0). With S, C the stress and tangent of the wrapped law and
 * S_wr, C_wr the tension-field stress and (consistent) tangent, the points of the three states carry
 * - taut:     S, C;
 * - wrinkled: S_eff = (1-Eta) S_wr + Eta S and C_eff = (1-Eta) C_wr + Eta C;
 * - slack:    sigma S and sigma C, with sigma = SlackMultiplier if SlackMultiplier > 0, else Eta
 * (Zhang & Kiendl 2026, Eqs. 35 and 40; the same definition as option Eta of \ref gsMaterialMatrixTFTSplit). The
 * strain at slack points is scaled by sigma as well; at wrinkled points it is the elastic strain E + gamma n_1 and
 * is not blended (as in gsMaterialMatrixTFTSplit). C_eff is the derivative of S_eff, so the Newton iteration keeps
 * its rate. The state classification and theta, gamma do not depend on Eta.
 *
 * With option EtaTangentOnly, Eta enters the tangent only: C_eff and sigma C as above, while stress and strain stay
 * those of Eta = 0 (S_wr at wrinkled points, SlackMultiplier S at slack points). The converged solution is then an
 * exact tension field, since the residual is the pure tension-field stress, and Eta changes the Newton iteration
 * only. The tangent is no longer the derivative of the stress, so the convergence is at best linear.
 *
 * Complexity: S and C are quantities that the wrinkled and slack updates already hold, so the blend costs one 3-vector
 * or 3x3 axpy per point and no material evaluation. For an effective Eta of 0 no blend is evaluated at all, so the
 * pure tension-field quantities are returned bit for bit.
 *
 * @tparam     dim   The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T     Real type
 * @tparam     linear  The wrapped material is linear in E (\ref gsMaterialMatrixLinear). The tangent then
 *                     has no dC/dE term, and the reference metric is the mid-surface metric at every z.
 *
 * @ingroup    MaterialMatrix
 *
 */
template <  short_t dim,
            class T,
            bool linear=false
         >
class gsMaterialMatrixTFT : public gsMaterialMatrixBaseDim<dim,T>
{
public:
    /// Shared pointer for gsMaterialMatrixTFT
    typedef memory::shared_ptr< gsMaterialMatrixTFT > Ptr;

    /// Unique pointer for gsMaterialMatrixTFT
    typedef memory::unique_ptr< gsMaterialMatrixTFT > uPtr;

    /// Shared pointer to a \ref gsMaterialMatrixBase
    typedef typename gsMaterialMatrixBase<T>::Ptr material_ptr;

    /// Base class
    using Base = gsMaterialMatrixBaseDim<dim,T>;

    /// Shared pointer to a \ref gsFunctionSet
    typedef typename gsFunctionSet<T>::Ptr function_ptr;

    // Define clone functions
    GISMO_CLONE_FUNCTION(gsMaterialMatrixTFT)

public:

    // /**
    //  * @brief      Constructs a TFT material matrix using a reference to another material matrix
    //  *
    //  * @param      materialMatrix  The original material matrix
    //  */
    // gsMaterialMatrixTFT(typename gsMaterialMatrixBase<T>::uPtr materialMatrix)
    // :
    // gsMaterialMatrixTFT(memory::make_shared(materialMatrix->clone().release()))
    // {}

    /**
     * @brief      Constructs a TFT material matrix using a reference to another material matrix
     *
     * @param      materialMatrix  The original material matrix
     */
    gsMaterialMatrixTFT(gsMaterialMatrixBase<T> * materialMatrix)
    :
    gsMaterialMatrixTFT(memory::make_shared_not_owned(materialMatrix))
    {}

    /**
     * @brief      Constructs a TFT material matrix using a reference to another material matrix
     *
     * @param      materialMatrix  The original material matrix
     */
    gsMaterialMatrixTFT(gsMaterialMatrixBase<T> & materialMatrix)
    :
    gsMaterialMatrixTFT(memory::make_shared(materialMatrix.clone().release()))
    {}

    /**
     * @brief      Constructs a TFT material matrix using a reference to another material matrix
     *
     * @param      materialMatrix  The original material matrix
     */
    gsMaterialMatrixTFT(const material_ptr & materialMatrix)
    :
    m_materialMat(give(materialMatrix))
    {
        if (materialMatrix.get()->hasUndeformed())
            this->setUndeformed(materialMatrix.get()->getUndeformed());
        if (materialMatrix.get()->hasDeformed())
            this->setDeformed(materialMatrix.get()->getDeformed());

        this->setThickness(materialMatrix.get()->getThickness());

        m_options.addReal("SlackMultiplier","Multiplies the original value of the matrix for the slack state",0);
        m_options.addReal("Eta","Residual stiffness: fraction of the unmodified stress and tangent kept at wrinkled points, slack multiplier at slack points when SlackMultiplier is 0 (Zhang & Kiendl 2026, Eqs. 35, 40)",0);
        m_options.addSwitch("EtaTangentOnly","Apply Eta to the tangent only; stress and strain remain those of the pure tension field",false);
        m_options.addSwitch("Explicit","Explicit iterations; use tension field from a fixed deformed geometry that does not change when calling setDeformed",false);
    }

    /// See \ref gsMaterialMatrixBase for details
    /// Here, we use the integrated matrix and vector to make our computations!
    inline enum MatIntegration isMatIntegrated() const override
    // {return m_materialMat->isMatIntegrated(); gsWarn<<"Change this to a constant one, like for linear. So that no things out of the mid-plane are needed"; }
    {return MatIntegration::Linear;} // Multiplies with z t or its moment

    /// See \ref gsMaterialMatrixBase for details
    /// Here, we use the integrated matrix and vector to make our computations!
    inline enum MatIntegration isVecIntegrated() const override
    // {return m_materialMat->isVecIntegrated(); gsWarn<<"Change this to a constant one, like for linear. So that no things out of the mid-plane are needed"; }
    {return MatIntegration::Linear;} // Multiplies with z t or its moment

    // bool isLinear() { return isLinear_impl(); }
    // typename util::enable_if<U::Linear,bool>::type
    // eval_impl(const U & u, const index_t k)
    // { return m_materialMat; }

    /// See \ref gsMaterialMatrixBase for details
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const override
    { m_materialMat->density_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const override
    { m_materialMat->thickness_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_matrix;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_vector;

    /// Principal stresses of the modified (tension-field) stress, in every state: ascending eigenvalues
    /// (2 rows, row 0 = minor) of the pencil S^{ab} n_b = lambda G^{ab} n_b with G the reference metric of the point:
    /// the fibre metric G(z) for linear=false and the mid-surface metric A at every z for linear=true.
    ///
    /// @param[in]  patch  The patch
    /// @param[in]  u      Evaluation points
    /// @param[in]  z      Evaluation thickness
    /// @param[in]  out    (passed to \ref eval3D_stress)
    ///
    /// @return     2 x (u.cols()*z.rows()) matrix
    ///
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const override;
    using Base::eval3D_strain;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;
    using Base::eval3D_stress;

    /// Cauchy stress from the modified (tension-field) stress: S'/J with J the wrapped material's eval3D_detF
    gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;
    using Base::eval3D_CauchyStress;

    /// Principal Cauchy stresses from the modified (tension-field) stress: ascending eigenvalues
    /// (2 rows, row 0 = minor) of (1/J) g_cov S'^con, g_cov the deformed metric
    gsMatrix<T> eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;
    using Base::eval3D_CauchyPStress;

    /// Principal stress frame of the modified (tension-field) stress in every state: 9 rows, the column-major
    /// 3x3 matrix [d_min | d_max | N] of Cartesian unit vectors (ascending, column 0 = minor),
    /// N = normalise(G_1 x G_2) the unit normal of the reference configuration;
    /// d_min and d_max are Cartesian unit vectors in its tangent plane. The reference basis is the mid-surface
    /// basis for linear=true in dim 3 and the fibre basis otherwise
    gsMatrix<T> eval3D_pstressDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_pstressDir;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_tensionfield;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_theta(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_theta;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_gamma(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_gamma;

    /// Unit vector along the wrinkles on the deformed surface (3 rows); (0,0,0) unless the point is wrinkled
    gsMatrix<T> eval3D_wrinkleDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_wrinkleDir;

    /// Wrinkling intensity max(0,-E_min) of the total strain (1 row, E_min the minor principal value); 0 unless wrinkled
    gsMatrix<T> eval3D_wrinkleIntensity(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_wrinkleIntensity;

    /// Parametrisation-independent wrinkle strain gamma G^{ab} n_a n_b (1 row); NaN unless wrinkled
    gsMatrix<T> eval3D_gammaHat(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_gammaHat;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstretch(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override
    { return m_materialMat->eval3D_pstretch( patch,u,z ); }

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstretchDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const override
    { return m_materialMat->eval3D_pstretchDir( patch,u,z ); }

    /// See \ref gsMaterialMatrixBase for details
    void setParameters(const std::vector<gsFunctionSet<T>*> &pars) override
    { m_materialMat->setParameters(pars); }
    using Base::setParameters;

    /// See \ref gsMaterialMatrixBase for details
    void info() const override
    { m_materialMat->info(); }

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

    /// Updates the reference to the deformed patches used for TFT (explicit only)
    void updateDeformed(const gsFunctionSet<T> * deformed)
    {
        function_ptr fun = memory::make_shared_not_owned(deformed);
        m_defpatches0 = fun;
    }

    /// Updates the reference to the deformed patches used for TFT (explicit only)
    void updateDeformed(const function_ptr & deformed)
    {
        m_defpatches0 = deformed;
    }

    /// Computes theta
    gsMatrix<T> eval_theta(const gsMatrix<T> & Cs, const gsMatrix<T> & Ns, const gsMatrix<T> & Es) const;

    /// See \ref gsMaterialMatrixBase for details
    std::ostream &print(std::ostream &os) const override;

    /// See \ref gsMaterialMatrixBase for details
    const gsMaterialMatrixBase<T> * material() const override { return m_materialMat.get(); }

    /// See \ref gsMaterialMatrixBase for details
    gsMaterialMatrixBase<T> * material() override { return m_materialMat.get(); }

protected:
    template <bool _linear>
    typename std::enable_if< _linear, gsMatrix<T> >::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    template <bool _linear>
    typename std::enable_if<!_linear, gsMatrix<T> >::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    gsMatrix<T> _compute_TF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;
    gsMatrix<T> _compute_TF(const index_t patch, const gsVector<T> & u, const T & z)           const;

    /// Wrinkling state (1 taut, 0 wrinkled, -1 slack) for linear=true, classified on the mid-surface from the pair
    /// (N,E) of the wrapped law with the criterion selected by its TensionField option. Independent of z.
    gsMatrix<T> _compute_TF_mid(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    gsMatrix<T> _compute_E(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & E) const;

    gsMatrix<T> _compute_S(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S) const;

    gsMatrix<T> _compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S) const;

    gsMatrix<T> _compute_C(const T theta, const gsMatrix<T> & C, const gsMatrix<T> & S, const gsMatrix<T> & dC) const;


    /// Tangent C (9x1), stress N (3x1) and engineering strain E (3x1) of the wrapped material at (u,z), the triple
    /// that the wrinkle (eval_theta, gamma, modified stress and tangent) is built from.
    ///
    /// For linear=true the triple is taken on the mid-surface (z is not used): N is the wrapped VectorN stress, which
    /// is the stress of the mid-surface strain E, so that (N,E) is a consistent pair for the wrinkle equations. For
    /// linear=false the wrapped law is queried at height z.
    void _wrinkleState(const index_t patch, const gsVector<T> & u, const T z, gsMatrix<T> & C, gsMatrix<T> & N, gsMatrix<T> & E) const;

    /// theta and gamma at one wrinkled point (the (C,N,E) triple of _wrinkleState: mid-surface for linear=true, at (u,z) otherwise; eval_theta root; _compute_gamma)
    void _wrinkleThetaGamma(const index_t patch, const gsVector<T> & u, const T z, T & theta, T & gamma) const;

    /// Reference/deformed 2x2 metrics and 3x2 Cartesian covariant bases of point k at height z (after _computePoints).
    /// The bases are the mid-surface ones for linear=true (dim 3) and the fibre ones otherwise.
    void _pointGeometry(const index_t k, const T z, gsMatrix<T> & Gori, gsMatrix<T> & Gdef, gsMatrix<T> & Gb, gsMatrix<T> & gb) const;

    /// Option Eta, checked to lie in [0,1].
    T _eta() const;

    /// Factor applied to stress, tangent and strain at slack points for the effective residual fraction \a etaEff:
    /// SlackMultiplier if etaEff == 0 or SlackMultiplier > 0, else etaEff.
    T _slackFactor(const T etaEff) const;

    T _compute_gamma(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;
    gsVector<T> _theta_interval(const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;
    bool _check_theta_full(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;
    bool _check_theta_gamma(const T & theta, const gsMatrix<T> & C, const gsVector<T> & N, const gsVector<T> & E) const;

protected:

    mutable material_ptr m_materialMat;
    using Base::m_patches;
    using Base::m_defpatches;
    function_ptr m_defpatches0;
    using Base::m_options;

    using Base::m_data;

};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixTFT.hpp)
#endif
