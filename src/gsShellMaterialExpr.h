/** @file gsShellMaterialExpr.h

    @brief Lightweight, ZERO-COPY expression views onto the six moment blocks
           (A/B/C/D/N/M) of a \ref gsShellMaterialProvider per-element cache.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

// Only meaningful together with gsShellMaterialProvider, i.e. when the
// gsPhaseFieldFracture module is enabled. The macro lives in the generated
// gsCore/gsConfigExt.h, pulled in by any gismo header above the guard
// (pattern of gsShellMaterialProvider.h:17-25).
#include <gsCore/gsLinearAlgebra.h>

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsExpressions/gsExpressions.h>

// MaterialOutput lives in gsKLShell/src/gsMaterialMatrixUtils.h, but that header
// and gsMaterialMatrixBase.h include EACH OTHER (Utils:18 <-> Base:19). Entering
// the cycle through Utils parses the Base class body before <gsIO/gsOptionList.h>
// and before the MatIntegration enum, which does not compile; entering through
// Base is the order every other gsKLShell header uses. Do not "simplify" this to
// a direct Utils include.
#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsShellMaterialProvider.h>

namespace gismo
{
namespace expr
{

/**
 * @brief   Compile-time descriptor of ONE moment block inside a
 *          \ref gsShellMaterialProvider result column.
 *
 * The 42 rows of column \a k are
 * \f$[\,A(9);\,B(9);\,C(9);\,D(9);\,N(3);\,M(3)\,]\f$, the 3x3 blocks stored
 * COLUMN-MAJOR (gsShellMaterialProvider.h:206-217):
 *
 * | out                    | offset | shape | request bit    |
 * |------------------------|--------|-------|----------------|
 * | MaterialOutput::MatrixA | 0      | 3 x 3 | \c ShellReq_A |
 * | MaterialOutput::MatrixB | 9      | 3 x 3 | \c ShellReq_B |
 * | MaterialOutput::MatrixC | 18     | 3 x 3 | \c ShellReq_C |
 * | MaterialOutput::MatrixD | 27     | 3 x 3 | \c ShellReq_D |
 * | MaterialOutput::VectorN | 36     | 3 x 1 | \c ShellReq_N |
 * | MaterialOutput::VectorM | 39     | 3 x 1 | \c ShellReq_M |
 *
 * The offsets are never spelled as literals here: they are read from the
 * provider's own enum, so the table cannot drift away from the producer.
 *
 * The primary template is deliberately left UNDEFINED: only the six outputs the
 * provider actually integrates have a specialization, so instantiating a view
 * on any other \ref MaterialOutput is a COMPILE error rather than a silent
 * out-of-range map.
 *
 * @tparam out The requested moment
 * @tparam dim The shell embedding dimension (2 = planar, 3 = surface)
 * @tparam T   Real type
 */
template<enum MaterialOutput out, short_t dim, class T>
struct shellMomentBlock;

template<short_t dim, class T>
struct shellMomentBlock<MaterialOutput::MatrixA,dim,T>
{
    typedef gsShellMaterialProvider<dim,T> Provider;
    enum { Offset  = Provider::OffsetA, Rows = 3, Cols = 3,
           Request = Provider::ShellReq_A };
    static const char * name() { return "shellMatA"; }
};

template<short_t dim, class T>
struct shellMomentBlock<MaterialOutput::MatrixB,dim,T>
{
    typedef gsShellMaterialProvider<dim,T> Provider;
    enum { Offset  = Provider::OffsetB, Rows = 3, Cols = 3,
           Request = Provider::ShellReq_B };
    static const char * name() { return "shellMatB"; }
};

template<short_t dim, class T>
struct shellMomentBlock<MaterialOutput::MatrixC,dim,T>
{
    typedef gsShellMaterialProvider<dim,T> Provider;
    enum { Offset  = Provider::OffsetC, Rows = 3, Cols = 3,
           Request = Provider::ShellReq_C };
    static const char * name() { return "shellMatC"; }
};

template<short_t dim, class T>
struct shellMomentBlock<MaterialOutput::MatrixD,dim,T>
{
    typedef gsShellMaterialProvider<dim,T> Provider;
    enum { Offset  = Provider::OffsetD, Rows = 3, Cols = 3,
           Request = Provider::ShellReq_D };
    static const char * name() { return "shellMatD"; }
};

template<short_t dim, class T>
struct shellMomentBlock<MaterialOutput::VectorN,dim,T>
{
    typedef gsShellMaterialProvider<dim,T> Provider;
    enum { Offset  = Provider::OffsetN, Rows = 3, Cols = 1,
           Request = Provider::ShellReq_N };
    static const char * name() { return "shellVecN"; }
};

template<short_t dim, class T>
struct shellMomentBlock<MaterialOutput::VectorM,dim,T>
{
    typedef gsShellMaterialProvider<dim,T> Provider;
    enum { Offset  = Provider::OffsetM, Rows = 3, Cols = 1,
           Request = Provider::ShellReq_M };
    static const char * name() { return "shellVecM"; }
};

/**
 * @brief   A lightweight, ZERO-COPY view of ONE moment block of a
 *          \ref gsShellMaterialProvider.
 *
 * Replaces \c reshape(getCoeff(gsMaterialMatrixIntegrate<...>),3,3) in the
 * Kirchhoff-Love assembler algebra: the six views share ONE batched per-element
 * sweep instead of paying the full geometric precomputation six times
 * (gsKLShell issue #28).
 *
 * ### What eval() returns
 * A \c gsAsConstMatrix mapped straight onto the provider's 42 x N cache -- no
 * copy, no temporary. The view is natively 3x3 (A/B/C/D) or 3x1 (N/M), so the
 * assembler must NOT wrap it in \c reshape: \c reshape_expr hard-asserts the
 * total size (gsExpressions/reshape_expr.h:50) and a 42-row coefficient would
 * trip it.
 *
 * ### Protocol (mirrors the solids \c materialView_expr)
 *  - \ref parse registers the provider symbol and BOTH geometry maps, then
 *    hands the CALLING thread's map-data addresses to the provider through
 *    \c bindKinematics. This re-binding happens on EVERY parse, unconditionally:
 *    \c gsExprHelper::cleanUp destroys and re-allocates the map data, so an
 *    address kept from a previous parse is dangling. A stale bind is in fact the
 *    ONE scenario that can decouple the provider's point count from the maps'.
 *  - \c eval calls \c ensureFilled, which runs the per-element sweep on the
 *    FIRST view evaluation of the element and early-returns for every later
 *    call (of this and of the five sibling views).
 *
 * ### Two hard constraints inherited from the provider
 * 1. The provider is registered as a PLAIN variable (\c getCoeff / \c getVar of
 *    the function set), NEVER inside a composition: a composition lands in
 *    \c m_cdata and is evaluated at PHYSICAL points, which would then be
 *    forwarded to \c precomputeParameters -- which expects PARAMETRIC ones
 *    (gsExprHelper.h:458-470). This view only ever calls \c _coeff.parse().
 * 2. \c NEED_DERIV is never set on the PROVIDER symbol (only on the two maps).
 *    The provider is a \c gsFunctionSet without an analytic derivative, so
 *    \c gsFunction::deriv_into would FINITE-DIFFERENCE it, re-entering
 *    \c eval_into -- hence \c beginElement -- at perturbed points and thrashing
 *    the per-element cache. \c symbol_expr::parse sets \c NEED_VALUE|NEED_ACTIVE
 *    and nothing else, and no flag is added here.
 *
 * ### OpenMP: the maps are held BY VALUE
 * \c gsGeometryMap and \c gsFeVariable resolve \c data() through a member
 * pointer that \c gsExprHelper::add re-binds to the CALLING thread's slot. An
 * expression that stored a map by reference would therefore have every thread
 * rebind the SAME shared object during its own parse (measured: ~60-100 % crash
 * rate at OMP_NUM_THREADS=4 in the solids twin). The maps are consequently owned
 * per instance and the copy constructor DEEP-COPIES them, so that each thread's
 * expression copy owns distinct maps whose \c data() resolves to that thread's
 * own \c gsMapData -- which is exactly what \c bindKinematics must inject.
 *
 * A consequence worth knowing: \ref parse binds the view's OWN map copies, not
 * the \c gsGeometryMap objects the caller passed in. Both resolve to the same
 * per-function-set slot of \c gsExprHelper::m_mdata, so the two agree as soon as
 * the caller's maps are part of the parsed expression as well -- which they
 * always are in an assembly (the integrand uses them). A view is not a way to
 * get a map computed on its own.
 *
 * ### Complexity
 * O(1) per quadrature point: the sweep (O(N) geometry + ONE batched
 * condensation over N*nz points) runs once per element, shared by all six views.
 *
 * @tparam out The moment to view (MatrixA/B/C/D, VectorN, VectorM)
 * @tparam dim The shell embedding dimension (2 = planar, 3 = surface)
 * @tparam T   Real type
 *
 * @ingroup KLShell
 */
template<enum MaterialOutput out, short_t dim, class T>
class shellMaterialView_expr : public _expr<shellMaterialView_expr<out,dim,T> >
{
public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

private:
    typedef gsShellMaterialProvider<dim,T> Provider;
    typedef shellMomentBlock<out,dim,T>    Block;

    /// The provider registered as a plain variable. \c gsFeVariable::Nested_t
    /// is BY VALUE, so this instance owns its own symbol (per-thread copies).
    typename gsFeVariable<T>::Nested_t _coeff;

    /// The undeformed and deformed geometry maps, owned BY VALUE (see the class
    /// doc on OpenMP). \c gsGeometryMap has no accessible default constructor,
    /// hence \c unique_ptr rather than plain value members.
    memory::unique_ptr<gsGeometryMap<T> > _Gori;
    memory::unique_ptr<gsGeometryMap<T> > _Gdef;

    const Provider * _provider;

public:

    /**
     * @brief   Constructor.
     *
     * @param[in] coeff    The provider, registered as a PLAIN variable
     * @param[in] Gori     The UNDEFORMED midsurface map
     * @param[in] Gdef     The DEFORMED midsurface map (may be the same map as
     *                     @a Gori -- a linear assembly does exactly that)
     * @param[in] provider The provider itself (non-owning)
     */
    shellMaterialView_expr(const gsFeVariable<T>  & coeff,
                           const gsGeometryMap<T> & Gori,
                           const gsGeometryMap<T> & Gdef,
                           const Provider         * provider)
    :
    _coeff(coeff),
    _Gori(new gsGeometryMap<T>(Gori)),
    _Gdef(new gsGeometryMap<T>(Gdef)),
    _provider(provider)
    {
        GISMO_ENSURE(_provider!=nullptr,
                     "shellMaterialView: null provider; this view has no unbound fallback mode.");
    }

    /// Deep-copying copy constructor: every copy owns its OWN geometry maps.
    /// Required for per-thread data resolution -- see the class documentation.
    shellMaterialView_expr(const shellMaterialView_expr & other)
    :
    _coeff(other._coeff),
    _Gori(other._Gori ? new gsGeometryMap<T>(*other._Gori) : nullptr),
    _Gdef(other._Gdef ? new gsGeometryMap<T>(*other._Gdef) : nullptr),
    _provider(other._provider)
    { }

    /// Zero-copy map of this view's block at quadrature point @a k. The block
    /// is contiguous inside the (column-major) cache column, so it is mapped in
    /// place. O(1): the per-element sweep is shared by all six views.
    const gsAsConstMatrix<T> eval(const index_t k) const
    {
        // Non-requested rows are ZERO, not stale (the provider zeroes the whole
        // 42 x N buffer), so a mask/view mismatch would surface as a silently
        // vanishing stiffness contribution. This assert is the ONLY detector.
        // NOTE it is live in this build: gsDebug.h:105-110 gates GISMO_ASSERT on
        // #ifndef NDEBUG and CMAKE_CXX_FLAGS_RELWITHDEBINFO carries no -DNDEBUG.
        GISMO_ASSERT(_provider->requested((unsigned)Block::Request),
                     "The view "<<Block::name()<<" is evaluated, but its moment is not in the "
                     "provider's request mask ("<<_provider->requestMask()<<"): the block was "
                     "never integrated (see gsShellMaterialProvider::setRequested).");
        _provider->ensureFilled();
        const gsMatrix<T> & M = _provider->cachedMoments();
        GISMO_ASSERT(M.rows()==(index_t)Provider::TargetDim,
                     "The provider cache has "<<M.rows()<<" rows, expected "
                     <<(index_t)Provider::TargetDim<<".");
        GISMO_ASSERT(k>=0 && k<M.cols(),
                     "Point index "<<k<<" out of range: the element cached "<<M.cols()<<" points.");
        return gsAsConstMatrix<T>(M.data() + k*(index_t)Provider::TargetDim + (index_t)Block::Offset,
                                  (index_t)Block::Rows, (index_t)Block::Cols);
    }

    index_t rows() const { return (index_t)Block::Rows; }
    index_t cols() const { return (index_t)Block::Cols; }

    /**
     * @brief   Registers the provider symbol and the two geometry maps, then
     *          injects THIS thread's map data into the provider.
     *
     * @warning The \c bindKinematics call is unconditional and must stay so: the
     *          addresses handed over are invalidated by every
     *          \c gsExprHelper::cleanUp (which destroys and re-allocates the map
     *          data), so they have to be refreshed on every parse.
     */
    void parse(gsExprHelper<T> & evList) const
    {
        // The provider as a PLAIN variable: NEED_VALUE|NEED_ACTIVE only (see the
        // class doc, constraints 1 and 2). Its "evaluation" is the element
        // -boundary signal that stores the PARAMETRIC quadrature points.
        _coeff.parse(evList);

        // Both maps: registered so that their per-thread map data is computed
        // for every element, with the flag word gsShellKinematics and the
        // thickness evaluation need (gsShellMaterialProvider.h:123-131).
        evList.add(*_Gori);
        evList.add(*_Gdef);
        _Gori->data().flags |= _mapFlags();
        _Gdef->data().flags |= _mapFlags();

        // The CALLING thread's addresses, refreshed at every parse.
        _provider->bindKinematics(&_Gori->data(), &_Gdef->data());
    }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream & os) const { os << Block::name(); }

private:

    /// Flag word required on BOTH injected maps. NEED_VALUE is MANDATORY: the
    /// provider derives its point count from \c values[0].cols() (\c map.points
    /// is empty after \c gsExprHelper::precompute) and evaluates the thickness
    /// at the undeformed map values. NEED_JACOBIAN == NEED_DERIV; both are named
    /// because both appear in the two contracts. NEED_NORMAL is inert for
    /// \c dim==2 (gsFunction.hpp:749 computes normals only when
    /// \c tarDim == domDim+1), so one uniform word is safe in both dimensions.
    static unsigned _mapFlags()
    {
        return (unsigned)(NEED_VALUE | NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL)
             | (3==dim ? (unsigned)NEED_2ND_DER : 0u);
    }
};

/**
 * @brief   Factory for \ref shellMaterialView_expr.
 *
 * @a dim and @a T are deduced from @a provider; only the moment is explicit:
 * \code
 * auto A = shellMaterialView<MaterialOutput::MatrixA>(pv, G, Gdef, &provider);
 * \endcode
 *
 * @param[in] coeff    The provider registered as a PLAIN variable
 * @param[in] Gori     The undeformed midsurface map
 * @param[in] Gdef     The deformed midsurface map
 * @param[in] provider The provider (non-owning)
 */
template<enum MaterialOutput out, short_t dim, class T>
EIGEN_STRONG_INLINE
shellMaterialView_expr<out,dim,T>
shellMaterialView(const gsFeVariable<T>                 & coeff,
                  const gsGeometryMap<T>                & Gori,
                  const gsGeometryMap<T>                & Gdef,
                  const gsShellMaterialProvider<dim,T>  * provider)
{ return shellMaterialView_expr<out,dim,T>(coeff,Gori,Gdef,provider); }

// The same factory, spelled out once per output. These are one-line forwards:
// they exist only so that assembler algebra reads shellMatrixA(...) instead of
// shellMaterialView<MaterialOutput::MatrixA>(...).

/// Membrane stiffness A (3x3), see \ref shellMaterialView
template<short_t dim, class T> EIGEN_STRONG_INLINE
shellMaterialView_expr<MaterialOutput::MatrixA,dim,T>
shellMatrixA(const gsFeVariable<T> & coeff, const gsGeometryMap<T> & Gori,
             const gsGeometryMap<T> & Gdef, const gsShellMaterialProvider<dim,T> * provider)
{ return shellMaterialView<MaterialOutput::MatrixA>(coeff,Gori,Gdef,provider); }

/// Coupling stiffness B (3x3), see \ref shellMaterialView
template<short_t dim, class T> EIGEN_STRONG_INLINE
shellMaterialView_expr<MaterialOutput::MatrixB,dim,T>
shellMatrixB(const gsFeVariable<T> & coeff, const gsGeometryMap<T> & Gori,
             const gsGeometryMap<T> & Gdef, const gsShellMaterialProvider<dim,T> * provider)
{ return shellMaterialView<MaterialOutput::MatrixB>(coeff,Gori,Gdef,provider); }

/// Coupling stiffness C (3x3, equal to B), see \ref shellMaterialView
template<short_t dim, class T> EIGEN_STRONG_INLINE
shellMaterialView_expr<MaterialOutput::MatrixC,dim,T>
shellMatrixC(const gsFeVariable<T> & coeff, const gsGeometryMap<T> & Gori,
             const gsGeometryMap<T> & Gdef, const gsShellMaterialProvider<dim,T> * provider)
{ return shellMaterialView<MaterialOutput::MatrixC>(coeff,Gori,Gdef,provider); }

/// Bending stiffness D (3x3), see \ref shellMaterialView
template<short_t dim, class T> EIGEN_STRONG_INLINE
shellMaterialView_expr<MaterialOutput::MatrixD,dim,T>
shellMatrixD(const gsFeVariable<T> & coeff, const gsGeometryMap<T> & Gori,
             const gsGeometryMap<T> & Gdef, const gsShellMaterialProvider<dim,T> * provider)
{ return shellMaterialView<MaterialOutput::MatrixD>(coeff,Gori,Gdef,provider); }

/// Normal force N (3x1), see \ref shellMaterialView
template<short_t dim, class T> EIGEN_STRONG_INLINE
shellMaterialView_expr<MaterialOutput::VectorN,dim,T>
shellVectorN(const gsFeVariable<T> & coeff, const gsGeometryMap<T> & Gori,
             const gsGeometryMap<T> & Gdef, const gsShellMaterialProvider<dim,T> * provider)
{ return shellMaterialView<MaterialOutput::VectorN>(coeff,Gori,Gdef,provider); }

/// Bending moment M (3x1), see \ref shellMaterialView
template<short_t dim, class T> EIGEN_STRONG_INLINE
shellMaterialView_expr<MaterialOutput::VectorM,dim,T>
shellVectorM(const gsFeVariable<T> & coeff, const gsGeometryMap<T> & Gori,
             const gsGeometryMap<T> & Gdef, const gsShellMaterialProvider<dim,T> * provider)
{ return shellMaterialView<MaterialOutput::VectorM>(coeff,Gori,Gdef,provider); }

} // namespace expr
} // namespace gismo

#endif // gsPhaseFieldFracture_ENABLED
