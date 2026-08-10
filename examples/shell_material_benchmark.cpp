/** @file shell_material_benchmark.cpp

    @brief Baseline benchmark exposing the gsThinShellAssembler material redundancy,
           extended with the gsMaterialMatrix3D adapter end-to-end comparison.

    This driver is the baseline for gsKLShell issue #28 (material unification).
    The legacy gsThinShellAssembler material pipeline builds SIX independent
    gsMaterialMatrixIntegrate coefficients per assembly:
        - assembleMatrix: MatrixA, MatrixB, MatrixC, MatrixD  (4 outputs)
        - assembleVector: VectorN, VectorM                    (2 outputs)
    Each of these outputs independently re-runs
    gsMaterialMatrixBaseDim::_computePoints (both metrics + thickness + material
    parameters) on every element block, and -- for compressible nonlinear
    materials -- its own through-thickness sweep with per-point plane-stress
    C33 Newton iterations.

    This example MEASURES that redundancy and the effect of the same-input guard
    (issue #28 task 11) that eliminates it:
      * wall-clock assembly times (gsStopwatch);
      * two instrumented counters per assembly (gsMaterialMatrixBaseDim.h):
          - CALLS  : _computePoints INVOCATIONS (counted before the guard);
          - MISSES : full recomputations that passed the guard (the honest cost).

    Expected CALLS structure (single-threaded, unchanged by the guard):
      * per assembleMatrix ~ 4 x nElements  (A,B,C,D)
      * per assembleVector ~ 2 x nElements  (N,M)
      * per assemble()     ~ 2 x nElements (SvK) (measured; NH amplifies to ~8)
    A full Newton assembly (assembleMatrix + assembleVector) therefore INVOKES
    ~6 x nElements for SvK. With the same-input guard the MISSES collapse to
    ~1 x nElements per assembly kind: the six A/B/C/D/N/M coefficients share a
    single metric/thickness/parameter computation per element per assembly.
    The CALLS count is the metric-recomputation redundancy factor and looks
    (nearly) the same for the linear (SvK) and nonlinear (NH) material; the
    ADDITIONAL cost of the per-output through-thickness C33 Newton sweep (which is
    NOT inside _computePoints and so NOT removed by the guard) shows up in the
    SvK<->NH wall-time gap instead.

    ### Adapter mode (-a 1, issue #28 task 16, requires the PFF module)

    With `-a 1` the SAME assembler protocol is driven by the gsMaterialMatrix3D
    adapter, which wraps ONE 3D gsPhaseFieldFracture (PFF) constitutive law behind
    the classic gsMaterialMatrixBase interface (plane-stress condensation on the
    3D law, per-thread cross-output cache). `-a 0` reproduces the task-10/11 legacy
    output byte-for-byte (same lines) so the baselines stay comparable; ALL new
    output (the legacy-vs-adapter K/R agreement check and the adapter sweep/hit
    counters) lives under `-a 1`.

    Matched laws (task 15 pinned the twins numerically):
      * -m 0 : PFF gsLinearMaterial          (SvK oracle)
      * -m 1 : PFF gsNeoHookeQuadMaterial     (== legacy compressible NH; the
               QuadraticVolumetric alias U=(k/4)(J^2-1-2lnJ), matched to 6e-10).
    The PFF law's parameter functions are constants on the PARAMETRIC domain
    (dim 2 for a surface): the adapter evaluates them at the parametric points u,
    so they are built with domain dim 2 (NOT 3 -- a dim-3 gsConstantFunction would
    trip the u.rows()==domainDim assertion when sampled at 2-row parametric
    points). The legacy material keeps its dim-3 parameter functions (it samples at
    physical points). See task-16 report.

    Cross-mode agreement (-a 1, printed once before timing): K and R are assembled
    on the SAME deformed state with the legacy material and with the adapter, and
    ||K_legacy - K_adapter|| / ||K_legacy|| (and the same for R) is printed. These
    are NOT hard-asserted -- they are the honest end-to-end model-difference /
    parity numbers:
      * SvK  : an O(t^2/R^2) MODEL difference (legacy integrates a z-constant
               tangent analytically; the adapter Gauss-integrates the exact
               z-dependent metric). Task-15 measured MatrixA rel-diff 8.4e-6 at
               t=0.01 on a flatter fixture; here t=0.25 on R~25 gives O(t^2/R^2)
               ~ 1e-4, so a LARGER number is expected and correct.
      * NH   : the model twin is exact, so agreement is PARITY (<= ~1e-8) modulo
               the shared through-thickness quadrature.

    ### Provider mode (-a 2, issue #28 task 30, requires the PFF module)

    `-a 2` runs the Phase-4 pipeline: gsThinShellAssembler2 driven by ONE
    gsShellMaterialProvider (the six moments A/B/C/D/N/M come from a single
    batched per-element sweep, with a per-ROUTINE request mask) over the SAME
    PFF law, on the SAME problem, with the SAME through-thickness rule
    (NumGauss; the equality of the two values is ASSERTED, not assumed). It is
    a SELF-CONTAINED branch and returns before the -a 0 / -a 1 code, so their
    output is bit-for-bit untouched.

    Because the Phase-4 question is a COMPARISON, `-a 2` times BOTH the adapter
    path and the provider path in the SAME process:
      * every timed instance first runs an UNTIMED warm-up assembly (the
        cross-mode agreement pass below), so no timed rep pays first-touch cost;
      * the reps are INTERLEAVED and their order is ALTERNATED per rep, so that
        neither path can profit from running second (allocator/page warm-up and
        clock drift affect a ~20 ms measurement). `--sequential` switches to two
        back-to-back blocks instead: the control that shows whether interleaving
        biased the comparison (report both, they must agree);
      * mean, min and max over the reps are printed for both paths.

    The REQUIRED gate of task 30 is `assembleVector`: the ADAPTER's residual
    path pays a full constitutive sweep that its own matrix path already paid
    for (setDeformed bumps the config revision and splits the cross-output
    cache), whereas the PROVIDER's residual path requests only the stress
    moments N and M -- so matrixMomentFills() must be 0 after it, which this
    driver ASSERTS. The expected (non-required) gate is assembleMatrix <=
    adapter. Both verdicts are printed as one line each, whichever way they go.

    Geometry: single-patch Scordelis-Lo roof (a curved shell from the main
    gismo filedata; guaranteed on the default search path). This is the "quarter
    cylinder" case the task lists as acceptable; the redundancy measured here is
    geometry-independent.

    Measurement guidance: run single-threaded (OMP_NUM_THREADS=1); the counters are
    plain non-atomic globals and are only thread-exact for a single thread. The
    OMP_NUM_THREADS=4 run of `-m 1 -a 1` is a thread-STABILITY check only (per-
    thread caches make the sweep count scale roughly with the thread count). The
    same holds for the -a 2 counters (gsShellMaterialProvider::fillCount is a
    plain non-atomic member): the -a 2 counter identities are ASSERTED only when
    the team size is 1 and are printed-but-not-asserted otherwise.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>
#include <gsKLShell/gsKLShell.h>

#ifdef gsPhaseFieldFracture_ENABLED
#include <gsKLShell/src/gsMaterialMatrix3D.h>
// gsShellMaterialExpr.h (pulled in by gsThinShellAssembler2.h) only
// forward-declares gsExprHelper, so a real gsExprAssembler declaration must come
// first (same convention as the gsThinShellAssembler2 unit tests).
#include <gsAssembler/gsExprAssembler.h>
#include <gsKLShell/src/gsThinShellAssembler2.h>
#include <gsPhaseFieldFracture/materials/gsLinearMaterial.h>
#include <gsPhaseFieldFracture/materials/gsNeoHookeQuadMaterial.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

using namespace gismo;

#ifdef gsPhaseFieldFracture_ENABLED
/// Mean / min / max of a set of per-rep wall times (-a 2 only). The spread is
/// reported next to the mean because earlier measurements in this project
/// showed ~30% run-to-run scatter on the legacy path: a bare mean is not
/// enough to judge a head-to-head difference.
struct gsTimeStats
{
    real_t mean, min, max;
};

gsTimeStats gsComputeTimeStats(const std::vector<real_t> & v)
{
    GISMO_ASSERT(!v.empty(),"empty timing vector");
    gsTimeStats s;
    s.mean = 0.0;
    s.min  = v[0];
    s.max  = v[0];
    for (size_t i=0; i!=v.size(); ++i)
    {
        s.mean += v[i];
        if (v[i]<s.min) s.min = v[i];
        if (v[i]>s.max) s.max = v[i];
    }
    s.mean /= (real_t)v.size();
    return s;
}

/// ONE timed assembly call (-a 2 only).
/// @param path  1 = gsMaterialMatrix3D adapter + legacy gsThinShellAssembler,
///              2 = gsShellMaterialProvider + gsThinShellAssembler2
/// @param kind  0 = assemble(), 1 = assembleMatrix(def), 2 = assembleVector(def)
/// @return the wall time of that single call
real_t gsTimedAssembly(index_t path, index_t kind,
                       gsThinShellAssembler <3,real_t,true> & adp,
                       gsThinShellAssembler2<3,real_t,true> & prov,
                       const gsMultiPatch<real_t> & def)
{
    gsStopwatch clk;
    ThinShellAssemblerStatus st;
    real_t dt;
    if (path==1)
    {
        clk.restart();
        st = (kind==0 ? adp.assemble()
            : kind==1 ? adp.assembleMatrix(def)
                      : adp.assembleVector(def));
        dt = clk.stop();
    }
    else
    {
        clk.restart();
        st = (kind==0 ? prov.assemble()
            : kind==1 ? prov.assembleMatrix(def)
                      : prov.assembleVector(def));
        dt = clk.stop();
    }
    GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,
        "assembly failed (path "<<path<<", kind "<<kind<<"): status "<<(index_t)st);
    return dt;
}

/// One timed BLOCK of @a reps calls of routine @a kind on BOTH paths (-a 2).
///
/// Default (@a sequential == false): the reps are INTERLEAVED and their order is
/// ALTERNATED per rep, so that neither path can systematically profit from
/// running second (allocator/page state, clock drift). With @a sequential the
/// two paths run in separate back-to-back blocks instead -- the control that
/// shows whether interleaving biased the comparison. Both orderings are honest;
/// they are reported side by side.
void gsTimeBlock(index_t kind, index_t reps, bool sequential,
                 gsThinShellAssembler <3,real_t,true> & adp,
                 gsThinShellAssembler2<3,real_t,true> & prov,
                 const gsMultiPatch<real_t> & def,
                 std::vector<real_t> & t1, std::vector<real_t> & t2)
{
    t1.assign(reps,0.0);
    t2.assign(reps,0.0);
    if (sequential)
    {
        for (index_t i=0; i!=reps; ++i) t1[i] = gsTimedAssembly(1,kind,adp,prov,def);
        for (index_t i=0; i!=reps; ++i) t2[i] = gsTimedAssembly(2,kind,adp,prov,def);
    }
    else
    {
        for (index_t i=0; i!=reps; ++i)
        {
            if (i%2==0)
            {
                t1[i] = gsTimedAssembly(1,kind,adp,prov,def);
                t2[i] = gsTimedAssembly(2,kind,adp,prov,def);
            }
            else
            {
                t2[i] = gsTimedAssembly(2,kind,adp,prov,def);
                t1[i] = gsTimedAssembly(1,kind,adp,prov,def);
            }
        }
    }
}
#endif

int main(int argc, char *argv[])
{
    //! [Parse command line]
    index_t numRefine  = 3;    // -r
    index_t numElevate = 0;    // -e
    index_t material   = 0;    // -m : 0 = SvK linear, 1 = NH compressible
    index_t adapter    = 0;    // -a : 0 = legacy material, 1 = gsMaterialMatrix3D adapter
    index_t reps       = 3;    // -n : timing repetitions per assembly kind
    real_t  thickness  = 0.25; // ~ 0.01 * L for the Scordelis-Lo roof (L ~ 50)
    real_t  E_modulus  = 4.32e8;
    real_t  PoissonRatio = 0.0;
    bool    sequential = false; // --sequential : -a 2 measurement-bias control

    gsCmdLine cmd("gsKLShell material redundancy baseline benchmark (issue #28).");
    cmd.addInt ("r", "uniformRefine", "Number of uniform h-refinement steps", numRefine);
    cmd.addInt ("e", "degreeElevation", "Number of degree elevation steps", numElevate);
    cmd.addInt ("m", "material", "Material: 0 = SvK linear, 1 = Neo-Hookean compressible", material);
    cmd.addInt ("a", "adapter", "Material driver: 0 = legacy gsMaterialMatrix, 1 = gsMaterialMatrix3D adapter, 2 = gsThinShellAssembler2 + gsShellMaterialProvider (head-to-head vs 1)", adapter);
    cmd.addInt ("n", "reps", "Timing repetitions per assembly kind", reps);
    cmd.addReal("T", "thickness", "Shell thickness", thickness);
    cmd.addReal("E", "Emodulus", "Young's modulus", E_modulus);
    cmd.addReal("v", "PoissonRatio", "Poisson ratio", PoissonRatio);
    cmd.addSwitch("sequential", "(-a 2 only) time the two paths in SEPARATE back-to-back blocks "
                                "instead of interleaving them (measurement-bias control)", sequential);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ENSURE(reps>=1,"reps must be >= 1");
    GISMO_ENSURE(adapter>=0 && adapter<=2,
        "Unknown material driver "<<adapter<<" (use -a 0, -a 1 or -a 2)");

#ifndef gsPhaseFieldFracture_ENABLED
    GISMO_ENSURE(adapter==0,
        "The adapter path (-a 1) and the provider path (-a 2) require the "
        "gsPhaseFieldFracture module (gsMaterialMatrix3D / gsShellMaterialProvider). "
        "Reconfigure with the PFF module enabled.");
#endif

    // -a 1 : the timed assembler is the LEGACY one, driven by the adapter.
    // -a 2 : self-contained provider-vs-adapter head-to-head (returns early).
    const bool useAdapter  = (adapter==1);
    const bool useProvider = (adapter==2);
#ifndef gsPhaseFieldFracture_ENABLED
    GISMO_UNUSED(useAdapter);
    GISMO_UNUSED(useProvider);
#endif

    //! [Read geometry]
    // Single-patch curved shell from the main gismo filedata (Scordelis-Lo roof).
    gsMultiPatch<> mp;
    gsReadFile<>("surfaces/scordelis_lo_roof.xml", mp);
    GISMO_ENSURE(mp.nPatches()>0,"Failed to read geometry 'surfaces/scordelis_lo_roof.xml'");
    if (PoissonRatio==0.0 && material==1)
        PoissonRatio = 0.3; // NH compressible needs a physical Poisson ratio
    //! [Read geometry]

    //! [Refine and elevate]
    if (numElevate!=0)
        mp.degreeElevate(numElevate);
    for (index_t r=0; r<numRefine; ++r)
        mp.uniformRefine();

    // Representative deformed configuration: a small uniform stretch so that the
    // strain is genuinely non-zero (mandatory for the NH C33-Newton to be
    // representative), mirroring how the shell examples build mp_def.
    gsMultiPatch<> mp_def = mp;
    for (size_t p=0; p!=mp_def.nPatches(); ++p)
        mp_def.patch(p).coefs().array() *= 1.01;
    //! [Refine and elevate]

    gsMultiBasis<> dbasis(mp);
    const index_t nElements = dbasis.totalElements();

    //! [Boundary conditions]
    // Assembly-only benchmark: BCs only need to keep assembly well-defined.
    // Scordelis-Lo diaphragm conditions on the straight edges.
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2);
    //! [Boundary conditions]

    //! [Material and forcing]
    gsVector<> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<> force(tmp,3);

    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus), 3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsFunctionExpr<> rho("1.0",3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    // Legacy material (always built: -a 0 uses it for timing, -a 1 uses it as the
    // agreement-check reference).
    gsMaterialMatrixBase<real_t>::uPtr legacyMaterial;
    std::string matName;
    if (material==0)
    {
        // SvK linear
        legacyMaterial = memory::make_unique(new gsMaterialMatrixLinear<3,real_t>(mp,t,E,nu,rho));
        matName = "SvK (linear)";
    }
    else if (material==1)
    {
        // Neo-Hookean, compressible
        gsOptionList options;
        options.addInt   ("Material","Material model",   0); // Material::NH
        options.addSwitch("Compressibility","Compressible",true);
        options.addInt   ("Implementation","Implementation",1); // Analytical
        options.setInt("Material", (index_t)Material::NH);
        legacyMaterial = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        GISMO_ENSURE(legacyMaterial!=nullptr,"getMaterialMatrix returned null");
        matName = "Neo-Hookean (compressible)";
    }
    else
        GISMO_ERROR("Unknown material "<<material<<" (use -m 0 or -m 1)");

#ifdef gsPhaseFieldFracture_ENABLED
    // Adapter material (built for -a 1 and, as the head-to-head reference AND
    // the timed adapter path, for -a 2). The PFF law owns its constant
    // parameter functions on the parametric domain (dim 2 for a surface).
    memory::unique_ptr<gsMaterialBase<real_t> >   pffLaw;
    memory::unique_ptr<gsMaterialMatrix3D<3,real_t> > adapterMaterial;
    std::string adapterName;
    if (useAdapter || useProvider)
    {
        if (material==0)
        {
            pffLaw = memory::make_unique(new gsLinearMaterial<real_t>(E_modulus,PoissonRatio,2));
            adapterName = "gsMaterialMatrix3D[gsLinearMaterial]";
        }
        else // material==1
        {
            // The MATCHED NH twin (task 15): QuadraticVolumetric == legacy comp-NH.
            pffLaw = memory::make_unique(new gsNeoHookeQuadMaterial<real_t>(E_modulus,PoissonRatio,2));
            adapterName = "gsMaterialMatrix3D[gsNeoHookeQuadMaterial]";
        }
        adapterMaterial = memory::make_unique(new gsMaterialMatrix3D<3,real_t>(mp,t,*pffLaw));
    }
#endif

    // Which material drives the timed assembler. NOTE: for -a 0 the Material
    // header line is kept BYTE-IDENTICAL to the task-10/11 baseline (bare matName,
    // no suffix); the mode label is only added on the new -a 1 output.
    gsMaterialMatrixBase<real_t>* activeMaterial = legacyMaterial.get();
    std::string activeName = matName;
#ifdef gsPhaseFieldFracture_ENABLED
    if (useAdapter)
    {
        activeMaterial = adapterMaterial.get();
        activeName = matName + " [adapter " + adapterName + "]";
    }
#endif
    //! [Material and forcing]

#ifdef gsPhaseFieldFracture_ENABLED
    //! [Provider pipeline head-to-head (-a 2, task 30)]
    // SELF-CONTAINED branch: it RETURNS before the -a 0 / -a 1 code below, so
    // that their output stays bit-for-bit what tasks 10/11/16 recorded.
    if (useProvider)
    {
        ThinShellAssemblerStatus st;

#ifdef _OPENMP
        const int nThreads = omp_get_max_threads();
#else
        const int nThreads = 1;
#endif
        // gsShellMaterialProvider::fillCount is a plain non-atomic member: only
        // a team of ONE makes the counter identities exact integers. They are
        // therefore ASSERTED single-threaded and only printed otherwise (a
        // multi-threaded run of this mode is an agreement/stability check, and
        // its timings are meaningless for the same reason).
        const bool exactCounters = (nThreads==1);

        // The three drivers of the SAME problem on the SAME deformed state:
        //   mode 0 : legacy gsMaterialMatrix*       (MODEL reference, not timed)
        //   mode 1 : gsMaterialMatrix3D adapter     (TIMED -- the thing to beat)
        //   mode 2 : gsThinShellAssembler2/provider (TIMED -- the Phase-4 path)
        gsThinShellAssembler <3,real_t,true> legAsm(mp,dbasis,bc,force,legacyMaterial.get());
        gsThinShellAssembler <3,real_t,true> adpAsm(mp,dbasis,bc,force,adapterMaterial.get());
        gsThinShellAssembler2<3,real_t,true> newAsm(mp,dbasis,bc,force,t,*pffLaw);

        // Same through-thickness rule on both timed sides -- ASSERTED, not
        // assumed. gsMaterialMatrix3D registers no "NumGauss" of its own, so the
        // adapter path resolves it inside gsMaterialMatrixIntegrate through
        // askInt("NumGauss",4); gsThinShellAssembler2 owns the option itself.
        const index_t ngAdapter  = adapterMaterial->options().askInt("NumGauss",4);
        const index_t ngProvider = newAsm.options().getInt("NumGauss");
        GISMO_ENSURE(ngAdapter==ngProvider,
            "NumGauss mismatch: adapter "<<ngAdapter<<" vs provider "<<ngProvider
            <<" -- the two paths would not integrate the same problem.");

        gsInfo << "===================================================================\n";
        gsInfo << " gsKLShell provider pipeline head-to-head (issue #28, task 30)\n";
        gsInfo << "===================================================================\n";
        gsInfo << " Material     : "<< matName <<"\n";
        gsInfo << "   -a 1 driver: "<< adapterName <<" + gsThinShellAssembler\n";
        gsInfo << "   -a 2 driver: gsShellMaterialProvider["
               << (material==0 ? "gsLinearMaterial" : "gsNeoHookeQuadMaterial")
               << "] + gsThinShellAssembler2\n";
        gsInfo << " Refinement   : "<< numRefine <<" (elevate "<<numElevate<<")\n";
        gsInfo << " #patches     : "<< mp.nPatches() <<"\n";
        gsInfo << " #elements    : "<< nElements <<"\n";
        gsInfo << " degree       : "<< dbasis.minCwiseDegree() <<"\n";
        gsInfo << " #DoFs        : "<< newAsm.numDofs()
               << " (adapter path: "<< adpAsm.numDofs() <<")\n";
        gsInfo << " reps         : "<< reps <<"\n";
        gsInfo << " NumGauss     : "<< ngProvider <<" (adapter "<< ngAdapter
               <<") -- ASSERTED equal\n";
        gsInfo << " OMP threads  : "<< nThreads
               << (exactCounters ? "  (counter identities ASSERTED)"
                                 : "  (counters printed, NOT asserted; timings NOT meaningful)")
               << "\n";
        gsInfo << "-------------------------------------------------------------------\n";

        //! [Agreement pass -- doubles as the UNTIMED warm-up of all three instances]
        // assembleMatrix goes first everywhere: assembleVector may NOT be the
        // first assembly on a fresh gsThinShellAssembler2 (initVector(1) does
        // not size the system matrix; task-28 finding, inherited from legacy).
        st = legAsm.assembleMatrix(mp_def);
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"legacy assembleMatrix() failed");
        gsSparseMatrix<real_t> K0 = legAsm.matrix();
        st = legAsm.assembleVector(mp_def);
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"legacy assembleVector() failed");
        gsMatrix<real_t> R0 = legAsm.rhs();

        st = adpAsm.assembleMatrix(mp_def);
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"adapter assembleMatrix() failed");
        gsSparseMatrix<real_t> K1 = adpAsm.matrix();
        st = adpAsm.assembleVector(mp_def);
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"adapter assembleVector() failed");
        gsMatrix<real_t> R1 = adpAsm.rhs();

        st = newAsm.assembleMatrix(mp_def);
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"provider assembleMatrix() failed");
        gsSparseMatrix<real_t> K2 = newAsm.matrix();   // copy OUT: the next
        st = newAsm.assembleVector(mp_def);            // assembly overwrites it
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"provider assembleVector() failed");
        gsMatrix<real_t> R2 = newAsm.rhs();

        const real_t K1n = (K1.norm()==0 ? 1.0 : K1.norm());
        const real_t R1n = (R1.norm()==0 ? 1.0 : R1.norm());
        const real_t K0n = (K0.norm()==0 ? 1.0 : K0.norm());
        const real_t R0n = (R0.norm()==0 ? 1.0 : R0.norm());
        const real_t dK21 = (K2-K1).norm() / K1n;
        const real_t dR21 = (R2-R1).norm() / R1n;
        const real_t dK20 = (K2-K0).norm() / K0n;
        const real_t dR20 = (R2-R0).norm() / R0n;
        const real_t agreeTol = 1e-12;

        gsInfo << " Cross-mode agreement (SAME deformed state, mp_def = 1.01*coefs).\n";
        gsInfo << " This pass is also the UNTIMED WARM-UP of all three instances.\n";
        gsInfo << "   vs -a 1 (SAME law, PARITY gate <= "<< agreeTol <<"):\n";
        gsInfo << "     ||K_2 - K_1|| / ||K_1|| = "<< dK21
               << (dK21==0.0 ? "   (EXACTLY 0)" : "")
               << "   ["<< (dK21<=agreeTol ? "OK" : "FAIL") <<"]\n";
        gsInfo << "     ||R_2 - R_1|| / ||R_1|| = "<< dR21
               << (dR21==0.0 ? "   (EXACTLY 0)" : "")
               << "   ["<< (dR21<=agreeTol ? "OK" : "FAIL") <<"]\n";
        gsInfo << "   vs -a 0 (legacy material -- printed, NOT gated):\n";
        gsInfo << "     ||K_2 - K_0|| / ||K_0|| = "<< dK20 <<"\n";
        gsInfo << "     ||R_2 - R_0|| / ||R_0|| = "<< dR20 <<"\n";
        if (material==0)
            gsInfo << "     [SvK: an O(t^2/R^2) MODEL difference (legacy integrates a "
                      "z-constant tangent\n      analytically, both new paths Gauss-integrate "
                      "the exact z-dependent metric);\n      ~4e-4 expected at t=0.25, R~25. NOT an error.]\n";
        else
            gsInfo << "     [NH : matched QuadraticVolumetric twin => PARITY expected here too.]\n";
        gsInfo << "-------------------------------------------------------------------\n";
        //! [Agreement pass]

        //! [Per-ASSEMBLY counter identities (untimed, one call each)]
        // The timed blocks below assert fills == reps*nElements over a whole
        // block; that aggregate would also be satisfied by a routine filling 2N
        // on one rep and 0 on the next. These three untimed single calls pin the
        // identity PER ASSEMBLY, which is the actual claim ("six coefficients ->
        // one sweep per element per assembly"). They also extend the warm-up.
        const size_t nE = (size_t)nElements;
        newAsm.resetMaterialFills();
        st = newAsm.assemble();
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"provider assemble() failed");
        const size_t f1A = newAsm.materialFills();
        newAsm.resetMaterialFills();
        st = newAsm.assembleMatrix(mp_def);
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"provider assembleMatrix() failed");
        const size_t f1M = newAsm.materialFills();
        const size_t m1M = newAsm.matrixMomentFills();
        const size_t s1M = newAsm.stressMomentFills();
        newAsm.resetMaterialFills();
        st = newAsm.assembleVector(mp_def);
        GISMO_ENSURE(st==ThinShellAssemblerStatus::Success,"provider assembleVector() failed");
        const size_t f1V = newAsm.materialFills();
        const size_t m1V = newAsm.matrixMomentFills();
        const size_t s1V = newAsm.stressMomentFills();

        gsInfo << " Provider counters for ONE assembly (#elements = "<< nElements <<"):\n";
        gsInfo << "   assemble()      : fills "<< f1A <<"\n";
        gsInfo << "   assembleMatrix  : fills "<< f1M <<"   matrixMoments "<< m1M
               << "   stressMoments "<< s1M <<"\n";
        gsInfo << "   assembleVector  : fills "<< f1V <<"   matrixMoments "<< m1V
               << "   stressMoments "<< s1V <<"\n";
        if (exactCounters)
        {
            GISMO_ENSURE(f1A==nE,"ONE assemble(): materialFills "<<f1A<<" != nElements "<<nE);
            GISMO_ENSURE(f1M==nE,"ONE assembleMatrix(): materialFills "<<f1M<<" != nElements "<<nE);
            GISMO_ENSURE(m1M==nE,"ONE assembleMatrix(): matrixMomentFills "<<m1M<<" != nElements "<<nE);
            GISMO_ENSURE(s1M==nE,"ONE assembleMatrix(): stressMomentFills "<<s1M<<" != nElements "<<nE);
            GISMO_ENSURE(f1V==nE,"ONE assembleVector(): materialFills "<<f1V<<" != nElements "<<nE);
            GISMO_ENSURE(m1V==0, "MASK GATE: ONE assembleVector() integrated tangent moments ("
                                 <<m1V<<" != 0)");
            GISMO_ENSURE(s1V==nE,"ONE assembleVector(): stressMomentFills "<<s1V<<" != nElements "<<nE);
            gsInfo << " [ASSERTED] fills == #elements for ONE assembly of each routine,\n";
            gsInfo << " [ASSERTED] and matrixMomentFills() == 0 for ONE assembleVector\n";
        }
        else
            gsInfo << " [not asserted: "<< nThreads <<" threads, the counters are non-atomic]\n";
        gsInfo << "-------------------------------------------------------------------\n";
        //! [Per-ASSEMBLY counter identities]

        //! [Timed blocks: reps INTERLEAVED, order ALTERNATED per rep]
        std::vector<real_t> tA1, tA2;   // assemble()
        std::vector<real_t> tM1, tM2;   // assembleMatrix(deformed)
        std::vector<real_t> tV1, tV2;   // assembleVector(deformed)

        // --- assemble()  (linear system; the provider path uses it for the
        //     Newton predictor and it also SIZES the system for assembleVector)
        newAsm.resetMaterialFills();
        gsMaterialMatrix3DResetSweeps();
        gsMaterialMatrix3DResetHits();
        gsTimeBlock(0,reps,sequential,adpAsm,newAsm,mp_def,tA1,tA2);
        const size_t fillsA  = newAsm.materialFills();
        const size_t sweepA  = gsMaterialMatrix3DSweeps();
        const size_t hitA    = gsMaterialMatrix3DHits();

        // --- assembleMatrix(deformed) : the tangent. Provider mask = all six
        //     moments (N and M enter the geometric terms), so matrixMomentFills
        //     and stressMomentFills both grow by nElements per rep.
        newAsm.resetMaterialFills();
        gsMaterialMatrix3DResetSweeps();
        gsMaterialMatrix3DResetHits();
        gsTimeBlock(1,reps,sequential,adpAsm,newAsm,mp_def,tM1,tM2);
        const size_t fillsM  = newAsm.materialFills();
        const size_t matMomM = newAsm.matrixMomentFills();
        const size_t strMomM = newAsm.stressMomentFills();
        const size_t sweepM  = gsMaterialMatrix3DSweeps();
        const size_t hitM    = gsMaterialMatrix3DHits();

        // --- assembleVector(deformed) : the residual, and the REQUIRED gate.
        //     Provider mask = N|M only: NO tangent moment is integrated at all,
        //     so matrixMomentFills() must stay at 0 over the whole block.
        newAsm.resetMaterialFills();
        gsMaterialMatrix3DResetSweeps();
        gsMaterialMatrix3DResetHits();
        gsTimeBlock(2,reps,sequential,adpAsm,newAsm,mp_def,tV1,tV2);
        const size_t fillsV  = newAsm.materialFills();
        const size_t matMomV = newAsm.matrixMomentFills();
        const size_t strMomV = newAsm.stressMomentFills();
        const size_t sweepV  = gsMaterialMatrix3DSweeps();
        const size_t hitV    = gsMaterialMatrix3DHits();
        //! [Timed blocks]

        const gsTimeStats sA1 = gsComputeTimeStats(tA1), sA2 = gsComputeTimeStats(tA2);
        const gsTimeStats sM1 = gsComputeTimeStats(tM1), sM2 = gsComputeTimeStats(tM2);
        const gsTimeStats sV1 = gsComputeTimeStats(tV1), sV2 = gsComputeTimeStats(tV2);

        //! [Report]
        gsInfo << " Wall times over "<< reps <<" reps, after an untimed warm-up assembly on\n";
        gsInfo << " every instance (the agreement pass above).\n";
        if (sequential)
            gsInfo << " Rep ordering: SEQUENTIAL (--sequential): all adapter reps, then all\n"
                      " provider reps. Bias control for the default interleaved ordering.\n\n";
        else
            gsInfo << " Rep ordering: INTERLEAVED and order-ALTERNATED (rep i even: adapter\n"
                      " first; rep i odd: provider first), so neither path profits from\n"
                      " running second. Run with --sequential for the bias control.\n\n";
        gsInfo << "Assembly kind       adapter (-a 1) [s]                provider (-a 2) [s]               speedup\n";
        gsInfo << "                  mean       min        max        mean       min        max        (a1/a2)\n";
        gsInfo << "----------------------------------------------------------------------------------------------\n";
        gsInfo << std::left  << std::setw(16) << "assemble()"    << std::right
               << std::setw(11) << sA1.mean << std::setw(11) << sA1.min << std::setw(11) << sA1.max
               << std::setw(13) << sA2.mean << std::setw(11) << sA2.min << std::setw(11) << sA2.max
               << std::setw(13) << sA1.mean/sA2.mean <<"\n";
        gsInfo << std::left  << std::setw(16) << "assembleMatrix" << std::right
               << std::setw(11) << sM1.mean << std::setw(11) << sM1.min << std::setw(11) << sM1.max
               << std::setw(13) << sM2.mean << std::setw(11) << sM2.min << std::setw(11) << sM2.max
               << std::setw(13) << sM1.mean/sM2.mean <<"\n";
        gsInfo << std::left  << std::setw(16) << "assembleVector" << std::right
               << std::setw(11) << sV1.mean << std::setw(11) << sV1.min << std::setw(11) << sV1.max
               << std::setw(13) << sV2.mean << std::setw(11) << sV2.min << std::setw(11) << sV2.max
               << std::setw(13) << sV1.mean/sV2.mean <<"\n";
        gsInfo << std::left  << std::setw(16) << "Newton step" << std::right
               << std::setw(11) << sM1.mean+sV1.mean << std::setw(11) << "" << std::setw(11) << ""
               << std::setw(13) << sM2.mean+sV2.mean << std::setw(11) << "" << std::setw(11) << ""
               << std::setw(13) << (sM1.mean+sV1.mean)/(sM2.mean+sV2.mean) <<"\n";
        gsInfo << "----------------------------------------------------------------------------------------------\n";
        gsInfo << " ('Newton step' = assembleMatrix + assembleVector, the per-iteration cost.)\n";
        gsInfo << "-------------------------------------------------------------------\n";

        const size_t expectFills = (size_t)reps * (size_t)nElements;
        gsInfo << " gsShellMaterialProvider counters (-a 2), over the whole timed block:\n";
        gsInfo << "   reps x #elements = "<< reps <<" x "<< nElements <<" = "<< expectFills <<"\n";
        gsInfo << "   assemble()      : fills "<< fillsA
               << "   (per assembly per element: "<< (real_t)fillsA/expectFills <<")\n";
        gsInfo << "   assembleMatrix  : fills "<< fillsM
               << "   matrixMoments "<< matMomM <<"   stressMoments "<< strMomM
               << "   (per assembly per element: "<< (real_t)fillsM/expectFills <<")\n";
        gsInfo << "   assembleVector  : fills "<< fillsV
               << "   matrixMoments "<< matMomV <<"   stressMoments "<< strMomV
               << "   (per assembly per element: "<< (real_t)fillsV/expectFills <<")\n";
        gsInfo << " => ONE sweep per element per assembly delivers all six moments,\n";
        gsInfo << "    and the RESIDUAL path integrates NO tangent moment at all.\n";
        if (exactCounters)
        {
            GISMO_ENSURE(fillsA==expectFills,
                "provider assemble(): materialFills "<<fillsA<<" != reps*nElements "<<expectFills);
            GISMO_ENSURE(fillsM==expectFills,
                "provider assembleMatrix(): materialFills "<<fillsM<<" != reps*nElements "<<expectFills);
            GISMO_ENSURE(matMomM==expectFills,
                "provider assembleMatrix(): matrixMomentFills "<<matMomM<<" != reps*nElements "<<expectFills);
            GISMO_ENSURE(strMomM==expectFills,
                "provider assembleMatrix(): stressMomentFills "<<strMomM<<" != reps*nElements "<<expectFills);
            GISMO_ENSURE(fillsV==expectFills,
                "provider assembleVector(): materialFills "<<fillsV<<" != reps*nElements "<<expectFills);
            GISMO_ENSURE(matMomV==0,
                "MASK GATE: provider assembleVector() integrated tangent moments ("
                <<matMomV<<" != 0)");
            GISMO_ENSURE(strMomV==expectFills,
                "provider assembleVector(): stressMomentFills "<<strMomV<<" != reps*nElements "<<expectFills);
            gsInfo << " [ASSERTED] fills == reps*#elements for all three routines\n";
            gsInfo << " [ASSERTED] matrixMomentFills() == 0 after assembleVector  <-- MASK GATE\n";
        }
        else
            gsInfo << " [not asserted: "<< nThreads <<" threads, the counters are non-atomic]\n";

        gsInfo << "-------------------------------------------------------------------\n";
        gsInfo << " gsMaterialMatrix3D adapter cross-output cache (-a 1), same blocks:\n";
        gsInfo << "   assemble()      : sweeps "<< sweepA <<" hits "<< hitA
               << "   (per assembly per element: sweeps "<< (real_t)sweepA/expectFills <<")\n";
        gsInfo << "   assembleMatrix  : sweeps "<< sweepM <<" hits "<< hitM
               << "   (per assembly per element: sweeps "<< (real_t)sweepM/expectFills <<")\n";
        gsInfo << "   assembleVector  : sweeps "<< sweepV <<" hits "<< hitV
               << "   (per assembly per element: sweeps "<< (real_t)sweepV/expectFills <<")\n";
        gsInfo << " => the adapter ALSO sweeps once per element per assembly kind, but its\n";
        gsInfo << "    residual sweep is a FULL one (setDeformed bumps the config revision,\n";
        gsInfo << "    so the tangent it computed for assembleMatrix cannot be reused and\n";
        gsInfo << "    the plane-stress condensation still produces the condensed tangent).\n";
        gsInfo << "-------------------------------------------------------------------\n";

        const bool gateVector = (sV2.mean < sV1.mean);
        const bool gateMatrix = (sM2.mean <= sM1.mean);
        gsInfo << " VERDICT ["<< matName <<"]\n";
        gsInfo << "   REQUIRED  assembleVector : provider "<< sV2.mean <<" s  vs  adapter "
               << sV1.mean <<" s  ->  "
               << (gateVector ? "PROVIDER FASTER" : "PROVIDER SLOWER")
               << "  ("<< sV1.mean/sV2.mean <<"x)  ["<< (gateVector ? "PASS" : "FAIL") <<"]\n";
        gsInfo << "   expected  assembleMatrix : provider "<< sM2.mean <<" s  vs  adapter "
               << sM1.mean <<" s  ->  "
               << (gateMatrix ? "PROVIDER <= ADAPTER" : "PROVIDER SLOWER")
               << "  ("<< sM1.mean/sM2.mean <<"x)  ["<< (gateMatrix ? "PASS" : "FAIL") <<"]\n";
        if (!exactCounters)
            gsInfo << "   [WARNING: "<< nThreads <<" threads -- these timings are an agreement/"
                      "stability check only, NOT a measurement.]\n";
        gsInfo << "===================================================================\n";
        //! [Report]

        return EXIT_SUCCESS;
    }
    //! [Provider pipeline head-to-head]
#endif

    //! [Assembler]
    gsThinShellAssemblerBase<real_t>* assembler =
        new gsThinShellAssembler<3, real_t, true>(mp,dbasis,bc,force,activeMaterial);
    //! [Assembler]

    gsInfo << "===================================================================\n";
    gsInfo << " gsKLShell material redundancy baseline (issue #28)\n";
    gsInfo << "===================================================================\n";
    gsInfo << " Material     : "<< activeName <<"\n";
    gsInfo << " Refinement   : "<< numRefine <<" (elevate "<<numElevate<<")\n";
    gsInfo << " #patches     : "<< mp.nPatches() <<"\n";
    gsInfo << " #elements    : "<< nElements <<"\n";
    gsInfo << " degree       : "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << " #DoFs        : "<< assembler->numDofs() <<"\n";
    gsInfo << " reps         : "<< reps <<"\n";
    gsInfo << "-------------------------------------------------------------------\n";

    gsStopwatch clock;
    ThinShellAssemblerStatus status;

#ifdef gsPhaseFieldFracture_ENABLED
    //! [Cross-mode agreement check (-a 1 only, printed once before timing)]
    if (useAdapter)
    {
        gsThinShellAssemblerBase<real_t>* legAsm =
            new gsThinShellAssembler<3, real_t, true>(mp,dbasis,bc,force,legacyMaterial.get());
        gsThinShellAssemblerBase<real_t>* adpAsm =
            new gsThinShellAssembler<3, real_t, true>(mp,dbasis,bc,force,adapterMaterial.get());

        status = legAsm->assembleMatrix(mp_def);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"legacy assembleMatrix() failed");
        gsSparseMatrix<real_t> Kleg = legAsm->matrix();
        status = legAsm->assembleVector(mp_def);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"legacy assembleVector() failed");
        gsMatrix<real_t> Rleg = legAsm->rhs();

        status = adpAsm->assembleMatrix(mp_def);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"adapter assembleMatrix() failed");
        gsSparseMatrix<real_t> Kadp = adpAsm->matrix();
        status = adpAsm->assembleVector(mp_def);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"adapter assembleVector() failed");
        gsMatrix<real_t> Radp = adpAsm->rhs();

        const real_t Knorm = Kleg.norm();
        const real_t Rnorm = Rleg.norm();
        const real_t Kdiff = (Kleg - Kadp).norm() / (Knorm==0 ? 1.0 : Knorm);
        const real_t Rdiff = (Rleg - Radp).norm() / (Rnorm==0 ? 1.0 : Rnorm);

        gsInfo << " Cross-mode agreement (legacy vs adapter, SAME deformed state):\n";
        gsInfo << "   ||K_legacy - K_adapter|| / ||K_legacy|| = "<< Kdiff <<"\n";
        gsInfo << "   ||R_legacy - R_adapter|| / ||R_legacy|| = "<< Rdiff <<"\n";
        if (material==0)
            gsInfo << "   [SvK: O(t^2/R^2) MODEL difference expected ~1e-4 here "
                      "(t=0.25,R~25); task-15 measured 8.4e-6 at t=0.01. NOT an error.]\n";
        else
            gsInfo << "   [NH : matched QuadraticVolumetric twin => PARITY expected <= ~1e-8.]\n";
        gsInfo << "-------------------------------------------------------------------\n";

        delete legAsm;
        delete adpAsm;
    }
    //! [Cross-mode agreement check]
#endif

    //! [Timing: assemble() (linear / first assembly)]
    gsMaterialMatrixResetComputePointsCalls();
    gsMaterialMatrixResetComputePointsMisses();
#ifdef gsPhaseFieldFracture_ENABLED
    size_t sw_assemble = 0, hit_assemble = 0;
    if (useAdapter) { gsMaterialMatrix3DResetSweeps(); gsMaterialMatrix3DResetHits(); }
#endif
    clock.restart();
    status = assembler->assemble();
    real_t t_assemble = clock.stop();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"assemble() failed");
    const size_t c_assemble = gsMaterialMatrixComputePointsCalls();
    const size_t m_assemble = gsMaterialMatrixComputePointsMisses();
#ifdef gsPhaseFieldFracture_ENABLED
    if (useAdapter) { sw_assemble = gsMaterialMatrix3DSweeps(); hit_assemble = gsMaterialMatrix3DHits(); }
#endif
    //! [Timing: assemble() (linear / first assembly)]

    //! [Timing: assembleMatrix(deformed)]
    gsMaterialMatrixResetComputePointsCalls();
    gsMaterialMatrixResetComputePointsMisses();
#ifdef gsPhaseFieldFracture_ENABLED
    size_t sw_matrix_total = 0, hit_matrix_total = 0;
    if (useAdapter) { gsMaterialMatrix3DResetSweeps(); gsMaterialMatrix3DResetHits(); }
#endif
    real_t t_matrix = 0.0;
    for (index_t i=0; i!=reps; ++i)
    {
        clock.restart();
        status = assembler->assembleMatrix(mp_def);
        t_matrix += clock.stop();
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"assembleMatrix() failed");
    }
    const size_t c_matrix_total = gsMaterialMatrixComputePointsCalls();
    const size_t m_matrix_total = gsMaterialMatrixComputePointsMisses();
#ifdef gsPhaseFieldFracture_ENABLED
    if (useAdapter) { sw_matrix_total = gsMaterialMatrix3DSweeps(); hit_matrix_total = gsMaterialMatrix3DHits(); }
#endif
    t_matrix /= reps;
    const real_t c_matrix = (real_t)c_matrix_total / reps;
    const real_t m_matrix = (real_t)m_matrix_total / reps;
    //! [Timing: assembleMatrix(deformed)]

    //! [Timing: assembleVector(deformed)]
    gsMaterialMatrixResetComputePointsCalls();
    gsMaterialMatrixResetComputePointsMisses();
#ifdef gsPhaseFieldFracture_ENABLED
    size_t sw_vector_total = 0, hit_vector_total = 0;
    if (useAdapter) { gsMaterialMatrix3DResetSweeps(); gsMaterialMatrix3DResetHits(); }
#endif
    real_t t_vector = 0.0;
    for (index_t i=0; i!=reps; ++i)
    {
        clock.restart();
        status = assembler->assembleVector(mp_def);
        t_vector += clock.stop();
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"assembleVector() failed");
    }
    const size_t c_vector_total = gsMaterialMatrixComputePointsCalls();
    const size_t m_vector_total = gsMaterialMatrixComputePointsMisses();
#ifdef gsPhaseFieldFracture_ENABLED
    if (useAdapter) { sw_vector_total = gsMaterialMatrix3DSweeps(); hit_vector_total = gsMaterialMatrix3DHits(); }
#endif
    t_vector /= reps;
    const real_t c_vector = (real_t)c_vector_total / reps;
    const real_t m_vector = (real_t)m_vector_total / reps;
    //! [Timing: assembleVector(deformed)]

    //! [Report]
    gsInfo << "Assembly kind      mean time [s]     calls/call  misses/call   miss/elem\n";
    gsInfo << "-------------------------------------------------------------------------\n";
    gsInfo << std::left  << std::setw(18) << "assemble()"
           << std::right << std::setw(12) << t_assemble
           << std::setw(14) << c_assemble
           << std::setw(13) << m_assemble
           << std::setw(12) << (real_t)m_assemble/nElements <<"\n";
    gsInfo << std::left  << std::setw(18) << "assembleMatrix"
           << std::right << std::setw(12) << t_matrix
           << std::setw(14) << c_matrix
           << std::setw(13) << m_matrix
           << std::setw(12) << m_matrix/nElements <<"\n";
    gsInfo << std::left  << std::setw(18) << "assembleVector"
           << std::right << std::setw(12) << t_vector
           << std::setw(14) << c_vector
           << std::setw(13) << m_vector
           << std::setw(12) << m_vector/nElements <<"\n";
    gsInfo << "-------------------------------------------------------------------------\n";
    gsInfo << " _computePoints per element (calls = invocations, misses = actual recomputes):\n";
    gsInfo << "   assembleMatrix : calls "<< c_matrix/nElements
           << " -> misses "<< m_matrix/nElements <<"  (guard target: ~1)\n";
    gsInfo << "   assembleVector : calls "<< c_vector/nElements
           << " -> misses "<< m_vector/nElements <<"  (guard target: ~1)\n";
    gsInfo << "   assemble()     : calls "<< (real_t)c_assemble/nElements
           << " -> misses "<< (real_t)m_assemble/nElements <<"\n";
    gsInfo << " => a full Newton assembly (matrix+vector) INVOKES points "
           << c_matrix/nElements + c_vector/nElements
           << "x per element but with the same-input guard RECOMPUTES only "
           << m_matrix/nElements + m_vector/nElements <<"x.\n";
#ifdef gsPhaseFieldFracture_ENABLED
    if (useAdapter)
    {
        const real_t sw_matrix  = (real_t)sw_matrix_total  / reps;
        const real_t hit_matrix = (real_t)hit_matrix_total / reps;
        const real_t sw_vector  = (real_t)sw_vector_total  / reps;
        const real_t hit_vector = (real_t)hit_vector_total / reps;
        gsInfo << "-------------------------------------------------------------------------\n";
        gsInfo << " gsMaterialMatrix3D adapter cross-output cache (sweeps=misses, hits=reuse):\n";
        gsInfo << "   assemble()     : sweeps "<< sw_assemble
               << " hits "<< hit_assemble
               << "  (per elem: sweeps "<< (real_t)sw_assemble/nElements
               << " hits "<< (real_t)hit_assemble/nElements <<")\n";
        gsInfo << "   assembleMatrix : sweeps "<< sw_matrix
               << " hits "<< hit_matrix
               << "  (per elem: sweeps "<< sw_matrix/nElements
               << " hits "<< hit_matrix/nElements <<")   [A,B,C,D share 1 sweep]\n";
        gsInfo << "   assembleVector : sweeps "<< sw_vector
               << " hits "<< hit_vector
               << "  (per elem: sweeps "<< sw_vector/nElements
               << " hits "<< hit_vector/nElements <<")   [N,M share 1 sweep]\n";
        gsInfo << " => the six A/B/C/D/N/M integrators collapse onto 1 constitutive sweep\n"
               << "    per element per assembly kind; the rest are cache hits.\n";
    }
#endif
    gsInfo << "===================================================================\n";
    //! [Report]

    delete assembler;
    return EXIT_SUCCESS;
}
