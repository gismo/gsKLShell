/** @file gsCompositeMatrix_test.cpp

    @brief Provides unittests for gsCompositeMatrix, the plane-stress stiffness
           of an orthotropic ply.

    gsCompositeMatrix(Exx,Eyy,Gxy,nuxy,nuyx) requires the reciprocity relation
    nuyx*Exx = nuxy*Eyy up to a relative mismatch of 1%. Within the tolerance
    it returns an exactly symmetric matrix; beyond it, it throws
    std::runtime_error (GISMO_ENSURE) in every build type.

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_MATRIX_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

SUITE(gsCompositeMatrix_test)                 // The suite should have the same name as the file
{
    // Isotropic limit: Exx = Eyy = E, nuxy = nuyx = nu, Gxy = E/(2(1+nu)) gives the
    // plane-stress matrix E/(1-nu^2) [1 nu 0; nu 1 0; 0 0 (1-nu)/2].
    TEST(isotropic)
    {
        const real_t E = 210.0, nu = 0.3;
        gsMatrix<real_t> G = gsCompositeMatrix<real_t>(E,E,E/(2*(1+nu)),nu,nu);

        gsMatrix<real_t> ref(3,3);
        ref << 1,  nu, 0,
               nu, 1,  0,
               0,  0,  (1-nu)/2;
        ref *= E/(1-nu*nu);
        CHECK_MATRIX_CLOSE(ref,G,1e-12*E);
    }

    // Consistent orthotropic input, nu21 = nu12*E2/E1: G(0,1) = G(1,0) = nu12*E2/(1-nu12*nu21).
    TEST(orthotropic_consistent)
    {
        const real_t E1 = 140.0, E2 = 10.0, G12 = 5.0, nu12 = 0.3;
        const real_t nu21 = nu12*E2/E1;
        gsMatrix<real_t> G = gsCompositeMatrix<real_t>(E1,E2,G12,nu12,nu21);

        const real_t den = 1 - nu12*nu21;
        CHECK_CLOSE(E1/den,      G(0,0), 1e-12*E1);
        CHECK_CLOSE(E2/den,      G(1,1), 1e-12*E1);
        CHECK_CLOSE(G12,         G(2,2), 1e-12*E1);
        CHECK_CLOSE(nu12*E2/den, G(0,1), 1e-12*E1);
        CHECK_EQUAL(G(0,1), G(1,0));
        CHECK_EQUAL(0.0, G(0,2)); CHECK_EQUAL(0.0, G(2,0));
        CHECK_EQUAL(0.0, G(1,2)); CHECK_EQUAL(0.0, G(2,1));
    }

    // Datasheet-rounded nu21 = 0.0214 (exact: 0.0214285...): nu21*E1 = 2.996 vs nu12*E2 = 3,
    // a relative mismatch of 1.3e-3, is accepted and the off-diagonals become their mean.
    TEST(orthotropic_rounded)
    {
        const real_t E1 = 140.0, E2 = 10.0, G12 = 5.0, nu12 = 0.3, nu21 = 0.0214;
        gsMatrix<real_t> G = gsCompositeMatrix<real_t>(E1,E2,G12,nu12,nu21);
        CHECK_EQUAL(G(0,1), G(1,0));
        CHECK_CLOSE(0.5*(nu21*E1 + nu12*E2)/(1 - nu12*nu21), G(0,1), 1e-12*E1);
    }

    // Equal Poisson ratios with Exx != Eyy: nu21*E1 = 42 vs nu12*E2 = 3, rejected in every build type.
    TEST(orthotropic_inconsistent)
    {
        CHECK_THROW(gsCompositeMatrix<real_t>(140.0,10.0,5.0,0.3,0.3), std::runtime_error);
    }
}
