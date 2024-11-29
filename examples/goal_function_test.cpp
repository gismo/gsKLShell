/** @file commandLineArg_example.cpp

    @brief Tutorial on how to use command line parser in G+Smo.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <gismo.h>
#include <variant>
#include <gsKLShell/src/gsGoalFunctionExpr.h>
#include <gsKLShell/src/getMaterialMatrix.h>

using namespace gismo;

template <class _GF>
void assemble(gsExprAssembler<> & A, _GF gf)
{
  typedef gsExprAssembler<>::space space;
  space u = A.trialSpace(0);
  A.initSystem();
  A.assemble(gf.derivative(u));
  gsDebugVar(A.rhs());
}

// // Static dispatch function
// template<class _GF>
// void runMe(_GF & gf)
// {

// }

// Preferred
// Static dispatch function
template<class T, enum GoalFunction _GF>
void runMe(gsGoalFunction<T,_GF> & gf, gsExprAssembler<T> & A)
{
  gsExprEvaluator<> ev(A);
  gsExprAssembler<>::space u = A.trialSpace(0);

  gf.setAssembler(&A);

  A.initSystem();
  A.assemble(gf.derivative(u));
  gsDebugVar(A.rhs());
}


// template<class T>
// class gsMyGeometry : public gsControlledFunction<T>
// {

//   // G = [a*x, b*y, x^2+y^2]

//   gsMyGeometry()
//   {}

//   void eval_into()
//   {
//     mp.eval_into();
//   }

//   gsVector<> parameters()
//   {
//     mp.coefs()[slice]
//   }

// }

int main(int argc, char* argv[])
{
  index_t goal = 0;
  index_t comp = 9;
  gsCmdLine cmd("Tutorial Command Line Arguments");
  cmd.addInt("g", "goal", "goal function: 0 = , 1 = ",goal);
  cmd.addInt("c", "comp", "goal function component",comp);
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


  gsMultiPatch<> mp_ori;
  mp_ori.addPatch(gsNurbsCreator<>::BSplineSquare());
  mp_ori.embed(3);
  gsMultiPatch<> mp_def = mp_ori;
  mp_def.patch(0).coefs().col(1).array() += 1;
  gsMultiBasis<> mb(mp_ori);

  typedef gsExprAssembler<>::space space;
  gsExprAssembler<> A(1,1);
  A.setIntegrationElements(mb);
  A.getSpace(mb,3);

  // Make material matrix
  std::vector<gsFunctionSet<>*> parameters;
  gsMaterialMatrixBase<real_t>* materialMatrix;
  gsOptionList options;

  gsFunctionExpr<> E(std::to_string(1.0),3);
  gsFunctionExpr<> t(std::to_string(1.0),3);
  gsFunctionExpr<> nu(std::to_string(0.3),3);

  parameters.resize(2);
  parameters[0] = &E;
  parameters[1] = &nu;
  options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
  options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
  materialMatrix = getMaterialMatrix<3,real_t>(mp_ori,t,parameters,options);

  // GoalFunctionType gf;
//     gsInfo<<"Using GoalFunction::DisplacementNorm\n";
  // gf =

  // Static dispatch
  if      (goal == 0 && comp != 9)
  {
    gsGoalFunction<real_t,GoalFunction::DisplacementComponent> gf(mp_ori,mp_def,comp);
    runMe(gf,A);
  }
  else if (goal == 0 && comp == 9)
  {
    gsGoalFunction<real_t,GoalFunction::DisplacementNorm> gf(mp_ori,mp_def);
    runMe(gf,A);
  }
  else if (goal == 1 && comp == 9)
  {
    gsGoalFunction<real_t,GoalFunction::StretchNorm> gf(mp_ori,mp_def,materialMatrix);
    runMe(gf,A);
  }

  else
  {}

  return 0;
}
// //
//     // auto gf = (goal==0 && comp==9) ? gsGoalFunction<real_t,GoalFunction::DisplacementNorm>(A,mp_ori,mp_def)
//     //                                : gsGoalFunction<real_t,GoalFunction::DisplacementComponent>(A,mp_ori,mp_def,comp);

//     // gsGoalFunctionBase<real_t> *  gf1 = new gsGoalFunction<real_t,GoalFunction::DisplacementNorm>(A,mp_ori,mp_def);
//     // gsGoalFunctionBase<real_t> *  gf2 = new gsGoalFunction<real_t,GoalFunction::DisplacementComponent>(A,mp_ori,mp_def,comp);

//     // gsGoalFunctionBase<real_t> * gf;
//     // if (goal==0)
//     // {
//     //   if (comp==9)
//     //     gf = new gsGoalFunction<real_t,GoalFunction::DisplacementNorm>(A,mp_ori,mp_def);
//     //   else
//     //     gf = new gsGoalFunction<real_t,GoalFunction::DisplacementComponent>(A,mp_ori,mp_def,comp);
//     // }
//     // else if (goal ==1)
//     //     gf = new gsGoalFunction<real_t,GoalFunction::DisplacementComponent>(A,mp_ori,mp_def);

//     gsDebugVar(gf.function());
//     // gsDebugVar(gf->function2());
//     gsDebugVar(gf.derivative(u));


//     gsDebugVar(ev.integral(gf.function()));

//     assemble(A,gf);

    // delete gf;
