/** @file gsThinShellDWRHelper.h

    @brief Provides DWR assembly routines for the Kirchhoff-Love shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsThinShellAssemblerDWR.h>

namespace gismo
{

template <class T>
class gsThinShellDWRHelper

{
public:
    // empty constructor
    gsThinShellDWRHelper() {};


    gsThinShellDWRHelper(gsThinShellAssemblerDWRBase<T> * assembler)
    :
    m_assembler(assembler)
    {

    }

    void computeError(const gsVector<T> & U)
    {
        gsMultiPatch<T> primalL, deformed;

        m_assembler->constructSolutionL(U,deformed);
        m_assembler->constructMultiPatchL(U,primalL);

        this->computeError(deformed,primalL);
    }

    void computeError(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> primalL)
    {
        std::string empty;
        this->computeError(deformed,primalL,empty);
    }

    void computeError(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> primalL,
                        std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false)
    {
        gsMultiPatch<T> dualL, dualH;
        gsVector<T> solVector;
        gsSparseSolver<>::LU solver;

        gsInfo << "Assembling dual matrix (L)... "<< std::flush;
        m_assembler->assembleMatrixL(deformed);
        gsInfo << "done.\n";
        gsInfo << "Assembling dual vector (L)... "<< std::flush;
        m_assembler->assembleDualL(primalL,deformed);
        gsInfo << "done.\n";
        gsInfo << "Solving dual (L), size = "<<m_assembler->matrixL().rows()<<","<<m_assembler->matrixL().cols()<<"... "<< std::flush;
        solver.compute(m_assembler->matrixL());
        solVector = solver.solve(m_assembler->dualL());
        m_assembler->constructMultiPatchL(solVector,dualL);
        gsInfo << "done.\n";

        gsInfo << "Assembling dual matrix (H)... "<< std::flush;
        m_assembler->assembleMatrixH(deformed);
        gsInfo << "done.\n";
        gsInfo << "Assembling dual vector (H)... "<< std::flush;
        m_assembler->assembleDualH(primalL,deformed);
        gsInfo << "done.\n";
        gsInfo << "Solving dual (H), size = "<<m_assembler->matrixH().rows()<<","<<m_assembler->matrixH().cols()<<"... "<< std::flush;
        solver.compute(m_assembler->matrixH());
        solVector = solver.solve(m_assembler->dualH());
        m_assembler->constructMultiPatchH(solVector,dualH);
        gsInfo << "done.\n";

        m_assembler->computeErrorElements(dualL,dualH,deformed,filename,np,parametric,mesh);
    }

    T error() { return m_assembler->error(); }
    std::vector<T> absErrors() { return m_assembler->absErrors(); }

protected:
    // bool m_verbose;
    gsThinShellAssemblerDWRBase<T> * m_assembler;
};

} // namespace gismo
