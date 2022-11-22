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
    typedef typename gsThinShellAssemblerDWRBase<T>::bContainer bContainer;


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
        gsMatrix<T> points;
        bContainer bnds;
        std::string empty;
        this->computeError(deformed,primalL,bnds,points,true,empty);
    }

    void computeError(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> primalL, const bContainer & bnds, const gsMatrix<T> & points, bool interior)
    {
        std::string empty;
        this->computeError(deformed,primalL,bnds,points,interior,empty);
    }


    void computeError(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> primalL,
                        std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false)
    {
        gsMatrix<T> points;
        bContainer bnds;
        this->computeError(deformed,primalL,bnds,points,true,filename,np,parametric,mesh);
    }

    void computeError(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> primalL,
                        const bContainer & bnds, const gsMatrix<T> & points, bool interior,
                        std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false)
    {
        gsMultiPatch<T> dualL, dualH;
        gsVector<T> solVector;
        gsVector<T> rhsL(m_assembler->numDofsL());
        rhsL.setZero();
        gsVector<T> rhsH(m_assembler->numDofsH());
        rhsH.setZero();

        gsSparseSolver<real_t>::uPtr solver;
        #ifdef GISMO_WITH_PARDISO
            solver = gsSparseSolver<real_t>::get( "PardisoLU");
        #else
            solver = gsSparseSolver<real_t>::get( "LU");
        #endif

        gsInfo << "Assembling dual matrix (L)... "<< std::flush;
        m_assembler->assembleMatrixL(deformed);
        gsInfo << "done.\n";
        gsInfo << "Assembling dual vector (L)... "<< std::flush;
        if (interior)
        {
            m_assembler->assembleDualL(primalL,deformed);
            rhsL+=m_assembler->dualL();
        }
        if (bnds.size()!=0)
        {
            m_assembler->assembleDualL(bnds,primalL,deformed);
            rhsL+=m_assembler->dualL();
        }
        if (points.cols()!=0)
        {
            m_assembler->assembleDualL(points,primalL,deformed);
            rhsL+=m_assembler->dualL();
        }
        gsInfo << "done.\n";
        gsInfo << "Solving dual (L), size = "<<m_assembler->matrixL().rows()<<","<<m_assembler->matrixL().cols()<<"... "<< std::flush;
        solver->compute(m_assembler->matrixL());
        solVector = solver->solve(rhsL);
        m_assembler->constructMultiPatchL(solVector,dualL);
        gsInfo << "done.\n";

        gsInfo << "Assembling dual matrix (H)... "<< std::flush;
        m_assembler->assembleMatrixH(deformed);
        gsInfo << "done.\n";
        gsInfo << "Assembling dual vector (H)... "<< std::flush;
        if (interior)
        {
            m_assembler->assembleDualH(primalL,deformed);
            rhsH+=m_assembler->dualH();
        }
        if (bnds.size()!=0)
        {
            m_assembler->assembleDualH(bnds,primalL,deformed);
            rhsH+=m_assembler->dualH();
        }
        if (points.cols()!=0)
        {
            m_assembler->assembleDualH(points,primalL,deformed);
            rhsH+=m_assembler->dualH();
        }

        gsInfo << "done.\n";
        gsInfo << "Solving dual (H), size = "<<m_assembler->matrixH().rows()<<","<<m_assembler->matrixH().cols()<<"... "<< std::flush;
        solver->compute(m_assembler->matrixH());
        solVector = solver->solve(rhsH);
        m_assembler->constructMultiPatchH(solVector,dualH);
        gsInfo << "done.\n";

        m_errors = m_assembler->computeErrorElements(dualL,dualH,deformed,filename,np,parametric,mesh);
        m_error = m_assembler->error();
    }

    T error() const { return  m_error; }
    std::vector<T> errors() const { return m_errors; }
    std::vector<T> absErrors() const
    { 
        std::vector<T> result = m_errors;
        for (typename std::vector<T>::iterator it = result.begin(); it!=result.end(); it++)
            *it = std::abs(*it);
        return result;
    }

protected:
    // bool m_verbose;
    gsThinShellAssemblerDWRBase<T> * m_assembler;
    T m_error;
    std::vector<T> m_errors;
};

} // namespace gismo
