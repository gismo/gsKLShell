/** @file gsFunctionSum.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
*/

#include<gsCore/gsFunctionSet.h>

#pragma once

namespace gismo
{
template<class T> class gsFunctionPieceSum;

template<class T>
class gsFunctionSum : public gsFunctionSet<T>
{
public:
    gsFunctionSum(const std::vector<const gsFunctionSet<T> * > & functions)
    :
    m_functions(functions),
    m_size(functions.size())
    {
        _init(m_functions);
    }

    gsFunctionSum(const gsFunctionSet<T> * undeformed, const gsFunctionSet<T> * deformation)
    :
    m_size(2)
    {
        m_functions.push_back(undeformed);
        m_functions.push_back(deformation);

        _init(m_functions);
    }

    index_t nPieces() const {return m_size;}

    short_t domainDim() const
    { return m_functions[0]->domainDim(); }

    short_t targetDim() const
    { return m_functions[0]->targetDim(); }

    const gsFunctionPieceSum<T> & piece(const index_t k) const { return m_pieces[k]; }

    index_t size() const {return nPieces();}

    std::ostream &print(std::ostream &os) const
    {
        for (index_t p = 1; p!=m_size; p++)
            gsInfo<<"Piece "<<p<<":\n"<<m_pieces[p]<<"\n";
    }


    index_t funcSize() const {return m_functions.size();}

    const gsFunctionSet<T> * func(const index_t i) const { return m_functions[i]; }

private:
    void _init(const std::vector<const gsFunctionSet<T> * > & functions)
    {
        GISMO_ASSERT(functions.size()>0,"Must give one or more function sets.");
        for (index_t p = 1; p!=m_size; p++)
        {
            GISMO_ASSERT(functions[p-1]->nPieces()==functions[p]->nPieces(),"Number of pieces does not match for function "<<p-1<<" and "<<p<<"!");
            GISMO_ASSERT(functions[p-1]->domainDim()==functions[p]->domainDim(),"Domain dimension does not match for function "<<p-1<<" and "<<p<<"!");
            GISMO_ASSERT(functions[p-1]->targetDim()==functions[p]->targetDim(),"Target dimension does not match for function "<<p-1<<" and "<<p<<"!");
        }

        for (index_t p = 0; p!=m_size; p++)
            m_pieces.push_back(gsFunctionPieceSum<T>(this,p));
    }

protected:
    std::vector<const gsFunctionSet<T> * > m_functions;
    std::vector<gsFunctionPieceSum<T> > m_pieces;
    const index_t m_size;

};

template<class T>
class gsFunctionPieceSum :  public gsFunction<T>
{
    /// Shared pointer for gsFunctionPieceSum
    typedef memory::shared_ptr< gsFunctionPieceSum<T> > Ptr;

    /// Auto pointer for gsFunctionExpr
    typedef memory::unique_ptr< gsFunctionPieceSum<T> > uPtr;

public:
    gsFunctionPieceSum(const gsFunctionSum<T> * geom, const index_t index)
    :
    m_geom(geom),
    m_index(index)
    {

    }

    short_t domainDim() const
    { return m_geom->domainDim(); }

    short_t targetDim() const
    { return m_geom->targetDim(); }

    /// Evaluates the non-zero spline functions at value u.
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_geom->func(0)->piece(m_index).eval_into(u,result);
        for (index_t p = 1; p!=m_geom->funcSize(); p++)
        {
            m_geom->func(p)->piece(m_index).eval_into(u,tmp);
            result += tmp;
        }
    }

    /// Evaluates the (partial) derivatives of non-zero spline functions at (the columns of) u.
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        gsMatrix<T> tmp;
        m_geom->func(0)->piece(m_index).deriv_into(u,result);
        for (index_t p = 1; p!=m_geom->funcSize(); p++)
        {
            m_geom->func(p)->piece(m_index).deriv_into(u,tmp);
            result += tmp;
        }

    }

    /// Evaluates the (partial) derivatives of the nonzero spline functions at points \a u into \a result.
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        gsMatrix<T> tmp;
        m_geom->func(0)->piece(m_index).deriv2_into(u,result);
        for (index_t p = 1; p!=m_geom->funcSize(); p++)
        {
            m_geom->func(p)->piece(m_index).deriv2_into(u,tmp);
            result += tmp;
        }

    }

    /// @brief Evaluate the nonzero spline functions and their derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const gsMatrix<T> & u, int n, std::vector<gsMatrix<T> >& result) const
    {
        std::vector<gsMatrix<T> > tmp;
        m_geom->func(0)->piece(m_index).evalAllDers_into(u,n,result);
        for (index_t p = 1; p!=m_geom->funcSize(); p++)
        {
            m_geom->func(p)->piece(m_index).evalAllDers_into(u,n,tmp);
            for (index_t k = 0; result.size(); k++)
            {
                result[k] += tmp[k];
            }
        }
    }

    GISMO_CLONE_FUNCTION(gsFunctionPieceSum)

    std::ostream &print(std::ostream &os) const
    {
        for (index_t p = 0; p!= m_geom->funcSize(); p++)
            gsInfo<<" function "<<p<<":\n"<<m_geom->func(p)->piece(m_index)<<"\n";
    }

protected:
    const gsFunctionSum<T> * m_geom;
    const index_t m_index;

};

} // namespace gismo
