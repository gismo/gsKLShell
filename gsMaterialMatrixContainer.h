/** @file gsMaterialMatrixContainer.h

    @brief Provides a container for material matrices for thin shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{


/**
 * @brief      This class serves as the evaluator of material matrices, based on \ref gsMaterialMatrixBase
 *
 * @tparam     T     Real tyoe
 * @tparam     out   Output type (see \ref MaterialOutput)
 *
 * @ingroup    KLShell
 */
template <class T>
class gsMaterialMatrixContainer
{
public:
    typedef typename std::vector<gsMaterialMatrixBase<T> *> Container;

    typedef typename Container::iterator       iterator;
    typedef typename Container::const_iterator const_iterator;

    /// Shared pointer for gsMaterialMatrixContainer
    typedef memory::shared_ptr< gsMaterialMatrixContainer > Ptr;

    /// Unique pointer for gsMaterialMatrixContainer
    typedef memory::unique_ptr< gsMaterialMatrixContainer > uPtr;
public:

    /// Constructor
    gsMaterialMatrixContainer( index_t size = 0 )
    :
    m_size(size)
    {
        m_container.reserve(m_size);
    }

    gsMaterialMatrixContainer(gsMaterialMatrixBase<T> * mat)
    {
        m_container.push_back(mat);
    }

    gsMaterialMatrixContainer(const gsMaterialMatrixContainer & other)
    {
        m_container.resize(other.m_container.size() );
        cloneAll( other.m_container.begin(), other.m_container.end(),
                  m_container.begin() );
    }

    gsMaterialMatrixContainer(Container & funcs)
    {
        m_container.swap(m_container); // funcs are consumed
    }

    ~gsMaterialMatrixContainer()
    { }

    void add(gsMaterialMatrixBase<T> * mat)
    {
        m_container.push_back( mat );
    }

    gsMaterialMatrixBase<T> * piece(const index_t k) const
    {
        return m_container.at(k);
    }

    index_t size() const {return m_container.size();}

    std::ostream &print(std::ostream &os) const
    {
        os << "Piecewise Function with "<<m_container.size() <<" pieces.\n";
        return os;
    }

    friend std::ostream & operator<<(std::ostream & os, const gsMaterialMatrixContainer & pwf)
    {
        return pwf.print(os);
    }

    /// Clear (delete) all functions
    void clear()
    {
        freeAll(m_container);
        m_container.clear();
    }

protected:
    Container m_container;
    index_t m_size;

};

} // namespace