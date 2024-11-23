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

#include <gsKLShell/src/gsMaterialMatrixBase.h>

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
class gsMaterialMatrixContainer // change name to PtrContainer
{
public:
    typedef typename std::vector<typename gsMaterialMatrixBase<T>::Ptr> Container;

    typedef typename Container::iterator       iterator;
    typedef typename Container::const_iterator const_iterator;

    /// Shared pointer for gsMaterialMatrixContainer
    typedef memory::shared_ptr< gsMaterialMatrixContainer > Ptr;

    /// Unique pointer for gsMaterialMatrixContainer
    typedef memory::unique_ptr< gsMaterialMatrixContainer > uPtr;
public:

    /// Constructor
    gsMaterialMatrixContainer( index_t size = 0 );

    ~gsMaterialMatrixContainer();

    /// Add a material matrix by copying argument
    void add(const gsMaterialMatrixBase<T> & mat);

    ///\brief Add a material matrix from a gsMaterialMatrixBase<T> *
    void add(const gsMaterialMatrixBase<T> * mat);

    /// Set a material matrix by copying argument
    void set(const index_t i, const gsMaterialMatrixBase<T> & mat);

    ///\brief Set a material matrix from a gsMaterialMatrixBase<T> *
    void set(const index_t i, const gsMaterialMatrixBase<T> * mat);

    ///\brief Set a material matrix from a gsMaterialMatrixBase<T>::uPtr
    void set(const index_t i, const typename gsMaterialMatrixBase<T>::Ptr mat);

    gsMaterialMatrixBase<T> * piece(const index_t k) const;

    index_t size() const {return m_container.size();}

    std::ostream &print(std::ostream &os) const;

    friend std::ostream & operator<<(std::ostream & os, const gsMaterialMatrixContainer & mc)
    {
        return mc.print(os);
    }

    /// Clear all function pointers
    void clear();

protected:
    Container m_container;

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixContainer.hpp)
#endif
