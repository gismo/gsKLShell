/** @file gsMaterialMatrixContainer.hpp

    @brief Provides a container for material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
*/

#pragma once

// #include <gsKLShell/src/gsMaterialMatrixUtils.h>
// #include <gsKLShell/src/gsMaterialMatrixContainer.h>

using namespace gismo;

namespace gismo
{

template <class T>
gsMaterialMatrixContainer<T>::gsMaterialMatrixContainer( index_t size )
{
    m_container.resize(size);
    // To do: initialize with null pointers
}

    // gsMaterialMatrixContainer(const gsMaterialMatrixContainer & other)
    // {
    //     // for (index_t k=0; k!=other.m_container.size(); k++)
    //     //     add(memory::make_unique(other.m_container.at(k)));
    //     m_container = give(other.m_container);
    // }

template <class T>
gsMaterialMatrixContainer<T>::~gsMaterialMatrixContainer()
{
    // freeAll(m_container);
}
template <class T>
void gsMaterialMatrixContainer<T>::add(const gsMaterialMatrixBase<T> & mat)
{
    m_container.push_back( memory::make_shared(mat.clone().release()) );
}

template <class T>
void gsMaterialMatrixContainer<T>::add(const gsMaterialMatrixBase<T> * mat)
{
    m_container.push_back( memory::make_shared_not_owned(mat) );
}

template <class T>
void gsMaterialMatrixContainer<T>::set(const index_t i, const gsMaterialMatrixBase<T> & mat)
{
    m_container[i] = memory::make_shared(mat.clone().release());
}

template <class T>
void gsMaterialMatrixContainer<T>::set(const index_t i, const gsMaterialMatrixBase<T> * mat)
{
    m_container[i] = memory::make_shared_not_owned(mat);
}

template <class T>
void gsMaterialMatrixContainer<T>::set(const index_t i, const typename gsMaterialMatrixBase<T>::Ptr mat)
{
    m_container[i] = mat;
}

template <class T>
gsMaterialMatrixBase<T> * gsMaterialMatrixContainer<T>::piece(const index_t k) const
{
    return m_container.at(k).get();
}

template <class T>
std::ostream & gsMaterialMatrixContainer<T>::print(std::ostream &os) const
{
    os << "Material container with " << m_container.size() << " materials.\n";
    return os;
}

template <class T>
void gsMaterialMatrixContainer<T>::clear()
{
    m_container.clear();
}

} // end namespace
