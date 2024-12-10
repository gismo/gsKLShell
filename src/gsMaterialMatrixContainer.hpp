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
#include <gsKLShell/src/gsMaterialMatrixXml.hpp>

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


namespace internal
{
/// @brief get a MaterialMatrixContainer from XML data
///
/// \ingroup KLShell
template<class T>
class gsXml< gsMaterialMatrixContainer<T> >
{
private:
    gsXml() { }
    typedef gsMaterialMatrixContainer<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(gsMaterialMatrixContainer<T>);
    static std::string tag ()  { return "MaterialMatrixContainer"; }
    static std::string type () { return ""; }

    GSXML_GET_POINTER(Object);

    static void get_into(gsXmlNode * node,Object & obj)
    {
        const int size = atoi(node->first_attribute("size")->value());

        // Read material inventory
        int count = countByTag("MaterialMatrix", node);
        std::vector<typename gsMaterialMatrixBase<T>::Ptr> mat(count);
        for (gsXmlNode * child = node->first_node("MaterialMatrix"); child; child =
                child->next_sibling("MaterialMatrix"))
        {
            const int i = atoi(child->first_attribute("index")->value());
            mat[i] = memory::make_shared(gsXml<gsMaterialMatrixBase<T>>::get(child));
        }

        obj = gsMaterialMatrixContainer<T>(size);
        for (gsXmlNode * child = node->first_node("group"); child;
                child = child->next_sibling("group"))
        {
            const int mIndex = atoi(child->first_attribute("material")->value());
            std::istringstream group_str;
            group_str.str( child->value() );

            for(int patch; ( gsGetInt(group_str,patch)); )
                obj.set(patch,mat[mIndex]);
        }

    }

    static gsXmlNode * put (const Object & /* obj */,
                            gsXmlTree & /* data */)
    {
        GISMO_ERROR("Writing gsMaterialMatrixContainer to Xml is not implemented");
        // gsWarn<<"Writing gsMaterialMatrixContainer to Xml is not implemented\n";
        // gsXmlNode * result;
        // return result;
        // return putMaterialMatrixToXml< Object >( obj,data );
    }
};
}

} // end namespace
