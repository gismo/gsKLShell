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
    gsMaterialMatrixContainer( index_t size = 0 )
    {
        m_container.resize(size);
        // To do: initialize with null pointers
    }

    gsMaterialMatrixContainer(const gsMaterialMatrixContainer & other)
    {
        // for (index_t k=0; k!=other.m_container.size(); k++)
        //     add(memory::make_unique(other.m_container.at(k)));
        m_container = give(other.m_container);
    }

    ~gsMaterialMatrixContainer()
    {
        // freeAll(m_container);
    }

    /// Add a material matrix by copying argument
    void add(const gsMaterialMatrixBase<T> & mat)
    {
        m_container.push_back( memory::make_shared(mat.clone().release()) );
    }

    ///\brief Add a material matrix from a gsMaterialMatrixBase<T>::uPtr
    void add(const gsMaterialMatrixBase<T> * mat)
    {
        m_container.push_back( memory::make_shared_not_owned(mat) );
    }

    /// Set a material matrix by copying argument
    void set(const index_t i, const gsMaterialMatrixBase<T> & mat)
    {
        m_container[i] = memory::make_shared(mat.clone().release());
    }

    ///\brief Set a material matrix from a gsMaterialMatrixBase<T>::uPtr
    void set(const index_t i, const gsMaterialMatrixBase<T> * mat)
    {
        m_container[i] = memory::make_shared_not_owned(mat);
    }

    ///\brief Set a material matrix from a gsMaterialMatrixBase<T>::uPtr
    void set(const index_t i, const typename gsMaterialMatrixBase<T>::Ptr mat)
    {
        m_container[i] = mat;
    }

    gsMaterialMatrixBase<T> * piece(const index_t k) const
    {
        return m_container.at(k).get();
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

    /// Clear all function pointers
    void clear()
    {
        m_container.clear();
    }

protected:
    Container m_container;

};

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

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data)
    {
        GISMO_ERROR("Writing gsMaterialMatrixContainer to Xml is not implemented");
        // gsWarn<<"Writing gsMaterialMatrixContainer to Xml is not implemented\n";
        // gsXmlNode * result;
        // return result;
        // return putMaterialMatrixToXml< Object >( obj,data );
    }
};

} // namespace internal

} // namespace gismo
