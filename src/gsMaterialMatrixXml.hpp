/** @file gsMaterialMatrixXml.hpp

    @brief TO DO

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): TO DO
*/

#pragma once

#include <fstream>
#include <gsCore/gsFunctionExpr.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo {

namespace internal {

template<class Object>
Object getMaterialMatrixFromXml ( gsXmlNode * node)
{
    typedef typename Object::Scalar_t T;

    //gsWarn<<"Reading "<< gsXml<Object>::type() <<" Geometry..\n";
    assert ( ( !strcmp( node->name(),"MaterialMatrix") ) &&
             ( !strcmp(node->first_attribute("type")->value(), gsXml<Object>::type().c_str() ) ) );

    Object result = Object();

    gsXmlNode * tmp;

    tmp = node->first_node("Thickness");
    GISMO_ASSERT(tmp,"Thickness must be assigned!");
    gsFunctionExpr<T> thickness;
    getFunctionFromXml(tmp->first_node("Function"), thickness);
    result.setThickness(thickness);

    tmp = node->first_node("Density");
    gsFunctionExpr<T> density;
    bool hasDensity = tmp;
    if ( hasDensity )
    {
        getFunctionFromXml(tmp->first_node("Function"), density);
        result.setDensity(density);
    }

    gsXmlNode * parNode = node->first_node("Parameters");
    // Read function inventory
    gsFunctionExpr<T> fun;
    for (gsXmlNode * child = parNode->first_node("Function"); child; child =
            child->next_sibling("Function"))
    {
        const int i = atoi(child->first_attribute("index")->value());
        getFunctionFromXml(child, fun);
        result.setParameter(i,fun);
    }
    return result;
}

template<class Object>
gsXmlNode * putMaterialMatrixToXml ( Object const & obj, gsXmlTree & data)
{

    typedef typename Object::Scalar_t T;

    // Make a new XML Geometry node
    gsXmlNode * mm = internal::makeNode("MaterialMatrix", data);
    mm->append_attribute( makeAttribute("type",
                                        internal::gsXml<Object>::type().c_str(), data) );

    GISMO_ASSERT(obj.hasThickness(),"Thickness is not assigned");
    gsXmlNode * t = internal::makeNode("Thickness", data);
    gsXmlNode * tfun = putFunctionToXml<T>(obj.getThickness(), data, 0);
    t->append_node(tfun);
    mm->append_node(t);
    if (obj.hasDensity())
    {
        gsXmlNode * r = internal::makeNode("Density", data);
        gsXmlNode * rfun = putFunctionToXml<T>(obj.getDensity(), data, 0);
        r->append_node(rfun);
        mm->append_node(r);
    }

    gsXmlNode * p = internal::makeNode("Parameters", data);
    for (index_t k=0; k!=obj.numParameters(); k++)
    {
        gsXmlNode * pfun = putFunctionToXml<T>(obj.getParameter(k), data, k);
        p->append_node(pfun);
    }
    mm->append_node(p);

    return mm;
}

}// end namespace internal

}// end namespace gismo

//#undef GSXML_COMMON_FUNCTIONS
//#undef TMPLA2
//#undef TMPLA3
