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

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.h>


namespace gismo {

namespace internal {

/// Get a Geometry from XML data
template<class T>
class gsXml< gsMaterialMatrixBase<T> >
{
private:
    gsXml() { }
    typedef gsMaterialMatrixBase<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "MaterialMatrix"; }
    static std::string type () { return ""; }

    // GSXML_GET_INTO(Object)

    static void get_into(gsXmlNode * node, Object & obj)
    {
        obj = *get(node);
        // obj = *tmp;
    }

    static Object * get(gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"MaterialMatrix") ),
                      "Something went wrong, was waiting for a MaterialMatrix tag.\n" );

        gsXmlAttribute * gtype = node->first_attribute("type");
        if ( ! gtype )
        {
            gsWarn<< "MaterialMatrix without a type in the xml file\n";
            return NULL;
        }

        const gsXmlAttribute * att_TFT = node->first_attribute("TFT");
        bool TFT = false;
        if (NULL != att_TFT)
        {
            if ((strcmp("1", att_TFT->value()) == 0) ||
                (strcmp("true", att_TFT->value()) == 0) ||
                (strcmp("True", att_TFT->value()) == 0)   )
            TFT = true;
        }

        std::string s = gtype->value() ;
        Object * tmp = get_impl(node);
        if (TFT)
        {
            if      (   ( s == "Linear2"            ) )
                return new gsMaterialMatrixTFT<2,T,true>(tmp);
            else if (   ( s == "Linear3"            ) )
                return new gsMaterialMatrixTFT<3,T,true>(tmp);
            else if (   ( s == "CompressibleNH2"    ) ||
                        ( s == "IncompressibleNH2"  ) ||
                        ( s == "CompressibleNHe2"   ) ||
                        // ( s == "IncompressibleNHe2" ) ||
                        ( s == "CompressibleMR2"    ) ||
                        ( s == "IncompressibleMR2"  ) ||
                        ( s == "CompressibleOG2"    ) ||
                        ( s == "IncompressibleOG2"  ) )
                return new gsMaterialMatrixTFT<2,T,false>(tmp);
            else if (   ( s == "CompressibleNH3"    ) ||
                        ( s == "IncompressibleNH3"  ) ||
                        ( s == "CompressibleNHe3"   ) ||
                        // ( s == "IncompressibleNHe3" ) ||
                        ( s == "CompressibleMR3"    ) ||
                        ( s == "IncompressibleMR3"  ) ||
                        ( s == "CompressibleOG3"    ) ||
                        ( s == "IncompressibleOG3"  ) )
                return new gsMaterialMatrixTFT<3,T,false>(tmp);
            else
            {
                gsWarn<<"Material matrix for TFT model not recognised\n";
                return NULL;
            }
        }
        else
            return tmp;
    }

    static Object * get_impl(gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"MaterialMatrix") ),
                      "Something went wrong, was waiting for a MaterialMatrix tag.\n" );

        gsXmlAttribute * gtype = node->first_attribute("type");
        if ( ! gtype )
        {
            gsWarn<< "MaterialMatrix without a type in the xml file\n";
            return NULL;
        }

        std::string s = gtype->value() ;
        if ( s == "Linear2"    )
            return gsXml< gsMaterialMatrixLinear<2,T> >::get(node);
        if ( s == "Linear3"    )
            return gsXml< gsMaterialMatrixLinear<3,T> >::get(node);

        if ( s == "CompressibleNH2"    )
            return gsXml< gsMaterialMatrixNonlinear<2,T,11,true> >::get(node);
        if ( s == "CompressibleNH3"    )
            return gsXml< gsMaterialMatrixNonlinear<3,T,11,true> >::get(node);
        if ( s == "IncompressibleNH2"  )
            return gsXml< gsMaterialMatrixNonlinear<2,T,11,false> >::get(node);
        if ( s == "IncompressibleNH3"  )
            return gsXml< gsMaterialMatrixNonlinear<3,T,11,false> >::get(node);

        if ( s == "CompressibleNHe2"    )
            return gsXml< gsMaterialMatrixNonlinear<2,T,12,true> >::get(node);
        if ( s == "CompressibleNHe3"    )
            return gsXml< gsMaterialMatrixNonlinear<3,T,12,true> >::get(node);
        // if ( s == "IncompressibleNHe2"  )
        //     return gsXml< gsMaterialMatrixNonlinear<2,T,12,false> >::get(node);
        // if ( s == "IncompressibleNHe3"  )
        //     return gsXml< gsMaterialMatrixNonlinear<3,T,12,false> >::get(node);

        if ( s == "CompressibleMR2"    )
            return gsXml< gsMaterialMatrixNonlinear<2,T,13,true> >::get(node);
        if ( s == "CompressibleMR3"    )
            return gsXml< gsMaterialMatrixNonlinear<3,T,13,true> >::get(node);
        if ( s == "IncompressibleMR2"  )
            return gsXml< gsMaterialMatrixNonlinear<2,T,13,false> >::get(node);
        if ( s == "IncompressibleMR3"  )
            return gsXml< gsMaterialMatrixNonlinear<3,T,13,false> >::get(node);

        if ( s == "CompressibleOG2"    )
            return gsXml< gsMaterialMatrixNonlinear<2,T,34,true> >::get(node);
        if ( s == "CompressibleOG3"    )
            return gsXml< gsMaterialMatrixNonlinear<3,T,34,true> >::get(node);
        if ( s == "IncompressibleOG2"  )
            return gsXml< gsMaterialMatrixNonlinear<2,T,34,false> >::get(node);
        if ( s == "IncompressibleOG3"  )
            return gsXml< gsMaterialMatrixNonlinear<3,T,34,false> >::get(node);

        gsWarn<<"gsMaterialMatrixBase: get<MaterialMatrixBase<T>>: No known MaterialMatrix \""<<s<<"\". Error.\n";
        return NULL;
    }

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data)
    {
        const Object * ptr = & obj;
        if (dynamic_cast<const gsMaterialMatrixTFT<2,T,true> *>( ptr )
            ||
            dynamic_cast<const gsMaterialMatrixTFT<3,T,true> *>( ptr )
            ||
            dynamic_cast<const gsMaterialMatrixTFT<2,T,false> *>( ptr )
            ||
            dynamic_cast<const gsMaterialMatrixTFT<3,T,false> *>( ptr )
           )
        {
            gsXmlNode * tmp = put_impl(ptr->material(),data);
            tmp->append_attribute(internal::makeAttribute("TFT","true",data));
            return tmp;
        }
        else
            return put_impl(ptr,data);
    }

    static gsXmlNode * put_impl (const Object * obj,
                            gsXmlTree & data)
    {
        if ( const gsMaterialMatrixLinear<2,T> * mm =
             dynamic_cast<const gsMaterialMatrixLinear<2,T> *>( obj ) )
            return gsXml< gsMaterialMatrixLinear<2,T> >::put(*mm,data);
        if ( const gsMaterialMatrixLinear<3,T> * mm =
             dynamic_cast<const gsMaterialMatrixLinear<3,T> *>( obj ) )
            return gsXml< gsMaterialMatrixLinear<3,T> >::put(*mm,data);

        // CompressibleNH2
        if ( const gsMaterialMatrixNonlinear<2,T,11,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<2,T,11,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<2,T,11,true> >::put(*mm,data);
        // CompressibleNH3
        if ( const gsMaterialMatrixNonlinear<3,T,11,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<3,T,11,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<3,T,11,true> >::put(*mm,data);
        // IncompressibleNH2
        if ( const gsMaterialMatrixNonlinear<2,T,11,false> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<2,T,11,false> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<2,T,11,false> >::put(*mm,data);
        // IncompressibleNH3
        if ( const gsMaterialMatrixNonlinear<3,T,11,false> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<3,T,11,false> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<3,T,11,false> >::put(*mm,data);

        // CompressibleNHe2
        if ( const gsMaterialMatrixNonlinear<2,T,12,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<2,T,12,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<2,T,12,true> >::put(*mm,data);
        // CompressibleNHe3
        if ( const gsMaterialMatrixNonlinear<3,T,12,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<3,T,12,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<3,T,12,true> >::put(*mm,data);
        // // IncompressibleNHe2
        // if ( const gsMaterialMatrixNonlinear<2,T,12,false> * mm =
        //      dynamic_cast<const gsMaterialMatrixNonlinear<2,T,12,false> *>( obj ) )
        //     return gsXml< gsMaterialMatrixNonlinear<2,T,12,false> >::put(*mm,data);
        // // IncompressibleNHe3
        // if ( const gsMaterialMatrixNonlinear<3,T,12,false> * mm =
        //      dynamic_cast<const gsMaterialMatrixNonlinear<3,T,12,false> *>( obj ) )
        //     return gsXml< gsMaterialMatrixNonlinear<3,T,12,false> >::put(*mm,data);

        // CompressibleMR2
        if ( const gsMaterialMatrixNonlinear<2,T,13,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<2,T,13,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<2,T,13,true> >::put(*mm,data);
        // CompressibleMR3
        if ( const gsMaterialMatrixNonlinear<3,T,13,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<3,T,13,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<3,T,13,true> >::put(*mm,data);
        // IncompressibleMR2
        if ( const gsMaterialMatrixNonlinear<2,T,13,false> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<2,T,13,false> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<2,T,13,false> >::put(*mm,data);
        // IncompressibleMR3
        if ( const gsMaterialMatrixNonlinear<3,T,13,false> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<3,T,13,false> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<3,T,13,false> >::put(*mm,data);

        // CompressibleOG2
        if ( const gsMaterialMatrixNonlinear<2,T,34,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<2,T,34,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<2,T,34,true> >::put(*mm,data);
        // CompressibleOG3
        if ( const gsMaterialMatrixNonlinear<3,T,34,true> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<3,T,34,true> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<3,T,34,true> >::put(*mm,data);
        // IncompressibleOG2
        if ( const gsMaterialMatrixNonlinear<2,T,34,false> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<2,T,34,false> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<2,T,34,false> >::put(*mm,data);
        // IncompressibleOG3
        if ( const gsMaterialMatrixNonlinear<3,T,34,false> * mm =
             dynamic_cast<const gsMaterialMatrixNonlinear<3,T,34,false> *>( obj ) )
            return gsXml< gsMaterialMatrixNonlinear<3,T,34,false> >::put(*mm,data);

        gsWarn<<"gsMaterialMatrixBase: put<MaterialMatrixBase<T>>: No known MaterialMatrix "<< obj <<"Error.\n";
        return NULL;
    }
};

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
    gsXml<gsFunctionExpr<T> >::get_into(tmp->first_node("Function"), thickness);
    result.setThickness(thickness);

    tmp = node->first_node("Density");
    gsFunctionExpr<T> density;
    bool hasDensity = tmp;
    if ( hasDensity )
    {
        gsXml<gsFunctionExpr<T> >::get_into(tmp->first_node("Function"), density);
        result.setDensity(density);
    }

    gsXmlNode * parNode = node->first_node("Parameters");
    // Read function inventory
    gsFunctionExpr<T> fun;
    for (gsXmlNode * child = parNode->first_node("Function"); child; child =
            child->next_sibling("Function"))
    {
        const int i = atoi(child->first_attribute("index")->value());
        gsXml<gsFunctionExpr<T> >::get_into(child, fun);
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
