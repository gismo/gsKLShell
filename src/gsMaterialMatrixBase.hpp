/** @file gsMaterialMatrixNonlinear.hpp

    @brief Provides hyperelastic material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

/*
    To Do [updated 16-06-2020]:
    - Make beta (compressible materials) and material parameters universal for all integration points over the thickness. So get them out of the _dPsi functions etc and move them into the integration loops as global variables.

*/



#pragma once

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.h>

namespace gismo
{

// template <class T>
// std::ostream & gsMaterialMatrixBase<T>::print(std::ostream &os) const
// {
//     // TFT
//     if (const gsMaterialMatrixTFT<2,T,true> * mm =
//          dynamic_cast<const gsMaterialMatrixTFT<2,T,true> *>( this ) )
//         return mm->print(os);
//     if (const gsMaterialMatrixTFT<2,T,false> * mm =
//          dynamic_cast<const gsMaterialMatrixTFT<2,T,false> *>( this ) )
//         return mm->print(os);
//     if (const gsMaterialMatrixTFT<3,T,true> * mm =
//          dynamic_cast<const gsMaterialMatrixTFT<3,T,true> *>( this ) )
//         return mm->print(os);
//     if (const gsMaterialMatrixTFT<3,T,false> * mm =
//          dynamic_cast<const gsMaterialMatrixTFT<3,T,false> *>( this ) )
//         return mm->print(os);

//     // Linear
//     if ( const gsMaterialMatrixLinear<2,T> * mm =
//          dynamic_cast<const gsMaterialMatrixLinear<2,T> *>( this ) )
//         return mm->print(os);
//     if ( const gsMaterialMatrixLinear<3,T> * mm =
//          dynamic_cast<const gsMaterialMatrixLinear<3,T> *>( this ) )
//         return mm->print(os);

//     // CompressibleNH2
//     if ( const gsMaterialMatrix<2,T,11,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,11,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleNH3
//     if ( const gsMaterialMatrix<3,T,11,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,11,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNH2
//     if ( const gsMaterialMatrix<2,T,11,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,11,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNH3
//     if ( const gsMaterialMatrix<3,T,11,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,11,false> *>( this ) )
//         return mm->print(os);

//     // CompressibleNHe2
//     if ( const gsMaterialMatrix<2,T,12,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,12,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleNHe3
//     if ( const gsMaterialMatrix<3,T,12,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,12,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNHe2
//     if ( const gsMaterialMatrix<2,T,12,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,12,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNHe3
//     if ( const gsMaterialMatrix<3,T,12,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,12,false> *>( this ) )
//         return mm->print(os);

//     // CompressibleMR2
//     if ( const gsMaterialMatrix<2,T,13,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,13,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleMR3
//     if ( const gsMaterialMatrix<3,T,13,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,13,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleMR2
//     if ( const gsMaterialMatrix<2,T,13,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,13,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleMR3
//     if ( const gsMaterialMatrix<3,T,13,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,13,false> *>( this ) )
//         return mm->print(os);

//     // CompressibleOG2
//     if ( const gsMaterialMatrix<2,T,34,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,34,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleOG3
//     if ( const gsMaterialMatrix<3,T,34,true> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,34,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleOG2
//     if ( const gsMaterialMatrix<2,T,34,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<2,T,34,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleOG3
//     if ( const gsMaterialMatrix<3,T,34,false> * mm =
//          dynamic_cast<const gsMaterialMatrix<3,T,34,false> *>( this ) )
//         return mm->print(os);

//     os<<"gsMaterialMatrixBase (type not understood).\n";
//     return os;
// }

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
                        ( s == "IncompressibleNHe2" ) ||
                        ( s == "CompressibleMR2"    ) ||
                        ( s == "IncompressibleMR2"  ) ||
                        ( s == "CompressibleOG2"    ) ||
                        ( s == "IncompressibleOG2"  ) )
                return new gsMaterialMatrixTFT<2,T,false>(tmp);
            else if (   ( s == "CompressibleNH3"    ) ||
                        ( s == "IncompressibleNH3"  ) ||
                        ( s == "CompressibleNHe3"   ) ||
                        ( s == "IncompressibleNHe3" ) ||
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
            return gsXml< gsMaterialMatrix<2,T,11,true> >::get(node);
        if ( s == "CompressibleNH3"    )
            return gsXml< gsMaterialMatrix<3,T,11,true> >::get(node);
        if ( s == "IncompressibleNH2"  )
            return gsXml< gsMaterialMatrix<2,T,11,false> >::get(node);
        if ( s == "IncompressibleNH3"  )
            return gsXml< gsMaterialMatrix<3,T,11,false> >::get(node);

        if ( s == "CompressibleNHe2"    )
            return gsXml< gsMaterialMatrix<2,T,12,true> >::get(node);
        if ( s == "CompressibleNHe3"    )
            return gsXml< gsMaterialMatrix<3,T,12,true> >::get(node);
        if ( s == "IncompressibleNHe2"  )
            return gsXml< gsMaterialMatrix<2,T,12,false> >::get(node);
        if ( s == "IncompressibleNHe3"  )
            return gsXml< gsMaterialMatrix<3,T,12,false> >::get(node);

        if ( s == "CompressibleMR2"    )
            return gsXml< gsMaterialMatrix<2,T,13,true> >::get(node);
        if ( s == "CompressibleMR3"    )
            return gsXml< gsMaterialMatrix<3,T,13,true> >::get(node);
        if ( s == "IncompressibleMR2"  )
            return gsXml< gsMaterialMatrix<2,T,13,false> >::get(node);
        if ( s == "IncompressibleMR3"  )
            return gsXml< gsMaterialMatrix<3,T,13,false> >::get(node);

        if ( s == "CompressibleOG2"    )
            return gsXml< gsMaterialMatrix<2,T,34,true> >::get(node);
        if ( s == "CompressibleOG3"    )
            return gsXml< gsMaterialMatrix<3,T,34,true> >::get(node);
        if ( s == "IncompressibleOG2"  )
            return gsXml< gsMaterialMatrix<2,T,34,false> >::get(node);
        if ( s == "IncompressibleOG3"  )
            return gsXml< gsMaterialMatrix<3,T,34,false> >::get(node);

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
        if ( const gsMaterialMatrix<2,T,11,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,11,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,11,true> >::put(*mm,data);
        // CompressibleNH3
        if ( const gsMaterialMatrix<3,T,11,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,11,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,11,true> >::put(*mm,data);
        // IncompressibleNH2
        if ( const gsMaterialMatrix<2,T,11,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,11,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,11,false> >::put(*mm,data);
        // IncompressibleNH3
        if ( const gsMaterialMatrix<3,T,11,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,11,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,11,false> >::put(*mm,data);

        // CompressibleNHe2
        if ( const gsMaterialMatrix<2,T,12,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,12,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,12,true> >::put(*mm,data);
        // CompressibleNHe3
        if ( const gsMaterialMatrix<3,T,12,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,12,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,12,true> >::put(*mm,data);
        // IncompressibleNHe2
        if ( const gsMaterialMatrix<2,T,12,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,12,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,12,false> >::put(*mm,data);
        // IncompressibleNHe3
        if ( const gsMaterialMatrix<3,T,12,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,12,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,12,false> >::put(*mm,data);

        // CompressibleMR2
        if ( const gsMaterialMatrix<2,T,13,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,13,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,13,true> >::put(*mm,data);
        // CompressibleMR3
        if ( const gsMaterialMatrix<3,T,13,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,13,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,13,true> >::put(*mm,data);
        // IncompressibleMR2
        if ( const gsMaterialMatrix<2,T,13,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,13,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,13,false> >::put(*mm,data);
        // IncompressibleMR3
        if ( const gsMaterialMatrix<3,T,13,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,13,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,13,false> >::put(*mm,data);

        // CompressibleOG2
        if ( const gsMaterialMatrix<2,T,34,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,34,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,34,true> >::put(*mm,data);
        // CompressibleOG3
        if ( const gsMaterialMatrix<3,T,34,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,34,true> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,34,true> >::put(*mm,data);
        // IncompressibleOG2
        if ( const gsMaterialMatrix<2,T,34,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,34,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<2,T,34,false> >::put(*mm,data);
        // IncompressibleOG3
        if ( const gsMaterialMatrix<3,T,34,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,34,false> *>( obj ) )
            return gsXml< gsMaterialMatrix<3,T,34,false> >::put(*mm,data);

        gsWarn<<"gsMaterialMatrixBase: put<MaterialMatrixBase<T>>: No known MaterialMatrix "<< obj <<"Error.\n";
        return NULL;
    }
};

}// end namespace internal

} // end namespace
