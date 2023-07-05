/** @file gsMaterialMatrix.hpp

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

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrix.h>

namespace gismo
{

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

        if ( const gsMaterialMatrixLinear<2,T> * mm =
             dynamic_cast<const gsMaterialMatrixLinear<2,T> *>( ptr ) )
            return gsXml< gsMaterialMatrixLinear<2,T> >::put(*mm,data);
        if ( const gsMaterialMatrixLinear<3,T> * mm =
             dynamic_cast<const gsMaterialMatrixLinear<3,T> *>( ptr ) )
            return gsXml< gsMaterialMatrixLinear<3,T> >::put(*mm,data);

        // CompressibleNH2
        if ( const gsMaterialMatrix<2,T,11,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,11,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,11,true> >::put(*mm,data);
        // CompressibleNH3
        if ( const gsMaterialMatrix<3,T,11,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,11,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,11,true> >::put(*mm,data);
        // IncompressibleNH2
        if ( const gsMaterialMatrix<2,T,11,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,11,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,11,false> >::put(*mm,data);
        // IncompressibleNH3
        if ( const gsMaterialMatrix<3,T,11,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,11,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,11,false> >::put(*mm,data);

        // CompressibleNHe2
        if ( const gsMaterialMatrix<2,T,12,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,12,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,12,true> >::put(*mm,data);
        // CompressibleNHe3
        if ( const gsMaterialMatrix<3,T,12,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,12,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,12,true> >::put(*mm,data);
        // IncompressibleNHe2
        if ( const gsMaterialMatrix<2,T,12,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,12,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,12,false> >::put(*mm,data);
        // IncompressibleNHe3
        if ( const gsMaterialMatrix<3,T,12,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,12,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,12,false> >::put(*mm,data);

        // CompressibleMR2
        if ( const gsMaterialMatrix<2,T,13,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,13,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,13,true> >::put(*mm,data);
        // CompressibleMR3
        if ( const gsMaterialMatrix<3,T,13,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,13,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,13,true> >::put(*mm,data);
        // IncompressibleMR2
        if ( const gsMaterialMatrix<2,T,13,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,13,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,13,false> >::put(*mm,data);
        // IncompressibleMR3
        if ( const gsMaterialMatrix<3,T,13,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,13,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,13,false> >::put(*mm,data);

        // CompressibleOG2
        if ( const gsMaterialMatrix<2,T,34,true> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,34,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,34,true> >::put(*mm,data);
        // CompressibleOG3
        if ( const gsMaterialMatrix<3,T,34,true> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,34,true> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,34,true> >::put(*mm,data);
        // IncompressibleOG2
        if ( const gsMaterialMatrix<2,T,34,false> * mm =
             dynamic_cast<const gsMaterialMatrix<2,T,34,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<2,T,34,false> >::put(*mm,data);
        // IncompressibleOG3
        if ( const gsMaterialMatrix<3,T,34,false> * mm =
             dynamic_cast<const gsMaterialMatrix<3,T,34,false> *>( ptr ) )
            return gsXml< gsMaterialMatrix<3,T,34,false> >::put(*mm,data);

        gsWarn<<"gsMaterialMatrixBase: put<MaterialMatrixBase<T>>: No known MaterialMatrix "<< obj <<"Error.\n";
        return NULL;
    }
};

}// end namespace internal

} // end namespace
