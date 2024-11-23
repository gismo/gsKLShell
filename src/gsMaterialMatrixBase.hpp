/** @file gsMaterialMatrixBase.hpp

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

#include <gsIO/gsXml.h>
// #include <gsKLShell/src/gsMaterialMatrixBase.h>
// #include <gsKLShell/src/gsMaterialMatrixLinear.h>
// #include <gsKLShell/src/gsMaterialMatrixLinear.h>
// #include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
// #include <gsKLShell/src/gsMaterialMatrixTFT.h>

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
//     if ( const gsMaterialMatrixNonlinear<2,T,11,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,11,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleNH3
//     if ( const gsMaterialMatrixNonlinear<3,T,11,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,11,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNH2
//     if ( const gsMaterialMatrixNonlinear<2,T,11,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,11,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNH3
//     if ( const gsMaterialMatrixNonlinear<3,T,11,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,11,false> *>( this ) )
//         return mm->print(os);

//     // CompressibleNHe2
//     if ( const gsMaterialMatrixNonlinear<2,T,12,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,12,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleNHe3
//     if ( const gsMaterialMatrixNonlinear<3,T,12,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,12,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNHe2
//     if ( const gsMaterialMatrixNonlinear<2,T,12,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,12,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleNHe3
//     if ( const gsMaterialMatrixNonlinear<3,T,12,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,12,false> *>( this ) )
//         return mm->print(os);

//     // CompressibleMR2
//     if ( const gsMaterialMatrixNonlinear<2,T,13,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,13,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleMR3
//     if ( const gsMaterialMatrixNonlinear<3,T,13,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,13,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleMR2
//     if ( const gsMaterialMatrixNonlinear<2,T,13,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,13,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleMR3
//     if ( const gsMaterialMatrixNonlinear<3,T,13,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,13,false> *>( this ) )
//         return mm->print(os);

//     // CompressibleOG2
//     if ( const gsMaterialMatrixNonlinear<2,T,34,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,34,true> *>( this ) )
//         return mm->print(os);
//     // CompressibleOG3
//     if ( const gsMaterialMatrixNonlinear<3,T,34,true> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,34,true> *>( this ) )
//         return mm->print(os);
//     // IncompressibleOG2
//     if ( const gsMaterialMatrixNonlinear<2,T,34,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<2,T,34,false> *>( this ) )
//         return mm->print(os);
//     // IncompressibleOG3
//     if ( const gsMaterialMatrixNonlinear<3,T,34,false> * mm =
//          dynamic_cast<const gsMaterialMatrixNonlinear<3,T,34,false> *>( this ) )
//         return mm->print(os);

//     os<<"gsMaterialMatrixBase (type not understood).\n";
//     return os;
// }


} // end namespace
