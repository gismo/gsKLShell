/** @file gsMaterialMatrix.h

    @brief Provides material matrices for the thin shell class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixComposite.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

// template<short_t d, class T>
// gsMaterialMatrixBase<T> * getMaterialMatrix(
//                                 const gsMultiPatch<T>               & mp,
//                                 const gsMultiPatch<T>               & mp_def,
//                                 const gsFunction<T>                 & thickness,
//                                 const std::vector<gsFunction<T> *>  & parameters,
//                                 const gsFunction<T>                 & rho,
//                                 const gsOptionList                  & options
//                                 );

// template<short_t d, class T>
// gsMaterialMatrixBase<T> * getMaterialMatrix(
//                                 const gsMultiPatch<T>               & mp,
//                                 const gsFunction<T>                 & thickness,
//                                 const std::vector<gsFunction<T> *>  & parameters,
//                                 const gsOptionList                  & options
//                                 )
// {
//     return getMaterialMatrix<d,T>(mp,NULL,thickness,parameters,NULL,options);
// }

// template<short_t d, class T>
// gsMaterialMatrixBase<T> * getMaterialMatrix(
//                                 const gsMultiPatch<T>               & mp,
//                                 const gsFunction<T>                 & thickness,
//                                 const std::vector<gsFunction<T> *>  & parameters,
//                                 const gsFunction<T>                 & rho,
//                                 const gsOptionList                  & options
//                                 )
// {
//     return getMaterialMatrix<d,T>(mp,NULL,thickness,parameters,rho,options);
// }

// template<short_t d, class T>
// gsMaterialMatrixBase<T> * getMaterialMatrix(
//                                 const gsMultiPatch<T>               & mp,
//                                 const gsMultiPatch<T>               & mp_def,
//                                 const gsFunction<T>                 & thickness,
//                                 const std::vector<gsFunction<T> *>  & parameters,
//                                 const gsOptionList                  & options
//                                 )
// {
//     return getMaterialMatrix<d,T>(mp,mp_def,thickness,parameters,NULL,options);
// }

template<short_t d, class T>
gsMaterialMatrixBase<T> * getMaterialMatrix(
                                const gsMultiPatch<T>               & mp,
                                const gsMultiPatch<T>               & mp_def,
                                const gsFunction<T>                 & thickness,
                                const std::vector<gsFunction<T> *>  & parameters,
                                const gsFunction<T>                 & rho,
                                const gsOptionList                  & options
                                )
{
    index_t material        = options.askInt("Material",0); // SvK by default
    bool compressibility    = options.askSwitch("Compressibility",false);
    index_t implementation  = options.askInt("Implementation",1);

    enum Material       mat     = static_cast<enum Material>(material);
    enum Implementation impl    = static_cast<enum Implementation>(implementation);

    if      (mat==Material::SvK)
    {
        if (impl==Implementation::Composite)
                return new gsMaterialMatrixComposite<d,T>(mp,mp_def,thickness,parameters,rho);
        else
                return new gsMaterialMatrixLinear<d,T>(mp,mp_def,thickness,parameters,rho);
    }
    else
    {
        if (impl==Implementation::Composite)
            GISMO_ERROR("Hyperelastic composites not available");
        //--------------------------------------NH Incompressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,mp_def,thickness,parameters,rho);
        }
        //--------------------------------------NH Compressible--------------------------------------------------
        else if ((mat==Material::NH) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::NH) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        //
        //--------------------------------------NH_ext Incompressible--------------------------------------------------
        else if ((mat==Material::NH_ext) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Extended NH with not available");
        }
        //---------------------------------------NH_ext Compressible-------------------------------------------------
        else if ((mat==Material::NH_ext) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::NH_ext) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::NH_ext, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        //
        //--------------------------------------MR Incompressible--------------------------------------------------
        else if ((mat==Material::MR) && (impl==Implementation::Analytical) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,mp_def,thickness,parameters,rho);
        }
        //---------------------------------------MR Compressible-------------------------------------------------
        else if      ((mat==Material::MR) && (impl==Implementation::Analytical) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Analytical>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Generalized) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Generalized>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        else if ((mat==Material::MR) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::MR, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        //
        //--------------------------------------OG Incompressible--------------------------------------------------
        else if ((mat==Material::OG) && (impl==Implementation::Analytical) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Ogden with analytical implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Generalized) && (!compressibility))
        {
            GISMO_ERROR("Incompressible Ogden with generalized implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Spectral) && (!compressibility))
        {
            constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,false>(mp,mp_def,thickness,parameters,rho);
        }
        //---------------------------------------OG Compressible-------------------------------------------------
        else if      ((mat==Material::OG) && (impl==Implementation::Analytical) && (compressibility))
        {
            GISMO_ERROR("Compressible Ogden with analytical implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Generalized) && (compressibility))
        {
            GISMO_ERROR("Compressible Ogden with generalized implementation not available");
        }
        else if ((mat==Material::OG) && (impl==Implementation::Spectral) && (compressibility))
        {
            constexpr int id = encodeMat_id<Material::OG, Implementation::Spectral>::id;
            return new gsMaterialMatrix<d,real_t,id,true>(mp,mp_def,thickness,parameters,rho);
        }
        else
            return NULL;
    }

}



} // namespace
