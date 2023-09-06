/**
 * @file SHm2lIndexConv.cpp
 * @brief Source of the index converter for SH from M to L ordering
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/SHm2lIndexConv.hpp"

// Project includes
//
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Enums/Dimensions.hpp"

namespace QuICC {

namespace Parallel {

   void SHm2lIndexConv::init(const Resolution& res, const Dimensions::Transform::Id id)
   {
      this->mMinL = res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
   }

}
}
