/**
 * @file PMIndexConv.cpp
 * @brief Source of the index converter that splits array for positive and negative FFT frequencies
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/PMIndexConv.hpp"

// Project includes
//
#include "QuICC/Enums/DimensionTools.hpp"

namespace QuICC {

namespace Parallel {

   void PMIndexConv::init(const Resolution& res, const Dimensions::Transform::Id id)
   {
      /// \mhdBug Calculation of mSimN with double static_cast is not ideal

      this->mspTResFwd = res.cpu()->dim(id);
      this->mspTResBwd = res.cpu()->dim(Dimensions::jump(id,1));
      this->mSimN = res.sim().dim(static_cast<Dimensions::Simulation::Id>(static_cast<int>(id)+1), Dimensions::Space::SPECTRAL);
   }

}
}
