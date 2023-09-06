/** 
 * @file DimensionTools.cpp
 * @brief Source of utility tools for the dimension IDs
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Enums/DimensionTools.hpp"

// Project includes
//

namespace QuICC {

namespace Dimensions {

   Transform::Id jump(Transform::Id id, int step)
   {
      if(id == Dimensions::Transform::TRA1D && step == -1)
      {
         return Dimensions::Transform::SPECTRAL;
      }
      else if(id == Dimensions::Transform::SPECTRAL && step == 1)
      {
         return Dimensions::Transform::TRA1D;
      }
      else
      {
         return static_cast<Transform::Id>(static_cast<int>(id)+step);
      }
   }
}
}
