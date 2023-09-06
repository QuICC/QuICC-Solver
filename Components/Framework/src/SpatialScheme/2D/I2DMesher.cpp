/** 
 * @file I2DMesher.cpp
 * @brief Source of the 1D spatial scheme mesher
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SpatialScheme/2D/I2DMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   I2DMesher::I2DMesher(const GridPurpose::Id purpose)
      : IMesher(purpose)
   {
   }

   int I2DMesher::nPhys3D() const
   {
      throw std::logic_error("Is not valid in 2D Mesher");
   }

   int I2DMesher::nSpec3D() const
   {
      throw std::logic_error("Is not valid in 2D Mesher");
   }

   int I2DMesher::nDealias3D() const
   {
      throw std::logic_error("Is not valid in 2D Mesher");
   }

} // SpatialScheme
} // QuICC
