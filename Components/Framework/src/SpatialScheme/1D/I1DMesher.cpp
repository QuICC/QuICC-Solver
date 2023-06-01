/** 
 * @file I1DMesher.cpp
 * @brief Source of the 1D spatial scheme mesher
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SpatialScheme/1D/I1DMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   I1DMesher::I1DMesher(const GridPurpose::Id purpose)
      : IMesher(purpose)
   {
   }

   int I1DMesher::nPhys2D() const
   {
      throw std::logic_error("Is not valid in 1D Mesher");
   }

   int I1DMesher::nPhys3D() const
   {
      throw std::logic_error("Is not valid in 1D Mesher");
   }

   int I1DMesher::nSpec2D() const
   {
      throw std::logic_error("Is not valid in 1D Mesher");
   }

   int I1DMesher::nSpec3D() const
   {
      throw std::logic_error("Is not valid in 1D Mesher");
   }

   int I1DMesher::nDealias2D() const
   {
      throw std::logic_error("Is not valid in 1D Mesher");
   }

   int I1DMesher::nDealias3D() const
   {
      throw std::logic_error("Is not valid in 1D Mesher");
   }

} // SpatialScheme
} // QuICC
