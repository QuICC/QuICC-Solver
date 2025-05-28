/**
 * @file FdWLFMesher.cpp
 * @brief Source of the WLF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/FdWLFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Setup/FiniteDiff.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace SpatialScheme {

   FdWLFMesher::FdWLFMesher(const GridPurpose::Id purpose)
      : xLFMesher(purpose)
   {
   }

   void FdWLFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      xLFMesher::init(dims, options);

      int& N = this->mDims.at(0);
      int& nN_ = this->mNdealias;

      // radial spectral resolution
      nN_ = N+1;

      // Get dealiased transform size
      this->mNr = nN_;
   }

} // SpatialScheme
} // QuICC
