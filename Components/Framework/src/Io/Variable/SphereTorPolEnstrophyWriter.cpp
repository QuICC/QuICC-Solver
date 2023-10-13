/**
 * @file SphereTorPolEnstrophyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics enstrophy calculation for toroidal/poloidal field in a sphere
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/SphereTorPolEnstrophyWriter.hpp"

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "Types/Math.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereTorPolEnstrophyWriter::SphereTorPolEnstrophyWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyWriter(prefix, type)
   {
   }

   SphereTorPolEnstrophyWriter::~SphereTorPolEnstrophyWriter()
   {
   }

   void SphereTorPolEnstrophyWriter::init()
   {
      // Sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalTorPolEnstrophyWriter::init();
   }

}
}
}
