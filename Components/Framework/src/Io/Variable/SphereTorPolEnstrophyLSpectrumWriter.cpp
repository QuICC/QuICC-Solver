/**
 * @file SphereTorPolEnstrophyLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics enstrophy L spectrum calculation for toroidal/poloidal field in a sphere
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/SphereTorPolEnstrophyLSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereTorPolEnstrophyLSpectrumWriter::SphereTorPolEnstrophyLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyLSpectrumWriter(prefix, type)
   {
   }

   SphereTorPolEnstrophyLSpectrumWriter::~SphereTorPolEnstrophyLSpectrumWriter()
   {
   }

   void SphereTorPolEnstrophyLSpectrumWriter::init()
   {
      // Sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

      ISphericalTorPolEnstrophyLSpectrumWriter::init();
   }

}
}
}
