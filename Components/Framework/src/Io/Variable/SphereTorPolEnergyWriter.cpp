/**
 * @file SphereTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for toroidal/poloidal field in a sphere
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/SphereTorPolEnergyWriter.hpp"

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "Types/Math.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereTorPolEnergyWriter::SphereTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnergyWriter(prefix, type)
   {
   }

   SphereTorPolEnergyWriter::~SphereTorPolEnergyWriter()
   {
   }

   void SphereTorPolEnergyWriter::init()
   {
      // Sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalTorPolEnergyWriter::init();
   }

}
}
}
