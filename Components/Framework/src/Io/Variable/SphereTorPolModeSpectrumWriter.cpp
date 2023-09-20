/**
 * @file SphereTorPolModeSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for toroidal/poloidal field in a sphere
 */

// Configuration includes
//

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/SphereTorPolModeSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "Types/Constants.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereTorPolModeSpectrumWriter::SphereTorPolModeSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolModeSpectrumWriter(prefix, type)
   {
   }

   void SphereTorPolModeSpectrumWriter::init()
   {
      // Sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

      ISphericalTorPolModeSpectrumWriter::init();
   }

}
}
}
