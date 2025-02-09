/**
 * @file SphereTorPolLSpectrumWriter.cpp
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
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereTorPolLSpectrumWriter::SphereTorPolLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolLSpectrumWriter(prefix, type)
   {
   }

   SphereTorPolLSpectrumWriter::~SphereTorPolLSpectrumWriter()
   {
   }

   void SphereTorPolLSpectrumWriter::init()
   {
      // Sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalTorPolLSpectrumWriter::init();
   }

}
}
}
