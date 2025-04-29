/**
 * @file SphereScalarNSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics L power spectrum calculation for scalar field in a sphere
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
#include "QuICC/Io/Variable/SphereScalarNSpectrumWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereScalarNSpectrumWriter::SphereScalarNSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarNSpectrumWriter(prefix, type)
   {
   }

   SphereScalarNSpectrumWriter::~SphereScalarNSpectrumWriter()
   {
   }

   void SphereScalarNSpectrumWriter::init()
   {
      // Normalize by sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalScalarNSpectrumWriter::init();
   }

}
}
}
