/**
 * @file ShellTorPolNCoeffsPowerWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics N spectra calculation for Tor and Pol scalar fields in a spherical shell
 */

// Configuration includes
//

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Io/Variable/ShellTorPolNCoeffsPowerWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"

namespace QuICC {

namespace Io {

namespace Variable {
   ShellTorPolNCoeffsPowerWriter::ShellTorPolNCoeffsPowerWriter(const std::string& prefix, const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> pF)
      : ISphericalTorPolNScalarSpectrumWriter(prefix, type, std::vector<std::shared_ptr<QuICC::DenseSM::IGenericProfile>>(pF.begin(), pF.end()))
   {
   }

   ShellTorPolNCoeffsPowerWriter::~ShellTorPolNCoeffsPowerWriter() = default;

   void ShellTorPolNCoeffsPowerWriter::init()
   {
      // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ri = this->mPhysical.find(NonDimensional::Lower1d::id())->second->value();
      MHDFloat ro = this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();
      this->mVolume = (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalTorPolNScalarSpectrumWriter::init();
   }

}
}
}
