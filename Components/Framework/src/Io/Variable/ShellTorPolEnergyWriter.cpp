/**
 * @file ShellTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical shell
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
#include "QuICC/Io/Variable/ShellTorPolEnergyWriter.hpp"

// Project includes
//
#include "Types/Math.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ShellTorPolEnergyWriter::ShellTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnergyWriter(prefix, type)
   {
   }

   ShellTorPolEnergyWriter::~ShellTorPolEnergyWriter()
   {
   }

   void ShellTorPolEnergyWriter::init()
   {
      // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ri = this->mPhysical.find(NonDimensional::Lower1d::id())->second->value();
      MHDFloat ro = this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();
      this->mVolume = (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalTorPolEnergyWriter::init();
   }

}
}
}
