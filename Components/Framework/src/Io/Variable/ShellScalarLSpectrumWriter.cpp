/**
 * @file ShellScalarLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics L energy spectrum calculation for scalar field in a spherical shell
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
#include "QuICC/Io/Variable/ShellScalarLSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ShellScalarLSpectrumWriter::ShellScalarLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarLSpectrumWriter(prefix, type)
   {
   }

   ShellScalarLSpectrumWriter::~ShellScalarLSpectrumWriter()
   {
   }

   void ShellScalarLSpectrumWriter::init()
   {
      // Normalize by spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ri = this->mPhysical.find(NonDimensional::Lower1D::id())->second->value();
      MHDFloat ro = this->mPhysical.find(NonDimensional::Upper1D::id())->second->value();
      this->mVolume = (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalScalarLSpectrumWriter::init();
   }

}
}
}
