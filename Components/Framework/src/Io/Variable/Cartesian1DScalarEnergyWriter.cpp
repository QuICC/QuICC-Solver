/**
 * @file Cartesian1DScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Chebyshev energy calculation for scalar field in a plane layer
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Io/Variable/Cartesian1DScalarEnergyWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian1DScalarEnergyWriter::Cartesian1DScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : ICartesian1DScalarEnergyWriter(prefix, type)
   {
   }

   void Cartesian1DScalarEnergyWriter::init()
   {
      // Normalize by volume
      const auto zi = this->mPhysical.find(NonDimensional::Lower1d::id())->second->value();
      const auto zo = this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();
      const auto box2D = this->res().sim().boxScale(Dimensions::Simulation::SIM2D);
      const auto box3D = this->res().sim().boxScale(Dimensions::Simulation::SIM3D);
      this->mVolume = (zo - zi)/(box2D*box3D);

      ICartesian1DScalarEnergyWriter::init();
   }

} // Variable
} // Io
} // QuICC
