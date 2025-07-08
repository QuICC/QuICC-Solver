/**
 * @file ShellTorPolRSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics radial
 * power calculation for toroidal/poloidal field in a spherical shell
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/ShellTorPolRSpectrumWriter.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

ShellTorPolRSpectrumWriter::ShellTorPolRSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalTorPolRSpectrumWriter(prefix, type)
{}

void ShellTorPolRSpectrumWriter::init()
{
   // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
   MHDFloat ri =
      this->mPhysical.find(NonDimensional::Lower1d::id())->second->value();
   MHDFloat ro =
      this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();
   this->mVolume = (4.0 / 3.0) * Math::PI * (std::pow(ro, 3) - std::pow(ri, 3));

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   // Set the radial grid to match computation
   Internal::Array g, w;
   Polynomial::Quadrature::ChebyshevRule quad;
   int size = 2 * this->res().sim().dim(Dimensions::Simulation::SIM1D,
                     Dimensions::Space::PHYSICAL);
   quad.computeQuadrature(g, w, size, ri, ro);
   this->mGrid.resize(size / 2);
   for (int i = 0; i < size / 2; i++)
   {
      this->mGrid(i) = static_cast<MHDFloat>(g(2 * i));
   }

   ISphericalTorPolRSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
