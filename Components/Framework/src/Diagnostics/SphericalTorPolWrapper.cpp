/**
 * @file SphericalTorPolWrapper.cpp
 * @brief Source of the spherical toroidal/poloidal decomposition to velocity wrapper
 */

// Debug includes
//
#include <cassert>

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Diagnostics/SphericalTorPolWrapper.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Diagnostics {

   SphericalTorPolWrapper::SphericalTorPolWrapper(Framework::Selector::VariantSharedVectorVariable spTorPol)
      : mspTorPol(spTorPol)
   {
      if(std::visit([](auto&& p){return (p->dom(0).res().sim().ss().formulation() != VectorFormulation::TORPOL);}, this->mspTorPol))
      {
         throw std::logic_error("Spherical Toroidal/Poloidal wrapper cannot be used");
      }
   }

   SphericalTorPolWrapper::~SphericalTorPolWrapper()
   {
   }

   const Framework::Selector::PhysicalScalarField& SphericalTorPolWrapper::one() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).phys().comp(FieldComponents::Physical::R);}, this->mspTorPol);
   }

   const Framework::Selector::PhysicalScalarField& SphericalTorPolWrapper::two() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).phys().comp(FieldComponents::Physical::THETA);}, this->mspTorPol);
   }

   const Framework::Selector::PhysicalScalarField& SphericalTorPolWrapper::three() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).phys().comp(FieldComponents::Physical::PHI);}, this->mspTorPol);
   }

   const Resolution& SphericalTorPolWrapper::res() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).res();}, this->mspTorPol);
   }

}
}
