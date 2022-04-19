/**
 * @file CartesianTorPolWrapper.cpp
 * @brief Source of the cartesian toroidal/poloidal decomposition to velocity wrapper
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
#include "QuICC/Diagnostics/CartesianTorPolWrapper.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Diagnostics {

   CartesianTorPolWrapper::CartesianTorPolWrapper(Framework::Selector::VariantSharedVectorVariable spTorPol)
      : mspTorPol(spTorPol)
   {
      if(std::visit([](auto&& p){return (p->dom(0).res().sim().ss().formulation() != VectorFormulation::TORPOL);}, this->mspTorPol))
      {
         throw std::logic_error("Cartesian Toroidal/Poloidal wrapper cannot be used");
      }
   }

   CartesianTorPolWrapper::~CartesianTorPolWrapper()
   {
   }

   const Framework::Selector::PhysicalScalarField& CartesianTorPolWrapper::one() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).phys().comp(FieldComponents::Physical::Z);}, this->mspTorPol);
   }

   const Framework::Selector::PhysicalScalarField& CartesianTorPolWrapper::two() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).phys().comp(FieldComponents::Physical::X);}, this->mspTorPol);
   }

   const Framework::Selector::PhysicalScalarField& CartesianTorPolWrapper::three() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).phys().comp(FieldComponents::Physical::Y);}, this->mspTorPol);
   }

   const Resolution& CartesianTorPolWrapper::res() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspTorPol));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).res();}, this->mspTorPol);
   }

}
}
