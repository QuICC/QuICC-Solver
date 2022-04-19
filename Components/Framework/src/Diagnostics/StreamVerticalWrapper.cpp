/**
 * @file StreamVerticalWrapper.cpp
 * @brief Source of the streamfunction and vertical velocity to velocity wrapper
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
#include "QuICC/Diagnostics/StreamVerticalWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   StreamVerticalWrapper::StreamVerticalWrapper(Framework::Selector::VariantSharedScalarVariable spStream, Framework::Selector::VariantSharedScalarVariable spVertical)
      : mspStream(spStream), mspVertical(spVertical)
   {
   }

   StreamVerticalWrapper::~StreamVerticalWrapper()
   {
   }

   const Framework::Selector::PhysicalScalarField& StreamVerticalWrapper::one() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspVertical));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).phys();}, this->mspVertical);
   }

   const Framework::Selector::PhysicalScalarField& StreamVerticalWrapper::two() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspStream));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).grad().comp(FieldComponents::Physical::Y);}, this->mspStream);
   }

   const Framework::Selector::PhysicalScalarField& StreamVerticalWrapper::three() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspStream));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).grad().comp(FieldComponents::Physical::X);}, this->mspStream);
   }

   const Resolution& StreamVerticalWrapper::res() const
   {
      // Safety assert
      assert(std::visit([](auto&& p){return (p != nullptr);}, this->mspStream));

      return std::visit([](auto&& p)->auto&& {return p->dom(0).res();}, this->mspStream);
   }

}
}
