/** 
 * @file IModelBackend.cpp
 * @brief Source of the interface for a model backend
 */

// First includes
//

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Model/IModelBackend.hpp"

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/NonDimensional/Coordinator.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"

namespace QuICC {

namespace Model {

   IModelBackend::IModelBackend()
      : mUseGalerkin(false), mUseSplitEquation(false), mUseLinearized(false)
   {
   }

   std::vector<std::size_t> IModelBackend::fieldIds() const
   {
      std::vector<std::size_t> ids;
      std::vector<std::string> names = this->fieldNames();
      for(auto it = names.cbegin(); it != names.cend(); ++it)
      {
         ids.push_back(Hasher::makeId(*it));
      }

      return ids;
   }

   std::vector<std::size_t> IModelBackend::paramIds() const
   {
      std::vector<std::size_t> ids;
      std::vector<std::string> names = this->paramNames();
      for(auto it = names.cbegin(); it != names.cend(); ++it)
      {
         ids.push_back(Hasher::makeId(*it));
      }

      return ids;
   }

   std::map<std::string,MHDFloat> IModelBackend::automaticParameters(const std::map<std::string,MHDFloat>& cfg) const
   {
      return std::map<::std::string,MHDFloat>();
   }

   void IModelBackend::checkParamNames(const std::vector<std::string>& names) const
   {
      // Validate list of NonDimensional parameters
      for(auto it = names.cbegin(); it != names.cend(); ++it)
      {
         size_t nd = Hasher::makeId(*it);
         if(NonDimensional::Coordinator::map().count(nd) == 0)
         {
            throw std::logic_error("Requested NonDimensional <"+ *it +"> number is not implemented!");
         }
      }
   }

   void IModelBackend::checkFieldNames(const std::vector<std::string>& names) const
   {
      // Validate list of PhysicalNames
      for(auto it = names.cbegin(); it != names.cend(); ++it)
      {
         size_t pn = Hasher::makeId(*it);
         if(PhysicalNames::Coordinator::map().count(pn) == 0)
         {
            throw std::logic_error("Requested PhysicalName <"+ *it +"> is not implemented!");
         }
      }
   }

   bool IModelBackend::useGalerkin() const
   {
      return this->mUseGalerkin;
   }

   void IModelBackend::enableGalerkin(const bool flag)
   {
      this->mUseGalerkin = flag;
   }

   bool IModelBackend::useSplitEquation() const
   {
      return this->mUseSplitEquation;
   }

   void IModelBackend::enableSplitEquation(const bool tag)
   {
      this->mUseSplitEquation = tag;
   }

   bool IModelBackend::useLinearized() const
   {
      return this->mUseLinearized;
   }

   void IModelBackend::enableLinearized(const bool flag)
   {
      this->mUseLinearized = flag;
   }
}
}
