/**
 * @file EquationParameters.cpp
 * @brief Source of the implementation of the non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/EquationParameters.hpp"

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/NonDimensional/Coordinator.hpp"

namespace QuICC {

namespace Equations {

   EquationParameters::EquationParameters()
   {
   }

   EquationParameters::~EquationParameters()
   {
   }

   MHDFloat EquationParameters::nd(std::size_t id) const
   {
      return this->mND.find(id)->second->value();
   }

   const EquationParameters::NDMapType& EquationParameters::map() const
   {
      return this->mND;
   }

   EquationParameters::NDMapType& EquationParameters::rMap()
   {
      return this->mND;
   }

   std::vector<std::size_t> EquationParameters::ids() const
   {
      // Storage for the IDs
      std::vector<std::size_t> ids;

      for(auto it = this->mND.cbegin(); it != this->mND.cend(); it++)
      {
         ids.push_back(it->first);
      }

      return ids;
   }

   std::vector<std::string> EquationParameters::names() const
   {
      // Storage for the names
      std::vector<std::string> names;

      for(auto it = this->mND.cbegin(); it != this->mND.cend(); it++)
      {
         names.push_back(it->second->tag());
      }

      return names;
   }

   void EquationParameters::init(const std::map<std::string,MHDFloat>& parameters)
   {
      for(auto it = parameters.cbegin(); it != parameters.cend(); ++it)
      {
         std::size_t nd = Hasher::makeId(it->first);

         if(NonDimensional::Coordinator::map().count(nd) > 0)
         {
            this->mND.insert(std::make_pair(nd,NonDimensional::Coordinator::map().find(nd)->second->create(it->second)));
         } else
         {
            throw std::logic_error("NonDimensional number " + it->first + " not found!");
         }
      }
   }

}
}
