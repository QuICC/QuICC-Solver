/**
 * @file VariableRequirement.cpp
 * @brief Source of the variable requirements
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Variables/VariableRequirement.hpp"

// Project includes
//

namespace QuICC {

   VariableRequirement::VariableRequirement()
      : mNoField(false, Tools::ComponentAlias<FieldComponents::Spectral::Id>(),Tools::ComponentAlias<FieldComponents::Physical::Id>())
   {
   }

   VariableRequirement::~VariableRequirement()
   {
   }

   const FieldRequirement& VariableRequirement::field(const std::size_t id) const
   {
      if(this->mInfo.count(id) == 0)
      {
         return this->mNoField;
      } else
      {
         return this->mInfo.find(id)->second;
      }
   }

   FieldRequirement& VariableRequirement::rField(const std::size_t id)
   {
      if(this->mInfo.count(id) == 0)
      {
         throw std::logic_error("Tried to modify requirements of inexistant field!");
      }

      return this->mInfo.find(id)->second;
   }

   FieldRequirement& VariableRequirement::addField(const std::size_t id, const FieldRequirement& req)
   {
      auto f = this->mInfo.insert(std::make_pair(id,req));

      return (f.first)->second;
   }

   VariableRequirement::const_iterator VariableRequirement::cbegin() const
   {
      return this->mInfo.cbegin();
   }

   VariableRequirement::const_iterator VariableRequirement::cend() const
   {
      return this->mInfo.cend();
   }

   void VariableRequirement::merge(const VariableRequirement& req)
   {
      for(auto it = req.cbegin(); it != req.cend(); it++)
      {
         if(this->mInfo.count(it->first) == 0)
         {
            this->mInfo.insert(*it);
         } else
         {
            this->mInfo.find(it->first)->second.merge(it->second);
         }
      }
   }
}
