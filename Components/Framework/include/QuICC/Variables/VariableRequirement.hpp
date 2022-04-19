/**
 * @file VariableRequirement.hpp
 * @brief Implementation of a class to store requirements for the variables
 */

#ifndef QUICC_VARIABLEREQUIREMENT_HPP
#define QUICC_VARIABLEREQUIREMENT_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Variables/FieldRequirement.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a class to store requirements for the variables
    */
   class VariableRequirement
   {
      public:
         /// Typedef for the const iterator
         typedef  std::map<std::size_t,FieldRequirement>::const_iterator   const_iterator;

         /**
          * @brief Constructor
          */
         VariableRequirement();

         /**
          * @brief Destructor
          */
         ~VariableRequirement();

         /**
          * @brief Get field requirements
          */
         const FieldRequirement& field(const std::size_t id) const;

         /**
          * @brief Set field requirements
          */
         FieldRequirement& rField(const std::size_t id);

         /**
          * @brief Add field requirement
          */
         FieldRequirement& addField(const std::size_t id, const FieldRequirement& req);

         /**
          * @brief Const iterator to access all requirements
          */
         const_iterator cbegin() const;

         /**
          * @brief Const iterator to access all requirements
          */
         const_iterator cend() const;

         /**
          * @brief Merge the information of two requirements
          *
          * The operation is basically a OR on the need? calls
          */
         void merge(const VariableRequirement& req);

      protected:

      private:
         /**
          * @brief Not required field (all answers are false)
          */
         FieldRequirement  mNoField;

         /**
          * @brief Storage for all the information
          */
         std::map<std::size_t, FieldRequirement> mInfo;
   };

}

#endif // QUICC_VARIABLEREQUIREMENT_HPP
