/** 
 * @file IStateGeneratorBuilder.hpp
 * @brief Interface for a state generator builder
 */

#ifndef QUICC_MODEL_ISTATEGENERATORBUILDER_HPP
#define QUICC_MODEL_ISTATEGENERATORBUILDER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Model/IPhysicalModel.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for the implementation of a physical models
    */
   template <typename TState> class IStateGeneratorBuilder: public IPhysicalModel
   {
      public:
         /**
          * @brief Constructor
          */
         IStateGeneratorBuilder() = default;

         /**
          * @brief Destructor
          */
         virtual ~IStateGeneratorBuilder() = default;

         /**
          * @brief Add the initial state generation equations
          *
          * @param spGen   Shared generator object
          */
         virtual void addStates(std::shared_ptr<TState> spGen) = 0;

      protected:

      private:
   };

}
}

#endif // QUICC_MODEL_ISTATEGENERATORBUILDER_HPP
