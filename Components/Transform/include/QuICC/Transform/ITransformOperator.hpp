/**
 * @file ITransformOperator.hpp
 * @brief Interface for a generic transform operator
 */

#ifndef QUICC_TRANSFORM_ITRANSFORMOPERATOR_HPP
#define QUICC_TRANSFORM_ITRANSFORMOPERATOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Interface for a generic transform operator
    */
   class ITransformOperator
   {
      public:
         /**
          * @brief Constructor
          */
         ITransformOperator();

         /**
          * @brief Destructor
          */
         virtual ~ITransformOperator();

         /**
          * @brief
          */
         bool isInitialized() const;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

      protected:
         /**
          * @brief Get operator name from typeid
          */
         std::string opName() const;

         /**
          * @brief Set profiling tag based on operator name
          */
         void setProfileTag();

         /**
          * @brief Profiling tag
          */
         std::string mProfileTag;

         /**
          * @brief Need initialization?
          */
         mutable bool mIsInitialized;

      private:
   };

}
}

#endif // QUICC_TRANSFORM_ITRANSFORMOPERATOR_HPP
