/**
 * @file ITransformOperator.hpp
 * @brief Interface for a generic transform operator
 */

#ifndef QUICC_TRANSFORM_ITRANSFORMOPERATOR_HPP
#define QUICC_TRANSFORM_ITRANSFORMOPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Transform/TransformSetup.hpp"


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
         virtual ~ITransformOperator() = default;

         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         virtual void init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const;

         /**
          * @brief Initialise the fft transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         virtual void init(SharedTransformSetup spSetup) const;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const;

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
