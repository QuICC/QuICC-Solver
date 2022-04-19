/**
 * @file IOperator.hpp
 * @brief Interface for a generic backward transform operator 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_IOPERATOR_HPP
#define QUICC_TRANSFORM_BACKWARD_IOPERATOR_HPP

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   /**
    * @brief Interface for a generic backward transform operator 
    */
   class IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IOperator(const std::string tag);

         /**
          * @brief Destructor
          */
         virtual ~IOperator();

         /**
          * @brief Tag of the operator
          */
         std::string tag() const;
         
      protected:

      private:
         /**
          * @brief Tag of operator
          */
         const std::string mTag;
   };

   /// Typedef for shared_pointer operator
   typedef std::shared_ptr<IOperator> SharedIOperator;

}
}
}

#endif // QUICC_TRANSFORM_BACKWARD_IOPERATOR_HPP
