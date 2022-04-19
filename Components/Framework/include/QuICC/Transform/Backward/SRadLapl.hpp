/**
 * @file SRadLapl.hpp
 * @brief Backward projection operator SRadLapl 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SRADLAPL_HPP
#define QUICC_TRANSFORM_BACKWARD_SRADLAPL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Backward/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   /**
    * @brief Backward projection operator SRadLapl
    */
   class SRadLapl: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         SRadLapl();

         /**
          * @brief Destructor
          */
         virtual ~SRadLapl();

         /**
          * @brief Unique id
          */
         static const std::size_t& id();
      
      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();
   };

}
}
}

#endif // QUICC_TRANSFORM_BACKWARD_SRADLAPL_HPP
