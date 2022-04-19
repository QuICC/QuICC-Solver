/**
 * @file SLapl.hpp
 * @brief Backward projection operator SLapl
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SLAPL_HPP
#define QUICC_TRANSFORM_BACKWARD_SLAPL_HPP

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
    * @brief Backward projection operator SLapl
    */
   class SLapl: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         SLapl();

         /**
          * @brief Destructor
          */
         virtual ~SLapl();

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

#endif // QUICC_TRANSFORM_BACKWARD_SLAPL_HPP
