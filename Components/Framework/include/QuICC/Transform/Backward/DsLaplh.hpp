/**
 * @file DsLaplh.hpp
 * @brief Backward projection operator DsLaplh 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_DSLAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_DSLAPLH_HPP

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
    * @brief Backward projection operator DsLaplh
    */
   class DsLaplh: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         DsLaplh();

         /**
          * @brief Destructor
          */
         virtual ~DsLaplh();

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

#endif // QUICC_TRANSFORM_BACKWARD_DSLAPLH_HPP
