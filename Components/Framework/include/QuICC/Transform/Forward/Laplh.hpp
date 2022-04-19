/**
 * @file Laplh.hpp
 * @brief Forward projection operator Laplh 
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Forward/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   /**
    * @brief Forward projection operator Laplh
    */
   class Laplh: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh();

         /**
          * @brief Destructor
          */
         virtual ~Laplh();

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH_HPP
