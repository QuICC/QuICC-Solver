/**
 * @file Pm.hpp
 * @brief Forward projection operator Pm 
 */

#ifndef QUICC_TRANSFORM_FORWARD_PM_HPP
#define QUICC_TRANSFORM_FORWARD_PM_HPP

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
    * @brief Forward projection operator Pm
    */
   class Pm: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Pm();

         /**
          * @brief Destructor
          */
         virtual ~Pm();

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

#endif // QUICC_TRANSFORM_FORWARD_PM_HPP
