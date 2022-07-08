/**
 * @file Pol.hpp
 * @brief Forward projection operator Pol 
 */

#ifndef QUICC_TRANSFORM_FORWARD_POL_HPP
#define QUICC_TRANSFORM_FORWARD_POL_HPP

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
    * @brief Forward projection operator Pol
    */
   class Pol: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Pol();

         /**
          * @brief Destructor
          */
         virtual ~Pol();

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

#endif // QUICC_TRANSFORM_FORWARD_POL_HPP
