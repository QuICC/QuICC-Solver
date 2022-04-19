/**
 * @file Sin_1LaplhDphi.hpp
 * @brief Backward projection operator Sin_1LaplhDphi 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SIN_1LAPLHDPHI_HPP
#define QUICC_TRANSFORM_BACKWARD_SIN_1LAPLHDPHI_HPP

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
    * @brief Backward projection operator Sin_1LaplhDphi
    */
   class Sin_1LaplhDphi: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1LaplhDphi();

         /**
          * @brief Destructor
          */
         virtual ~Sin_1LaplhDphi();

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

#endif // QUICC_TRANSFORM_BACKWARD_SIN_1LAPLHDPHI_HPP
