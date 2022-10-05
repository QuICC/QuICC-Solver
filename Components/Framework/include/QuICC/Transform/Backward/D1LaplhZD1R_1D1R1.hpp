/**
 * @file D1LaplhZD1R_1D1R1.hpp
 * @brief Backward transform operator Backard::D1LaplhZD1R_1D1R1
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1LAPLHZD1R_1D1R1_HPP
#define QUICC_TRANSFORM_BACKWARD_D1LAPLHZD1R_1D1R1_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Backward/IRegisterId.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   /**
    * @brief Backward transform operator Backard::D1LaplhZD1R_1D1R1
    */
   class D1LaplhZD1R_1D1R1: public IRegisterId<D1LaplhZD1R_1D1R1>
   {
      public:
         /**
          * @brief Constructor
          */
         D1LaplhZD1R_1D1R1();

         friend class IRegisterId<D1LaplhZD1R_1D1R1>;

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();

         /**
          * @brief Formatted name
          */
         static std::string sFormatted();
   };

} // Backward
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_BACKWARD_D1LAPLHZD1R_1D1R1_HPP
