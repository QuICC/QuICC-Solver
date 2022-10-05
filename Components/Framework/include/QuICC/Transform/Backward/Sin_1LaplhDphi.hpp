/**
 * @file Sin_1LaplhDphi.hpp
 * @brief Backward transform operator Backard::Sin_1LaplhDphi
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SIN_1LAPLHDPHI_HPP
#define QUICC_TRANSFORM_BACKWARD_SIN_1LAPLHDPHI_HPP

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
    * @brief Backward transform operator Backard::Sin_1LaplhDphi
    */
   class Sin_1LaplhDphi: public IRegisterId<Sin_1LaplhDphi>
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1LaplhDphi();

         friend class IRegisterId<Sin_1LaplhDphi>;

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

#endif // QUICC_TRANSFORM_BACKWARD_SIN_1LAPLHDPHI_HPP
