/**
 * @file SRadLapl.hpp
 * @brief Backward transform operator Backard::SRadLapl
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SRADLAPL_HPP
#define QUICC_TRANSFORM_BACKWARD_SRADLAPL_HPP

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
    * @brief Backward transform operator Backard::SRadLapl
    */
   class SRadLapl: public IRegisterId<SRadLapl>
   {
      public:
         /**
          * @brief Constructor
          */
         SRadLapl();

         friend class IRegisterId<SRadLapl>;

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

#endif // QUICC_TRANSFORM_BACKWARD_SRADLAPL_HPP
