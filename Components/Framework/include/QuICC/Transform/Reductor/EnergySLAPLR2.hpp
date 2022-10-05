/**
 * @file EnergySLAPLR2.hpp
 * @brief Reductor transform operator Reductor::EnergySLAPLR2
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGYSLAPLR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGYSLAPLR2_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Reductor/IRegisterId.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   /**
    * @brief Reductor transform operator Reductor::EnergySLAPLR2
    */
   class EnergySLAPLR2: public IRegisterId<EnergySLAPLR2>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergySLAPLR2();

         friend class IRegisterId<EnergySLAPLR2>;

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

} // Reductor
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGYSLAPLR2_HPP
