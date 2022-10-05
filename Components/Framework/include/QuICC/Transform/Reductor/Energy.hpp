/**
 * @file Energy.hpp
 * @brief Reductor transform operator Reductor::Energy
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGY_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGY_HPP

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
    * @brief Reductor transform operator Reductor::Energy
    */
   class Energy: public IRegisterId<Energy>
   {
      public:
         /**
          * @brief Constructor
          */
         Energy();

         friend class IRegisterId<Energy>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGY_HPP
