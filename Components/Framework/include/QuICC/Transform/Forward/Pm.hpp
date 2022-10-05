/**
 * @file Pm.hpp
 * @brief Forward transform operator Forward::Pm
 */

#ifndef QUICC_TRANSFORM_FORWARD_PM_HPP
#define QUICC_TRANSFORM_FORWARD_PM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Forward/IRegisterId.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   /**
    * @brief Forward transform operator Forward::Pm
    */
   class Pm: public IRegisterId<Pm>
   {
      public:
         /**
          * @brief Constructor
          */
         Pm();

         friend class IRegisterId<Pm>;

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

} // Forward
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_FORWARD_PM_HPP
