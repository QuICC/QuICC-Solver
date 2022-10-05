/**
 * @file Pol.hpp
 * @brief Forward transform operator Forward::Pol
 */

#ifndef QUICC_TRANSFORM_FORWARD_POL_HPP
#define QUICC_TRANSFORM_FORWARD_POL_HPP

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
    * @brief Forward transform operator Forward::Pol
    */
   class Pol: public IRegisterId<Pol>
   {
      public:
         /**
          * @brief Constructor
          */
         Pol();

         friend class IRegisterId<Pol>;

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

#endif // QUICC_TRANSFORM_FORWARD_POL_HPP
