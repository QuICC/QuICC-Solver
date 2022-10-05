/**
 * @file Q.hpp
 * @brief Forward transform operator Forward::Q
 */

#ifndef QUICC_TRANSFORM_FORWARD_Q_HPP
#define QUICC_TRANSFORM_FORWARD_Q_HPP

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
    * @brief Forward transform operator Forward::Q
    */
   class Q: public IRegisterId<Q>
   {
      public:
         /**
          * @brief Constructor
          */
         Q();

         friend class IRegisterId<Q>;

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

#endif // QUICC_TRANSFORM_FORWARD_Q_HPP
