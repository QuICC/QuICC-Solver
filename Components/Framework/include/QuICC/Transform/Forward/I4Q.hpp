/**
 * @file I4Q.hpp
 * @brief Forward transform operator Forward::I4Q
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4Q_HPP
#define QUICC_TRANSFORM_FORWARD_I4Q_HPP

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
    * @brief Forward transform operator Forward::I4Q
    */
   class I4Q: public IRegisterId<I4Q>
   {
      public:
         /**
          * @brief Constructor
          */
         I4Q();

         friend class IRegisterId<I4Q>;

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

#endif // QUICC_TRANSFORM_FORWARD_I4Q_HPP
