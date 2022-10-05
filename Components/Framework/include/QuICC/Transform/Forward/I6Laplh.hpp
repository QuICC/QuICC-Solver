/**
 * @file I6Laplh.hpp
 * @brief Forward transform operator Forward::I6Laplh
 */

#ifndef QUICC_TRANSFORM_FORWARD_I6LAPLH_HPP
#define QUICC_TRANSFORM_FORWARD_I6LAPLH_HPP

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
    * @brief Forward transform operator Forward::I6Laplh
    */
   class I6Laplh: public IRegisterId<I6Laplh>
   {
      public:
         /**
          * @brief Constructor
          */
         I6Laplh();

         friend class IRegisterId<I6Laplh>;

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

#endif // QUICC_TRANSFORM_FORWARD_I6LAPLH_HPP
