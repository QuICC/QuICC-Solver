/**
 * @file Laplh.hpp
 * @brief Forward transform operator Forward::Laplh
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH_HPP

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
    * @brief Forward transform operator Forward::Laplh
    */
   class Laplh: public IRegisterId<Laplh>
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh();

         friend class IRegisterId<Laplh>;

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH_HPP
