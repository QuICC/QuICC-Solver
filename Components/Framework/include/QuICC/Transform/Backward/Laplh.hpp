/**
 * @file Laplh.hpp
 * @brief Backward transform operator Backard::Laplh
 */

#ifndef QUICC_TRANSFORM_BACKWARD_LAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_LAPLH_HPP

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
    * @brief Backward transform operator Backard::Laplh
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

} // Backward
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_BACKWARD_LAPLH_HPP
