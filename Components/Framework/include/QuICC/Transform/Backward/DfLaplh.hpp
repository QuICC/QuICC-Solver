/**
 * @file DfLaplh.hpp
 * @brief Backward transform operator Backard::DfLaplh
 */

#ifndef QUICC_TRANSFORM_BACKWARD_DFLAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_DFLAPLH_HPP

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
    * @brief Backward transform operator Backard::DfLaplh
    */
   class DfLaplh: public IRegisterId<DfLaplh>
   {
      public:
         /**
          * @brief Constructor
          */
         DfLaplh();

         friend class IRegisterId<DfLaplh>;

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

#endif // QUICC_TRANSFORM_BACKWARD_DFLAPLH_HPP
