/**
 * @file DsLaplh.hpp
 * @brief Backward transform operator Backard::DsLaplh
 */

#ifndef QUICC_TRANSFORM_BACKWARD_DSLAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_DSLAPLH_HPP

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
    * @brief Backward transform operator Backard::DsLaplh
    */
   class DsLaplh: public IRegisterId<DsLaplh>
   {
      public:
         /**
          * @brief Constructor
          */
         DsLaplh();

         friend class IRegisterId<DsLaplh>;

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

#endif // QUICC_TRANSFORM_BACKWARD_DSLAPLH_HPP
