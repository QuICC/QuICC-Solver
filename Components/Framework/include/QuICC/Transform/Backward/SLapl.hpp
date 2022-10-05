/**
 * @file SLapl.hpp
 * @brief Backward transform operator Backard::SLapl
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SLAPL_HPP
#define QUICC_TRANSFORM_BACKWARD_SLAPL_HPP

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
    * @brief Backward transform operator Backard::SLapl
    */
   class SLapl: public IRegisterId<SLapl>
   {
      public:
         /**
          * @brief Constructor
          */
         SLapl();

         friend class IRegisterId<SLapl>;

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

#endif // QUICC_TRANSFORM_BACKWARD_SLAPL_HPP
