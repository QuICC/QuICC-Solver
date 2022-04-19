/**
 * @file Streamfunction.hpp
 * @brief Streamfunction physical name 
 */

#ifndef QUICC_PHYSICALNAMES_STREAMFUNCTION_HPP
#define QUICC_PHYSICALNAMES_STREAMFUNCTION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalNames/IRegisterId.hpp"

namespace QuICC {

namespace PhysicalNames {

   /**
    * @brief Streamfunction physical name
    */
   class Streamfunction: public IRegisterId<Streamfunction>
   {
      public:
         /**
          * @brief Constructor
          */
         Streamfunction();

         friend class IRegisterId<Streamfunction>;
      
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

}
}

#endif // QUICC_PHYSICALNAMES_STREAMFUNCTION_HPP
