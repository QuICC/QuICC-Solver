/**
 * @file Elevator.hpp
 * @brief Elevator number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_ELEVATOR_HPP
#define QUICC_NONDIMENSIONAL_ELEVATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/BasicTypes.hpp"
#include "QuICC/NonDimensional/IRegisterId.hpp"

namespace QuICC {

namespace NonDimensional {

   /**
    * @brief Elevator number nondimensional number
    */
   class Elevator: public IRegisterId<Elevator>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Elevator number
          */
         Elevator(const MHDFloat value);

         friend class IRegisterId<Elevator>;

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

#endif // QUICC_NONDIMENSIONAL_ELEVATOR_HPP
