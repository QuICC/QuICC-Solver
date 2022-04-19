/**
 * @file FastMean.hpp
 * @brief fast mean number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_FASTMEAN_HPP
#define QUICC_NONDIMENSIONAL_FASTMEAN_HPP

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
    * @brief fast mean number nondimensional number
    */
   class FastMean: public IRegisterId<FastMean>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of fast mean number
          */
         FastMean(const MHDFloat value);

         friend class IRegisterId<FastMean>;

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

#endif // QUICC_NONDIMENSIONAL_FASTMEAN_HPP
