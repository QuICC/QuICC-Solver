/** 
 * @file Tools.hpp
 * @brief Definition of some useful constants and tools for polynomial transforms
 */

#ifndef QUICC_TRANSFORM_POLY_TOOLS_HPP
#define QUICC_TRANSFORM_POLY_TOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

   /**
    * @brief Contains some useful constants and tools for polynomial transforms
    */
   class Tools
   {
      public:
         /**
          * @brief Compute the dealiased size using the 3/2a-rule
          *
          * @param size Size to dealias
          */
         static int dealias(const int size);
         
      protected:

      private:
         /**
          * @brief Standard dealiasing factor (usually 3/2)
          */
         static const MHDFloat STD_DEALIASING;

         /**
          * @brief Empty constructor
          */
         Tools();

         /**
          * @brief Empty Destructor
          */
         ~Tools();

   };

}
}
}

#endif // QUICC_TRANSFORM_POLY_TOOLS_HPP
