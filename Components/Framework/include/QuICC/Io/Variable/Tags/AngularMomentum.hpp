/** 
 * @file AngularMomentum.hpp
 * @brief Definitions and names use by the AngularMomentum writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_ANGULARMOMENTUM
#define QUICC_IO_VARIABLE_TAGS_ANGULARMOMENTUM

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Variable {

namespace Tags {

   /**
    * @brief Definitions and names use by the AngularMomentum writer
    */
   class AngularMomentum
   {
      public:
         /**
          * @brief HEADER part for AngularMomentum file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for AngularMomentum file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of AngularMomentum file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of AngularMomentum file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         AngularMomentum();

         /**
         * @brief Destructor
         */
         ~AngularMomentum();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_ANGULARMOMENTUM
