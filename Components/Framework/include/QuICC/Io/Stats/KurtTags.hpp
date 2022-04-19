/** 
 * @file KurtTags.hpp
 * @brief Definitions and names use by the average writers
 */

#ifndef QUICC_IO_STATS_KURTTAGS_HPP
#define QUICC_IO_STATS_KURTTAGS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Stats {

   /**
    * @brief Definitions and names use by the energy writers
    */
   class KurtTags
   {
      public:
         /**
          * @brief HEADER part for energy files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for energy files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of energy files
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of energy files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         KurtTags();

         /**
         * @brief Destructor
         */
         ~KurtTags();
   };
}
}
}

#endif // QUICC_IO_STATS_KURTTAGS_HPP
