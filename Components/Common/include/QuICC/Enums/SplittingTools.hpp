/**
 * @file SplittingTools.hpp
 * @brief Definition of some useful enums for splitting algorithms 
 */

#ifndef QUICC_SPLITTINGTOOLS_HPP
#define QUICC_SPLITTINGTOOLS_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "QuICC/Enums/Splitting.hpp"

namespace QuICC {

   /**
    * @brief Namespace holding the splitting algorithm related enums
    */
   namespace Splitting {

      /**
       * @brief Convert string to Algorithm ID
       */
      Algorithms::Id getAlgorithmId(const std::string s);

      /**
       * @brief Convert string to Algorithm ID
       */
      Groupers::Id getGrouperId(const std::string s);
   }
}

#endif // QUICC_SPLITTINGTOOLS_HPP
