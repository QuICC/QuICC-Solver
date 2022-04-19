/** 
 * @file Hasher.cpp
 * @brief Source of the interface to STL hasher
 */

// System includes
//
#include <functional>

// External includes
//

// Class include
//
#include "QuICC/Hasher.hpp"

// Project includes
//

namespace QuICC {

   std::size_t Hasher::makeId(const std::string s)
   {
      std::hash<std::string> hash_fn;
      return hash_fn(s);
   }

   Hasher::Hasher()
   {
   }

   Hasher::~Hasher()
   {
   }

}
