/** 
 * @file Tools.cpp
 * @brief Source of the utility functions for transform tests
 */

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/TestSuite/Transform/Tools.hpp"

namespace QuICC {

namespace TestSuite {

namespace Transform {

   std::string filenameExt(const FilenameExtDescr& descr)
   {
      std::string s = "";
      for(auto&& e: descr)
      {
         s += "_" + std::get<0>(e);
         if(std::get<2>(e))
         {
            s+= std::to_string(static_cast<int>(std::get<1>(e)));
         } else
         {
            s+= std::to_string(std::get<1>(e));
         }
      }

      return s;
   }

}
}
}
