/** 
 * @file SplittingTools.cpp
 * @brief Source of utility tools for the splitting IDs
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>
#include <algorithm>
#include <string>

// External includes
//

// Class include
//
#include "QuICC/Enums/SplittingTools.hpp"

// Project includes
//

#include <iostream>
namespace QuICC {

namespace Splitting {

   Algorithms::Id getAlgorithmId(const std::string s)
   {
      std::string tag = s;
      std::transform(s.cbegin(), s.cend(), tag.begin(),
                                          [](unsigned char c) { return std::tolower(c); });
      if(tag == "serial")
      {
         return Algorithms::SERIAL;
      }
      else if(tag == "single1d")
      {
         return Algorithms::SINGLE1D;
      }
      else if(tag == "single1d")
      {
         return Algorithms::SINGLE1D;
      }
      else if(tag == "single2d")
      {
         return Algorithms::SINGLE2D;
      }
      else if(tag == "tubular")
      {
         return Algorithms::TUBULAR;
      }
      else if(tag == "coupled2d")
      {
         return Algorithms::COUPLED2D;
      }
      else
      {
         throw std::logic_error("Unknown splitting algorithm: " + s + "(" + tag + ")");
      }
   }

   Groupers::Id getGrouperId(const std::string s)
   {
      std::string tag = s;
      std::transform(s.cbegin(), s.cend(), tag.begin(),
                                          [](unsigned char c) { return std::tolower(c); });
      if(tag == "equation")
      {
         return Groupers::EQUATION;
      }
      else if(tag == "single1d")
      {
         return Groupers::SINGLE1D;
      }
      else if(tag == "single2d")
      {
         return Groupers::SINGLE2D;
      }
      else if(tag == "transform")
      {
         return Groupers::TRANSFORM;
      }
      else
      {
         throw std::logic_error("Unknown communication grouper: " + s + "(" + tag + ")");
      }
   }

}
}
