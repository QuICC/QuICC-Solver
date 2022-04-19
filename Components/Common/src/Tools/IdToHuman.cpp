/**
 * @file IdToHuman.cpp
 * @brief Source of enum id to human strings converters
 */
// Configuration includes
//

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Tools/IdToHuman.hpp"

// Project includes
//

namespace QuICC {

namespace Tools {

   std::string IdToHuman::toString(const FieldComponents::Physical::Id id)
   {
      switch(id)
      {
         case FieldComponents::Physical::X:
            return "X";
         case FieldComponents::Physical::Y:
            return "Y";
         case FieldComponents::Physical::Z:
            return "Z";
         case FieldComponents::Physical::R:
            return "R";
         case FieldComponents::Physical::THETA:
            return "Theta";
         case FieldComponents::Physical::PHI:
            return "Phi";
         case FieldComponents::Physical::SCALAR:
            return "";
         default:
            throw std::logic_error("Unknown ID to string conversion requested (FieldComponents::Physical)");
      }
   }

   std::string IdToHuman::toTag(const FieldComponents::Physical::Id id)
   {
      switch(id)
      {
         case FieldComponents::Physical::X:
            return "x";
         case FieldComponents::Physical::Y:
            return "y";
         case FieldComponents::Physical::Z:
            return "z";
         case FieldComponents::Physical::R:
            return "r";
         case FieldComponents::Physical::THETA:
            return "theta";
         case FieldComponents::Physical::PHI:
            return "phi";
         case FieldComponents::Physical::SCALAR:
            return "";
         default:
            throw std::logic_error("Unknown ID to tag conversion requested (FieldComponents::Physical)");
      }
   }

   std::string IdToHuman::toString(const FieldComponents::Spectral::Id id)
   {
      switch(id)
      {
         case FieldComponents::Spectral::X:
            return "X";
         case FieldComponents::Spectral::Y:
            return "Y";
         case FieldComponents::Spectral::Z:
            return "Z";
         case FieldComponents::Spectral::R:
            return "R";
         case FieldComponents::Spectral::THETA:
            return "Theta";
         case FieldComponents::Spectral::PHI:
            return "Phi";
         case FieldComponents::Spectral::TOR:
            return "Toroidal";
         case FieldComponents::Spectral::POL:
            return "Poloidal";
         case FieldComponents::Spectral::Q:
            return "Q";
         case FieldComponents::Spectral::S:
            return "S";
         case FieldComponents::Spectral::T:
            return "T";
         case FieldComponents::Spectral::SCALAR:
            return "";
         default:
            throw std::logic_error("Unknown ID to string conversion requested (FieldComponents::Spectral)");
      }
   }

   std::string IdToHuman::toTag(const FieldComponents::Spectral::Id id)
   {
      switch(id)
      {
         case FieldComponents::Spectral::X:
            return "x";
         case FieldComponents::Spectral::Y:
            return "y";
         case FieldComponents::Spectral::Z:
            return "z";
         case FieldComponents::Spectral::R:
            return "r";
         case FieldComponents::Spectral::THETA:
            return "theta";
         case FieldComponents::Spectral::PHI:
            return "phi";
         case FieldComponents::Spectral::TOR:
            return "tor";
         case FieldComponents::Spectral::POL:
            return "pol";
         case FieldComponents::Spectral::Q:
            return "q";
         case FieldComponents::Spectral::S:
            return "s";
         case FieldComponents::Spectral::T:
            return "t";
         case FieldComponents::Spectral::SCALAR:
            return "";
         default:
            throw std::logic_error("Unknown ID to tag conversion requested (FieldComponents::Spectral)");
      }
   }

}
}
