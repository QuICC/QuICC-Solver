/**
 * @file HumanToId.cpp
 * @brief Source of human strings to enum id converters
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
#include "QuICC/Tools/HumanToId.hpp"

// Project includes
//
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Tools {

   FieldComponents::Spectral::Id HumanToId::toComp(const std::string& id)
   {
      if(id == IdToHuman::toTag(FieldComponents::Spectral::SCALAR))
      {
         return FieldComponents::Spectral::SCALAR;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::X))
      {
         return FieldComponents::Spectral::X;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::Y))
      {
         return FieldComponents::Spectral::Y;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::Z))
      {
         return FieldComponents::Spectral::Z;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::R))
      {
         return FieldComponents::Spectral::R;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::THETA))
      {
         return FieldComponents::Spectral::THETA;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::PHI))
      {
         return FieldComponents::Spectral::PHI;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::Q))
      {
         return FieldComponents::Spectral::Q;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::S))
      {
         return FieldComponents::Spectral::S;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::T))
      {
         return FieldComponents::Spectral::T;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::TOR))
      {
         return FieldComponents::Spectral::TOR;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::POL))
      {
         return FieldComponents::Spectral::POL;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::NOTUSED))
      {
         return FieldComponents::Spectral::NOTUSED;
      } else
      {
         throw std::logic_error("Unknown string to ID conversion requested (FieldComponents::Spectral)");
      }
   }
}
}
