/** 
 * @file Formatter.cpp
 * @brief Source of the implementation of a few formatting functions
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "QuICC/Tools/Formatter.hpp"

// Project includes
//

namespace QuICC {

namespace Tools {

namespace Formatter {

   void printNewline(std::ostream& stream)
   {
      stream << std::endl;
   }

   void printLine(std::ostream& stream, char fill)
   {
      stream << std::setfill(fill) << std::setw(TEXT_WIDTH) << "" << std::endl;
   }

   void printHalfLine(std::ostream& stream, char fillLeft, char fillRight)
   {
      stream << std::setfill(fillLeft) << std::setw(TEXT_WIDTH/2) << "" << std::setfill(fillRight) << std::setw(TEXT_WIDTH/2) << "" << std::endl;
   }

   void printCentered(std::ostream& stream, const std::string& text, char fill, int base)
   {
      int pad = TEXT_WIDTH - std::max(static_cast<int>(text.length()),base) - 2;

      stream << std::setfill(fill) << std::setw(text.length() + pad/2 + 2) << " " + text + " " << std::setw(pad/2 + pad%2) << "" << std::endl;
   }

   void printLeftCentered(std::ostream& stream, const std::string& text, char fill, int base, bool hasNewline)
   {
      int pad = TEXT_WIDTH - std::max(static_cast<int>(text.length()),base) - 2;

      stream << std::setfill(fill) << std::setw(text.length() + pad/2 + 2) << " " + text;
      if(hasNewline)
      {
         printNewline(stream);
      }
   }
} // Formatter
} // Tools
} // QuICC
