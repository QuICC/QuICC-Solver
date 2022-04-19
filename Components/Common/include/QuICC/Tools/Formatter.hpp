/**
 * @file Formatter.hpp
 * @brief Implementation of a few useful output formating static functions
 */

#ifndef QUICC_TOOLS_FORMATTER_HPP
#define QUICC_TOOLS_FORMATTER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <iostream>
#include <iomanip>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Tools {

namespace Formatter {

   /**
    * @brief Base width for stdout output
    */
   static const int TEXT_WIDTH = 50;

   /**
    * @brief Output a newline
    */
   void printNewline(std::ostream& stream);

   /**
    * @brief Output a line of characters
    */
   void printLine(std::ostream& stream, char fill);

   /**
    * @brief Output a line of a combination of two characters
    */
   void printHalfLine(std::ostream& stream, char fillLeft, char fillRight);

   /**
    * @brief Output a centered text with fill characters
    */
   void printCentered(std::ostream& stream, const std::string& text, char fill = ' ', int base = -1);

   /**
    * @brief Output a centered text with fill characters on the left
    */
   void printLeftCentered(std::ostream& stream, const std::string& text, char fill = ' ', int base = -1, bool hasNewLine = true);

   /**
    * @brief Format width for floating point value
    */
   inline auto ioFW(const int prec, const int w = 6)
   {
      return std::setw(prec + w);
   }

   /**
    * @brief Format width for integer value
    */
   inline auto ioIW(const int w = 5)
   {
      return std::setw(w);
   }

} // Formatter
} // Tools
} // QuICC 

#endif // QUICC_TOOLS_FORMATTER_HPP
