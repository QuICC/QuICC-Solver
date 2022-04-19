/** 
 * @file StageTimer.cpp
 * @brief Source of the stage timer implementation
 */

// Configuration includes
//

// System includes
//
#include <iostream>
#include <cmath>

// External includes
//

// Class include
//
#include "QuICC/Timers/StageTimer.hpp"

// Project includes
//
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

   bool StageTimer::sDoesIo = true;

   StageTimer::StageTimer(const int digits)
      : mcDigits(std::pow(10,digits)), mLevel(0), mTimer(false)
   {
   }

   StageTimer::~StageTimer()
   {
   }

   void StageTimer::allowIo(const bool doesIo)
   {
      StageTimer::sDoesIo = doesIo;
   }

   void StageTimer::stage(const std::string& msg)
   {
      if(StageTimer::sDoesIo)
      {
         Tools::Formatter::printHalfLine(std::cout, '/', '\\');
         Tools::Formatter::printCentered(std::cout, "... " + msg + " ...", '*');
         Tools::Formatter::printLine(std::cout, '-');
         Tools::Formatter::printNewline(std::cout);
      }
   }

   void StageTimer::completed(const std::string& msg)
   {
      if(StageTimer::sDoesIo)
      {
         Tools::Formatter::printLine(std::cout, '-');
         Tools::Formatter::printCentered(std::cout, msg, '*');
         Tools::Formatter::printHalfLine(std::cout, '\\', '/');
         Tools::Formatter::printNewline(std::cout);
      }
   }

   void StageTimer::msg(const std::string& msg, const int space)
   {
      if(StageTimer::sDoesIo)
      {
         std::cout << std::setfill(' ') << std::setw(space) << "" << msg << std::endl;
      }
   }

   void StageTimer::start(const std::string& msg, const int level)
   {
      this->mLevel = level;

      if(StageTimer::sDoesIo)
      {
         if(this->mLevel == 0)
         {
            StageTimer::msg( "(-- " + msg + " --)", 4 + this->mLevel*4);
         } else
         {
            StageTimer::msg( "- " + msg, 4 + this->mLevel*4);
         }
         this->mTimer.start();
      }
   }

   void StageTimer::done()
   {
      if(StageTimer::sDoesIo)
      {
         this->mTimer.stop();

         std::stringstream ss;
         ss << std::ceil(this->mcDigits*this->mTimer.time())/this->mcDigits;
         StageTimer::msg("    done (" + ss.str() + " s)", 4 + this->mLevel*4);
         Tools::Formatter::printNewline(std::cout);
      }
   }

}
