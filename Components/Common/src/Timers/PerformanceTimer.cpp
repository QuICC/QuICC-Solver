/** 
 * @file PerformanceTimer.cpp
 * @brief Source of the performance timer implementation
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
#include "QuICC/Timers/PerformanceTimer.hpp"

// Project includes
//
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

   bool PerformanceTimer::sDoesIo = true;

   PerformanceTimer::PerformanceTimer(const int nTimer, const std::string& name, const int digits)
      : mcDigits(std::pow(10,digits)), mShowOnDelete(true), mName("Perf timer")
   {
      this->init(nTimer);

      std::string n = name;
      if(n == "")
      {
         n = "Perf timer";
      }
      this->setName(n);
   }

   PerformanceTimer::PerformanceTimer(const std::map<int,std::string>& tag, const std::string& name, const int digits)
      : mcDigits(std::pow(10,digits)), mShowOnDelete(true), mName("Perf timer")
   {
      for(const auto& t: tag)
      {
         this->add(t.first, t.second);
      }

      std::string n = name;
      if(n == "")
      {
         n = "Perf timer";
      }
      this->setName(n);
   }

   PerformanceTimer::~PerformanceTimer()
   {
      this->show();
   }

   void PerformanceTimer::add(const int id, const std::string& tag)
   {
      std::string t = tag;
      if(tag == "")
      {
         t = std::to_string(id);
      }
      this->mTimer.insert(std::make_pair(id,TimerMacro(false)));
      this->mTag.insert(std::make_pair(id, t));
   }

   void PerformanceTimer::setTag(const int id, const std::string& tag)
   {
      this->mTag.at(id) = tag;
   }

   void PerformanceTimer::init(const int nTimer)
   {
      if(nTimer > 0)
      {
         for(int i = 0; i < nTimer; i++)
         {
            this->add(i);
         }
      }
   }

   void PerformanceTimer::setName(const std::string& name)
   {
      this->mName = name;
   }

   void PerformanceTimer::allowIo(const bool doesIo)
   {
      PerformanceTimer::sDoesIo = doesIo;
   }

   void PerformanceTimer::start(const int id)
   {
      if(PerformanceTimer::sDoesIo)
      {
         if(this->mTimer.count(id) == 0)
         {
            this->add(id);
         }
         this->mTimer.at(id).start();
      }
   }

   void PerformanceTimer::stop(const int id)
   {
      if(PerformanceTimer::sDoesIo)
      {
         if(this->mTimer.count(id) > 0)
         {
            this->mTimer.at(id).stop();
         } else
         {
            std::logic_error("Tried to stop nonexistent timer");
         }
      }
   }

   double PerformanceTimer::time(const int id)
   {
      return this->mTimer.at(id).time();
   }

   void PerformanceTimer::show()
   {
      if(PerformanceTimer::sDoesIo)
      {
         Tools::Formatter::printLine(std::cout, '-');
         Tools::Formatter::printCentered(std::cout, this->mName, '*');
         Tools::Formatter::printHalfLine(std::cout, '\\', '/');
         Tools::Formatter::printNewline(std::cout);
         for(const auto& t: this->mTimer)
         {
            std::stringstream ss;
            ss << std::ceil(this->mcDigits*t.second.time())/this->mcDigits;
            std::cout << "    " << this->mTag.at(t.first) << ": " << ss.str() << std::endl;
            Tools::Formatter::printNewline(std::cout);
         }
         Tools::Formatter::printLine(std::cout, '+');
      }

      this->mShowOnDelete = false;
   }

}
