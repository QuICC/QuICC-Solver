/** 
 * @file CflWriter.cpp 
 * @brief Source of the implementation of the CFL writer
 */

// Configuration includes
//

// System includes
//
#include <iomanip>
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Io/Ascii/CflWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Ascii/CflTags.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   CflWriter::CflWriter()
      : IAsciiWriter(CflTags::NAME, CflTags::EXTENSION, CflTags::HEADER, CflTags::TYPE, CflTags::VERSION, IAsciiWriter::EXTEND), mTime(-1.0), mTimestep(-1.0), mSteps(0.0), mChanged(false), mDt(2,1), mcIoHigh(14), mcIoLow(7), mcIoIntW(7), mNeedFancy(true)
   {
      mDt.setConstant(-1.0);
   }

   CflWriter::~CflWriter()
   {
   }

   void CflWriter::setSimTime(const MHDFloat time, const Matrix& dt, const MHDFloat steps)
   {
      this->mTime = time;

      this->mChanged = (dt(0,0) != this->mDt(0,0) || static_cast<long int>(this->mSteps) % 100 == 0);
      this->mDt = dt;

      this->mSteps = steps;
   }

   void CflWriter::writeContent()
   {
      this->fancyHeader();

      // pre write
      this->preWrite();

      using Tools::Formatter::ioFW;
      using Tools::Formatter::ioIW;

      // Check if the workflow allows IO to be performed
      if(this->mChanged && QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(this->mcIoHigh);
         this->mFile << std::left << ioFW(this->mcIoHigh) << this->mTime << "\t";

         this->mFile << std::left << ioFW(this->mcIoHigh) << this->mDt(0,0) << "\t";
         for(int j = 1; j < this->mDt.rows(); j++)
         {
            this->mFile << std::left << ioFW(this->mcIoHigh);
            if(this->mDt(j,0) == -1)
            {
               this->mFile << std::setfill(' ') << std::string(3, ' ') + "global";
            } else if(this->mDt(j,0) == -100)
            {
               this->mFile << std::setfill(' ') << std::string(3, ' ') + "fixed";
            } else if(this->mDt(j,0) == -101)
            {
               this->mFile<< std::setfill(' ') << std::string(3, ' ') + "max";
            } else if(this->mDt(j,0) == -102)
            {
               this->mFile << std::setfill(' ') << std::string(3, ' ') + "min";
            } else
            {
               this->mFile << this->mDt(j,0);
            }
            this->mFile << "\t";
         }
         this->mFile << std::setprecision(this->mcIoLow);
         for(int i = 1; i < this->mDt.cols(); i++)
         {
            this->mFile << std::left << ioFW(this->mcIoLow) << this->mDt(0,i) << "\t";
            for(int j = 1; j < this->mDt.rows(); j++)
            {
               this->mFile << std::left << ioFW(this->mcIoLow);
               if(this->mDt(j,i) == -1)
               {
                  this->mFile << std::setfill(' ') << std::string(3, ' ') + "global";
               } else
               {
                  this->mFile << this->mDt(j,i);
               }
               this->mFile << "\t";
            }
         }
         this->mFile << std::right << ioIW(this->mcIoIntW) << this->mSteps << std::endl;
      }

      // post write
      this->postWrite();
   }

   void CflWriter::fancyHeader()
   {
      if(this->mNeedFancy)
      {
         using Tools::Formatter::ioFW;
         using Tools::Formatter::ioIW;

         this->mFile << std::left << std::setfill(' ') << ioFW(this->mcIoHigh) << "# time" << "\t";
         this->mFile << ioFW(this->mcIoHigh) << std::string(2, ' ') + "dt" << "\t";
         this->mFile << ioFW(this->mcIoHigh) << std::string(5, ' ') + "loc" << "\t";

         std::stringstream ss;
         for(int i = 1; i < this->mDt.cols(); i++)
         {
            std::stringstream ss;
            ss << std::string(2, ' ') + "dt_" << i;
            this->mFile << ioFW(this->mcIoLow) << ss.str() << "\t";
            ss.str("");
            ss << std::string(3, ' ') + "loc_" << i;
            this->mFile << ioFW(this->mcIoLow) << ss.str() << "\t";
         }

         this->mFile << std::right << ioIW(this->mcIoIntW) << "steps";
         this->mFile << std::endl;

         this->mNeedFancy = false;
      }
   }

}
}
}
