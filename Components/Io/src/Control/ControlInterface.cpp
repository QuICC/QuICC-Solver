/**
 * @file ControlInterface.cpp
 * @brief Source of the external control interface implementation
 */

// Configuration includes
//

// System includes
//
#include <iostream>
#include <sstream>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Control/ControlInterface.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/RuntimeStatus/GoOn.hpp"
#include "QuICC/RuntimeStatus/Stop.hpp"
#include "Environment/MpiTypes.hpp"

namespace QuICC {

namespace Io {

namespace Control {

   ControlInterface::ControlInterface(const bool allowedIo)
      : ControlFile(), mAllowedIo(allowedIo), mStatus(RuntimeStatus::GoOn::id())
   {
      this->init();
   }

   ControlInterface::~ControlInterface()
   {
      this->finalize();
   }

   std::size_t  ControlInterface::status() const
   {
      return this->mStatus;
   }

   void ControlInterface::init()
   {
      // Check if the framework allows IO to be performed
      if(this->mAllowedIo)
      {
         // Create file
         this->create();

         // Close file
         this->close();
      }
   }

   void ControlInterface::read()
   {
      // Reset the current status
      this->mStatus = RuntimeStatus::GoOn::id();

      // Check if the framework allows IO to be performed
      if(this->mAllowedIo)
      {
         // Open file
         this->open();

         // Update current status
         this->update();

         // Close file
         this->close();
      }

      // Share the status over all MPI processes
      #ifdef QUICC_MPI
         std::size_t flag = this->mStatus;
         int ierr = MPI_Bcast(&flag, 1, Environment::MpiTypes::type<std::size_t>(), QuICCEnv().ioRank(), MPI_COMM_WORLD);
         QuICCEnv().check(ierr, 555);
         this->mStatus = flag;

         // Synchronize
         QuICCEnv().synchronize();
      #endif // QUICC_MPI
   }

   void ControlInterface::finalize()
   {
      // Check if the framework allows IO to be performed
      if(this->mAllowedIo)
      {
         // Close the file
         this->close();

         // Delete the file
         this->deleteFile();
      }
   }

   void ControlInterface::create()
   {
      // Check if the framework allows IO to be performed
      if(this->mAllowedIo)
      {
         // Create file
         this->mFile.open(this->filename().c_str(), std::fstream::out);

         // Check if creation was successful
         if(! this->mFile.is_open())
         {
            throw std::logic_error("Couldn't create control file " + this->filename() + "!");
         }
      }
   }

   void ControlInterface::open()
   {
      // Open file
      this->mFile.open(filename().c_str());

      // Check if file could be opened
      if(! this->mFile.is_open())
      {
         // File doesn't exits abort run
         this->mStatus = RuntimeStatus::Stop::id();

      // File exist keep going
      } else
      {
         this->mStatus = RuntimeStatus::GoOn::id();
      }
   }

   void ControlInterface::close()
   {
      // Check if the framework allows IO to be performed
      if(this->mAllowedIo)
      {
         this->mFile.close();
      }
   }

   void ControlInterface::deleteFile()
   {
      // Delete the file
      std::remove(this->filename().c_str());
   }

   void ControlInterface::update()
   {
      // Check if status requests stop of simulation
      if(this->mStatus == RuntimeStatus::Stop::id())
      {
         // Check if the framework allows IO to be performed
         if(this->mAllowedIo)
         {
            // Produce a nice looking output to std output
            Tools::Formatter::printLine(std::cout, '#');
            Tools::Formatter::printCentered(std::cout, "User requested abort!", '#');
            Tools::Formatter::printLine(std::cout, '#');
         }
      }
   }

}
}
}
