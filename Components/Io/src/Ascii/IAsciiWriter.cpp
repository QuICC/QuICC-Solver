/**
 * @file IAsciiWriter.cpp
 * @brief Source of the interfact to a general ASCII writer
 */

// Configuration includes
//

// System includes
//
#include <sstream>
#include <iomanip>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Ascii/IAsciiWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   const int IAsciiWriter::msIDWidth = 4;

   IAsciiWriter::IAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const IAsciiWriter::WriteMode mode)
      : AsciiFile(name, ext, header, type, version), mIsInitialized(false), mCounter(0), mEvery(1), mBaseName(name), mMode(mode)
   {
   }

   IAsciiWriter::~IAsciiWriter()
   {
   }

   void IAsciiWriter::overwriteOutput()
   {
      if(this->mIsInitialized)
      {
         throw std::logic_error("Cannot change output mode of " + this->filename() + " after initialization!");
      } else
      {
         this->mMode = IAsciiWriter::OVERWRITE;
      }
   }

   void IAsciiWriter::numberOutput()
   {
      if(this->mIsInitialized)
      {
         throw std::logic_error("Cannot change output mode of " + this->filename() + " after initialization!");
      } else
      {
         this->mMode = IAsciiWriter::NUMBER;
      }
   }

   void IAsciiWriter::extendOutput()
   {
      if(this->mIsInitialized)
      {
         throw std::logic_error("Cannot change output mode of " + this->filename() + " after initialization!");
      } else
      {
         this->mMode = IAsciiWriter::EXTEND;
      }
   }

   bool IAsciiWriter::isActive() const
   {
      return ((this->mCounter % this->mEvery) == 0);
   }

   void IAsciiWriter::onlyEvery(const int n)
   {
      if(n < 1)
      {
         throw std::logic_error("multiplier needs to be > 0 for onlyEvery");
      }

      this->mEvery = n;
   }

   void IAsciiWriter::init()
   {
      this->mIsInitialized = true;

      if(this->mMode == IAsciiWriter::EXTEND)
      {
         this->create();
      }
   }

   void IAsciiWriter::finalize()
   {
      if(this->mMode == IAsciiWriter::EXTEND)
      {
         this->end();
      }
   }

   void IAsciiWriter::preWrite()
   {
      if(this->mMode == IAsciiWriter::OVERWRITE || this->mMode == IAsciiWriter::NUMBER)
      {
         this->create();
      }
   }

   void IAsciiWriter::write()
   {
      if(this->isActive())
      {
         this->writeContent();
      }

      if(this->mMode == IAsciiWriter::NUMBER)
      {
         // Increment the file counter
         ++this->mCounter;
      }
   }

   void IAsciiWriter::postWrite()
   {
      if(this->mMode == IAsciiWriter::OVERWRITE || this->mMode == IAsciiWriter::NUMBER)
      {
         this->end();
      }
   }

   void IAsciiWriter::open()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Get file handle
         this->mFile.open(this->filename().c_str());

         // Check that opening file was successful
         if(! this->mFile.is_open())
         {
            throw std::logic_error("Couldn't open ASCII file " + this->filename() + "!");
         }
      }
   }

   void IAsciiWriter::openDebug()
   {
      // Get file handle
      this->mFile.open(this->filename().c_str());

      // Check that opening file was successful
      if(! this->mFile.is_open())
      {
         throw std::logic_error("Couldn't open ASCII file " + this->filename() + "!");
      }
   }

   void IAsciiWriter::close()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile.close();
      }
   }

   void IAsciiWriter::closeDebug()
   {
      this->mFile.close();
   }

   void IAsciiWriter::create()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         if(this->mMode == IAsciiWriter::NUMBER)
         {
            this->updateName();
         }

         // Create file
         this->open();

         // Add header information
         this->mFile << this->header() << std::endl;
         this->mFile << this->type() << std::endl;
         this->mFile << this->version() << std::endl;
      }
   }

   void IAsciiWriter::end()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Close the file
         this->close();
      }
   }

   void IAsciiWriter::updateName()
   {
      // Create stringstream
      std::ostringstream   oss;

      // Create zerofilled string out of counter value
      oss << std::setfill('0') << std::setw(msIDWidth) << this->mCounter;

      // Reset filename including counter zerofilled value
      this->resetName(this->mBaseName + oss.str());
   }

}
}
}
