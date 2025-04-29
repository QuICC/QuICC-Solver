/**
 * @file IAsciiReader.cpp
 * @brief Source of the interface to a general ASCII reader
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Ascii/IAsciiReader.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   IAsciiReader::IAsciiReader(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : AsciiFile(name, ext, header, type, version)
   {
   }

   IAsciiReader::~IAsciiReader()
   {
   }

   void IAsciiReader::init()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // open the file
         this->open();

         // Check that the file is compatible with the reader
         this->checkCompatibility();
      }
   }

   void IAsciiReader::open()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Get handle to file
         this->mFile.open(this->filename().c_str());

         // Check that opening was a success
         if(! this->mFile.is_open())
         {
            throw std::logic_error("Couldn't open ASCII file " + this->filename() + "!");
         }
      }
   }

   void IAsciiReader::checkCompatibility()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         std::string fileHead;
         std::string fileType;
         std::string fileVers;

         // Read header if present
         getline(this->mFile, fileHead);

         // Check reading header was successful
         if(this->mFile.good())
         {
            // Check header
            if(fileHead.compare(this->header()) == 0)
            {
               // Read type if present
               getline(this->mFile, fileType);

               // Check reading type was successful
               if(this->mFile.good())
               {
                  // Check type
                  if(fileType.compare(this->type()) == 0)
                  {
                     // Read version if present
                     getline(this->mFile, fileVers);

                     // Check that reading version was successful
                     if(this->mFile.good())
                     {
                        // Check version
                        if(fileVers.compare(this->version()) == 0)
                        {
                           //
                           // Compatibility check was successful
                           //
                        } else
                        {
                           throw std::logic_error("Wrong ASCII file version!");
                        }
                     } else
                     {
                        throw std::logic_error("Missing ASCII file version!");
                     }
                  } else
                  {
                     throw std::logic_error("Wrong ASCII file type!");
                  }
               } else
               {
                  throw std::logic_error("Missing ASCII file type!");
               }
            } else
            {
               throw std::logic_error("Wrong ASCII file header!");
            }
         } else
         {
            throw std::logic_error("Missing ASCII file header!");
         }
      }
   }

   void IAsciiReader::finalize()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Close file
         this->close();
      }
   }

   void IAsciiReader::close()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile.close();
      }
   }

}
}
}
