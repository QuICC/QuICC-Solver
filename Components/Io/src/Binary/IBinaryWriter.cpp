/**
 * @file IBinaryWriter.cpp
 * @brief Source of the general binary writer
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
#include "QuICC/Io/Binary/IBinaryWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Binary {

   IBinaryWriter::IBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : BinaryFile(name, ext, header, type, version)
   {
   }

   IBinaryWriter::~IBinaryWriter()
   {
   }

   void IBinaryWriter::open()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Get file handle
         this->mFile.open(this->filename().c_str(), std::ios::out | std::ios::binary);

         // Check that opening file was successful
         if(! this->mFile.is_open())
         {
            throw std::logic_error("Couldn't open file " + this->filename() + "!");
         }
      }
   }

   void IBinaryWriter::close()
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
