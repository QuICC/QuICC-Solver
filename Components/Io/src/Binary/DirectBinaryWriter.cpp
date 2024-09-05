/**
 * @file DirectBinaryWriter.cpp
 * @brief Source of the direct access binary writer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Binary/DirectBinaryWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Binary {

   DirectBinaryWriter::DirectBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IBinaryWriter(name, ext, header, type, version)
   {
   }

   DirectBinaryWriter::~DirectBinaryWriter()
   {
   }

   void DirectBinaryWriter::init()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Create file
         this->open();
      }
   }

   void DirectBinaryWriter::finalize()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Close the file
         this->close();
      }
   }

   std::ofstream& DirectBinaryWriter::file()
   {
      return this->mFile;
   }

}
}
}
