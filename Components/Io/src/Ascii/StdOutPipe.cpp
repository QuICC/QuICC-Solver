/** 
 * @file StdOutPipe.cpp
 * @brief Source of the implementation of the standard output pipe
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "QuICC/Io/Ascii/StdOutPipe.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Ascii {

   const std::string StdOutPipe::HEADER = "StandardOutput";

   const std::string StdOutPipe::TYPE = "Flat";

   const std::string StdOutPipe::VERSION = "1.0";

   const std::string StdOutPipe::BASENAME = "_stdout";

   const std::string StdOutPipe::EXTENSION = "";

   StdOutPipe::StdOutPipe(std::string name)
      : IAsciiWriter(name + StdOutPipe::BASENAME, StdOutPipe::EXTENSION, StdOutPipe::HEADER, StdOutPipe::TYPE, StdOutPipe::VERSION, IAsciiWriter::EXTEND)
   {
   }

   StdOutPipe::~StdOutPipe()
   {
   }

   void StdOutPipe::init()
   {
      // Initialise parent
      IAsciiWriter::init();

      // Backup std::cout buffer
      this->mpCoutBuffer = std::cout.rdbuf();

      // Redirect std::cout to this file
      std::cout.rdbuf(this->mFile.rdbuf());
   }

   void StdOutPipe::writeContent()
   {
      // Nothing has to be written explicitly
      // Once this file is initialise any std::cout will be written to this file
   }

   void StdOutPipe::finalize()
   {
      // Put std::cout buffer back in place
      std::cout.rdbuf(this->mpCoutBuffer);

      // Finalise the parent
      IAsciiWriter::finalize();
   }
}
}
}
