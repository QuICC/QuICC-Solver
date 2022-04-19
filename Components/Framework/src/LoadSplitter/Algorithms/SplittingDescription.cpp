/** 
 * @file SplittingDescription.cpp
 * @brief Source of the base of the implementation of the load splitting algorithm description
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   SplittingDescription::SplittingDescription()
   {
      #ifdef QUICC_DEBUG
         auto pVtp = std::make_shared<Io::Xml::VtpWriter>("Data_distribution_TRAB1D");
         pVtp->init();
         this->vtpFiles.push_back(pVtp);
         pVtp = std::make_shared<Io::Xml::VtpWriter>("Data_distribution_TRAB2D");
         pVtp->init();
         this->vtpFiles.push_back(pVtp);
         pVtp = std::make_shared<Io::Xml::VtpWriter>("Data_distribution_TRAB3D");
         pVtp->init();
         this->vtpFiles.push_back(pVtp);
      #endif //QUICC_DEBUG
   }

   SplittingDescription::~SplittingDescription()
   {
   }

}
}
