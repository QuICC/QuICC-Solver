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

   // Store description in debug mode
   #ifdef QUICC_DEBUG
      bool SplittingDescription::store = true;
   #else
      bool SplittingDescription::store = false;
   #endif //QUICC_DEBUG

   std::string SplittingDescription::nameStage(const int stageId) const
   {
      std::string name = "Data_distribution";

      switch(stageId)
      {
         case 0:
            name += "_TRAB1D";
            break;
         case 1:
            name += "_TRAB2D";
            break;
         case 2:
            name += "_TRAB3D";
            break;
         case 3:
            name += "_SPECTRAL";
            break;
         default:
            name += "_UNKNOWN";
            break;
      }

      return name;
   }

   void SplittingDescription::addStage(const int stageId, const SharedTransformResolution spTRes, const int cpuId)
   {
      // Only store information if enabled
      if(SplittingDescription::store)
      {
         if(this->vtpFiles.count(stageId) == 0)
         {
            auto fname = this->nameStage(stageId);
            auto pVtp = std::make_shared<Io::Xml::VtpWriter>(fname);
            pVtp->init();
            this->vtpFiles.emplace(stageId, pVtp);
         }

         this->vtpFiles.at(stageId)->representResolution(spTRes, cpuId);
      }
   }

   void SplittingDescription::save()
   {
      for(auto&& spF: this->vtpFiles)
      {
         spF.second->write();
         spF.second->finalize();
      }
   }

}
}
