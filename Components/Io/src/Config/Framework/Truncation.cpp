/** 
 * @file Truncation.cpp
 * @brief Source of the implementation of the truncation node of the configuration
 */

// Configuration includes
//

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Framework/Truncation.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   const std::string Truncation::PARENTTAG = "truncation";

   Truncation::Truncation(const int dim, const std::vector<bool>& isPeriodicBox)
      : IConfigurationNode(Truncation::PARENTTAG)
   {
      this->init(dim, isPeriodicBox);
   }

   Truncation::~Truncation()
   {
   }

   void Truncation::init(const int dim, const std::vector<bool>& isPeriodicBox)
   {
      // Safety assert
      assert(isPeriodicBox.size() == static_cast<size_t>(dim));

      // Add first dimension truncation tags
      if(dim > 0)
      {
         this->iTags().addTag("dim1D", -1);

         if(isPeriodicBox.at(0))
         {
            this->fTags().addTag("kc1D", -1);
            this->fTags().addTag("box1D", -1);
         }
      }

      // Add second dimension truncation tags
      if(dim > 1)
      {
         this->iTags().addTag("dim2D", -1);

         if(isPeriodicBox.at(1))
         {
            this->fTags().addTag("kc2D", -1);
            this->fTags().addTag("box2D", -1);
         }
      }

      // Add third dimension truncation tags
      if(dim > 2)
      {
         this->iTags().addTag("dim3D", -1);

         if(isPeriodicBox.at(2))
         {
            this->fTags().addTag("kc3D", -1);
            this->fTags().addTag("box3D", -1);
         }
      }
   }

   void Truncation::checkData()
   {
   }

}
}
}
}
