/** 
 * @file Transform.cpp
 * @brief Source of the implementation of the setup transform node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Setup/Transform.hpp"

// Project includes
//

namespace QuICC {
   
namespace Io {
   
namespace Config {
   
namespace Setup {

   const std::string Transform::PARENTTAG = "transform";

   Transform::Transform(const int dim)
      : IConfigurationNode(Transform::PARENTTAG)
   {
      this->init(dim);
   }

   Transform::~Transform()
   {
   }

   void Transform::init(const int dim)
   {
      if(dim > 0)
      {
         this->sTags().addTag("dim1D", "default");
      }

      if(dim > 1)
      {
         this->sTags().addTag("dim2D", "default");
      }

      if(dim > 2)
      {
         this->sTags().addTag("dim3D", "default");
      }
   }

   void Transform::checkData()
   {
   }

}
}
}
}
