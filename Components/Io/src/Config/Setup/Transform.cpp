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
         this->iTags().addTag("dim1D", 0);
      }

      if(dim > 1)
      {
         this->iTags().addTag("dim2D", 0);
      }

      if(dim > 2)
      {
         this->iTags().addTag("dim3D", 0);
      }
   }

   void Transform::checkData()
   {
   }

}
}
}
}
