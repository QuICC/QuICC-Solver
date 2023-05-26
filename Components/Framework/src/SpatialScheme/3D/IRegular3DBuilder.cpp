/** 
 * @file IRegular3DBuilder.cpp
 * @brief Source of the regular 3D scheme builder
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SpatialScheme/3D/IRegular3DBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D1D.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D2D.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D3D.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"

namespace QuICC {

namespace SpatialScheme {

   IRegular3DBuilder::IRegular3DBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : I3DBuilder(dim, purpose, options), mI(dim(0)), mJ(dim(1)), mK(dim(2))
   {
      assert(dim.size() == 3);
   }

   ArrayI IRegular3DBuilder::resolution() const
   {
      ArrayI space(this->dims());
      space << this->mI, this->mJ, this->mK;

      return space;
   }

   std::shared_ptr<Tools::IBase> IRegular3DBuilder::truncationTools(const Dimensions::Transform::Id transId) const
   {
      // Setup truncation tools
      std::shared_ptr<Tools::IBase> spTools;

      if(transId == Dimensions::Transform::TRA1D || transId == Dimensions::Transform::SPECTRAL)
      {
         spTools = std::make_shared<Tools::Uniform3D1D>();
      }
      else if(transId == Dimensions::Transform::TRA2D)
      {
         spTools = std::make_shared<Tools::Uniform3D2D>();
      }
      else if(transId == Dimensions::Transform::TRA3D)
      {
         spTools = std::make_shared<Tools::Uniform3D3D>();
      }

      return spTools;
   }

} // SpatialScheme
} // QuICC
