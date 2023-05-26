/**
 * @file IRegularSHlBuilder.cpp
 * @brief Source of the Regular basis + Spherical Harmonics scheme implementation with l spectral ordering
 */

// System includes
//
#include <cassert>
#include <set>
#include <stdexcept>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/3D/IRegularSHlBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"
#include "QuICC/SpatialScheme/Tools/UniformSH.hpp"
#include "QuICC/SpatialScheme/Tools/TriangularSH.hpp"
#include "QuICC/SpatialScheme/Tools/TrapezoidalSH.hpp"
#include "QuICC/SpatialScheme/Tools/SpectralTriangularSH.hpp"
#include "QuICC/SpatialScheme/Tools/SpectralTrapezoidalSH.hpp"
#include "QuICC/SpatialScheme/Tools/SHUniform2D.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D3D.hpp"
#include "QuICC/Resolutions/Tools/RegularSHlIndexCounter.hpp"
#include "QuICC/Resolutions/Tools/TriangularSHlIndexCounter.hpp"
#include "QuICC/Resolutions/Tools/TrapezoidalSHlIndexCounter.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"
#include "QuICC/Transform/Setup/Triangular.hpp"
#include "QuICC/Transform/Setup/Trapezoidal.hpp"

namespace QuICC {

namespace SpatialScheme {

   IRegularSHlBuilder::IRegularSHlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : I3DBuilder(dim, purpose, options), mI(dim(0)), mL(dim(1)), mM(dim(2))
   {
      assert(dim.size() == 3);
   }

   ArrayI IRegularSHlBuilder::resolution() const
   {
      ArrayI space(this->dims());
      space << this->mI, this->mL, this->mM;

      return space;
   }

   void IRegularSHlBuilder::addIndexCounter(SharedResolution spRes)
   {
      std::shared_ptr<IndexCounter> spCounter;

      const auto& opt = this->mOptions.at(0);

      // Uniform truncation
      if(std::find(opt.begin(), opt.end(), Transform::Setup::Uniform::id()) != opt.end())
      {
         spCounter = std::make_shared<RegularSHlIndexCounter>(spRes->spSim(), spRes->spCpu());
      }
      // Strict triangular truncation
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::Triangular::id()) != opt.end())
      {
         spCounter = std::make_shared<TriangularSHlIndexCounter>(spRes->spSim(), spRes->spCpu());
      }
      // Trapezoidal truncation
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::Trapezoidal::id()) != opt.end())
      {
         spCounter = std::make_shared<TrapezoidalSHlIndexCounter>(spRes->spSim(), spRes->spCpu());
      }
      else
      {
         throw std::logic_error("Unknown truncation ID");
      }

      spRes->setIndexCounter(spCounter);
   }

   std::shared_ptr<Tools::IBase> IRegularSHlBuilder::truncationTools(const Dimensions::Transform::Id transId) const
   {
      // Setup truncation tools
      std::shared_ptr<Tools::IBase> spTools;

      if(transId == Dimensions::Transform::SPECTRAL)
      {
         const auto& opt = this->mOptions.at(0);

         // Uniform truncation
         if(std::find(opt.begin(), opt.end(), Transform::Setup::Uniform::id()) != opt.end())
         {
            spTools = std::make_shared<Tools::UniformSH>();
         }
         // Strict triangular truncation
         else if(std::find(opt.begin(), opt.end(), Transform::Setup::Triangular::id()) != opt.end())
         {
            spTools = std::make_shared<Tools::SpectralTriangularSH>();
         }
         // Trapezoidal truncation
         else if(std::find(opt.begin(), opt.end(), Transform::Setup::Trapezoidal::id()) != opt.end())
         {
            spTools = std::make_shared<Tools::SpectralTrapezoidalSH>();
         }
         else
         {
            throw std::logic_error("Unknown truncation ID");
         }
      }
      else if(transId == Dimensions::Transform::TRA1D)
      {
         const auto& opt = this->mOptions.at(0);

         // Uniform truncation
         if(std::find(opt.begin(), opt.end(), Transform::Setup::Uniform::id()) != opt.end())
         {
            spTools = std::make_shared<Tools::UniformSH>();
         }
         // Strict triangular truncation
         else if(std::find(opt.begin(), opt.end(), Transform::Setup::Triangular::id()) != opt.end())
         {
            spTools = std::make_shared<Tools::TriangularSH>();
         }
         // Trapezoidal truncation
         else if(std::find(opt.begin(), opt.end(), Transform::Setup::Trapezoidal::id()) != opt.end())
         {
            spTools = std::make_shared<Tools::TrapezoidalSH>();
         }
         else
         {
            throw std::logic_error("Unknown truncation ID");
         }
      }
      else if(transId == Dimensions::Transform::TRA2D)
      {
         spTools = std::make_shared<Tools::SHUniform2D>(this->dim(transId, Dimensions::Data::DATB1D));
      }
      else if(transId == Dimensions::Transform::TRA3D)
      {
         spTools = std::make_shared<Tools::Uniform3D3D>();
      }

      return spTools;
   }

} // SpatialScheme
} // QuICC
