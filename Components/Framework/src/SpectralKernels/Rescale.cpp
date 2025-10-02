/**
 * @file Rescale.cpp
 * @brief Source of converter from Worland basis to another Worland type basis
 */

// System includes
//

// Class include
//
#include "QuICC/SpectralKernels/Rescale.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   Rescale::Rescale(const bool isComplex)
      : ISpectralKernel(isComplex)
   {
   }

   void Rescale::init(const FieldComponents::Spectral::Id comp, const MHDFloat scale)
   {
      this->mComp = comp;
      this->mScale = scale;
   }

   MHDVariant Rescale::compute(const int i, const int j, const int k) const
   {
      if(this->mIsComplex)
      {
         return MHDComplex(0,0);
      }
      else
      {
         return 0.0;
      }
   }

   void Rescale::apply(const std::size_t timeId)
   {
      if(timeId == SolveTiming::After::id())
      {
         const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);

         const auto hasMOrdering= this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

         auto applyImpl = [&](auto&& comp)
         {
            if(hasMOrdering)
            {
               // Loop over harmonic order m
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
               {
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     for(int i = 0; i < tRes.dim<Dimensions::Data::DATB1D>(j,k); i++)
                     {
                        auto val = this->mScale*comp.point(i, j, k);
                        comp.setPoint(val, i, j, k);
                     }
                  }
               }
            }
            else
            {
               throw std::logic_error("Not yet implemented");
            }
         };

         if(this->mScalars.size() == 1)
         {
            auto& s = this->mScalars.begin()->second;
            std::visit(
                  [&](auto&& field)
                  {
                     auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                     applyImpl(comp);
                  }, s);
         }
         else if(this->mVectors.size() == 1)
         {
            auto& v = this->mVectors.begin()->second;
            std::visit(
                  [&](auto&& field)
                  {
                     auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                     applyImpl(comp);
                  }, v);
         }
         else
         {
            throw std::logic_error("Missing field");
         }
      }
   }

} // Kernel
} // Spectral
} // QuICC
