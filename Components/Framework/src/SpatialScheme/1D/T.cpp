/**
 * @file T.cpp
 * @brief Source of the Chebyshev(FFT) scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/1D/T.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Coordinator.hpp"
#include "QuICC/SpatialScheme/1D/TBuilder.hpp"
#include "QuICC/Transform/CartesianChebyshevTransform.hpp"
#include "QuICC/Equations/Tools/FullyCoupled.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string T::sTag = "T";

   const std::string T::sFormatted = "T";

   const std::size_t T::sId = registerId<T>(T::sTag);

   T::T(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 1, T::sId, T::sTag, T::sFormatted)
   {
      // Enable basic scheme features
      this->enable(Feature::CartesianGeometry);
      this->enable(Feature::RegularSpectrum);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering123);
      this->enable(Feature::RealSpectrum);
   }

   T::~T()
   {
   }

   std::shared_ptr<IBuilder> T::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<TBuilder>(dim, this->purpose(), needInterpretation);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> T::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            auto spTT = std::make_shared<Transform::CartesianChebyshevTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::CartesianChebyshevTransform::SetupType>(spSetup);
            if(!spST)
            {
               throw std::logic_error("Incompatible transform setup given");
            }
            spTT->init(spST);
            spTransform = spTT;
            break;
         }
         default:
            throw std::logic_error("Tried to initialize too many transforms");
      }

      return spTransform;
   }

   std::shared_ptr<Parallel::IIndexConv> T::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         default:
            throw std::logic_error("Tried to initialize impossible index converter");
      }

      return spConv;
   }

   std::shared_ptr<Equations::Tools::ICoupling> T::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SINGLE)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::FullyCoupled>();
      return spCoupling;
   }

   T::VariantTransformDataPointer T::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<T::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   T::VariantTransformDataPointer T::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<T::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   T::ScalarVariable T::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      T::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   T::VectorVariable T::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      T::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

}
}
