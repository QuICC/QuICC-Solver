/**
 * @file Transform.cpp
 * @brief Source of the implementation of the associated Legendre transform
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Transform.hpp"
#include "Types/Internal/Math.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

   Array Transform::generateGrid(const int size)
   {
      // Initialise grid storage
      Internal::Array igrid(size);
      Internal::Array itmp(size);

      Polynomial::Quadrature::LegendreRule quad;
      quad.computeQuadrature(igrid, itmp, size);

      return igrid.array().acos().matrix().cast<MHDFloat>();
   }

   Transform::Transform()
   {
   }

   Transform::~Transform()
   {
   }

   void Transform::init(Transform::SharedSetupType spSetup)
   {

      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Initialise the quadrature grid and weights
      this->initQuadrature();

      // Lock setup
      this->mspSetup->lock();
   }

   void Transform::requiredOptions(std::set<std::size_t>&, const Dimensions::Transform::Id) const
   {
      //
      // No possible options
      //
   }

   void Transform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>&, const Dimensions::Transform::Id)
   {
      //
      // No possible options
      //
   }

   Array Transform::meshGrid() const
   {
      if(this->mIGrid.size() == 0 || this->mThGrid.size() == 0 || this->mIWeights.size() == 0)
      {
         throw std::logic_error("Transform has not been initialised!");
      }

      return this->mThGrid;
   }

   void Transform::initQuadrature()
   {
      if(this->mspSetup->purpose() == GridPurpose::SIMULATION)
      {
         // Set the grid and weights
         Polynomial::Quadrature::LegendreRule quad;
         quad.computeQuadrature(this->mIGrid, this->mIWeights, this->mspSetup->fwdSize());

      } else if(this->mspSetup->purpose() == GridPurpose::VISUALIZATION)
      {
         Internal::Array iGrid, iWeights;

         // Set the grid and weights
         Polynomial::Quadrature::LegendreRule quad;
         quad.computeQuadrature(iGrid, iWeights, this->mspSetup->fwdSize()-2);

         this->mIGrid.resize(this->mspSetup->fwdSize());
         this->mIWeights = Internal::Array::Zero(this->mspSetup->fwdSize());
         this->mIGrid.segment(1, iGrid.size()) = iGrid;
         this->mIGrid(0) = -1;
         this->mIGrid(iGrid.size()+1) = 1;
         this->mIWeights.segment(1, iWeights.size()) = iWeights;
      }

      // Compute theta grid
      this->mThGrid = this->mIGrid.array().acos().cast<MHDFloat>();

      // Normalise weights by 2*pi for spherical harmonics
      this->mIWeights.array() *= 2.0*Internal::Math::PI;
   }

   void Transform::transform(MatrixZ& rOut, const MatrixZ& in, const ITransformOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup, this->mIGrid, this->mIWeights);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      auto it = this->mOps.find(id);

      if(it != this->mOps.end())
      {
         this->transform(rOut, in, *(it->second));
      } else
      {
         throw std::logic_error("Requested ALegendre transform operator is not available");
      }
   }

   MHDFloat Transform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);

      // Storage for the grid and weight
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mIGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mIWeights.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mThGrid.size();

      for(auto it = this->mOps.cbegin(); it != this->mOps.cend(); ++it)
      {
         mem += it->second->requiredStorage();
      }
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
