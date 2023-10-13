/**
 * @file Transform.cpp
 * @brief Source of the implementation of the Worland transform
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Transform.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

   Array Transform::generateGrid(const int size)
   {
      // Initialise grid storage
      Internal::Array igrid(size);
      Internal::Array itmp(size);

      Polynomial::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, itmp, size);

      return igrid.cast<MHDFloat>();
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

   void Transform::addOperator(const MapFunctor& f)
   {
      f(this->mOps);
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
      if(this->mIGrid.size() == 0 || this->mIWeights.size() == 0)
      {
         throw std::logic_error("Transform has not been initialised!");
      }

      return this->mIGrid.cast<MHDFloat>();
   }

   void Transform::initQuadrature()
   {
      if(this->mspSetup->purpose() == GridPurpose::SIMULATION)
      {
         // Set the grid and weights
         Polynomial::Quadrature::WorlandRule quad;
         quad.computeQuadrature(this->mIGrid, this->mIWeights, this->mspSetup->fwdSize());

      } else if(this->mspSetup->purpose() == GridPurpose::VISUALIZATION)
      {
         Internal::Array iGrid, iWeights;

         // Set the grid and weights
         Polynomial::Quadrature::WorlandRule quad;
         quad.computeQuadrature(iGrid, iWeights, this->mspSetup->fwdSize()-2);

         this->mIGrid.resize(this->mspSetup->fwdSize());
         this->mIWeights = Internal::Array::Zero(this->mspSetup->fwdSize());
         this->mIGrid.segment(1, iGrid.size()) = iGrid;
         this->mIGrid(0) = 0;
         this->mIGrid(iGrid.size()+1) = 1;
         this->mIWeights.segment(1, iWeights.size()) = iWeights;
      }
   }

   void Transform::transform(MatrixZ& rOut, const MatrixZ& in, const IWorlandOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup, this->mIGrid, this->mIWeights);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(Matrix& rOut, const MatrixZ& in, const IWorlandOperator& op)
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
         throw std::logic_error("Requested Worland transform operator is not avaible");
      }
   }

   void Transform::transform(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      auto it = this->mOps.find(id);

      if(it != this->mOps.end())
      {
         this->transform(rOut, in, *(it->second));
      } else
      {
         throw std::logic_error("Requested Worland transform operator is not avaible");
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

      for(auto it = this->mOps.cbegin(); it != this->mOps.cend(); ++it)
      {
         mem += it->second->requiredStorage();
      }
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

} // Worland
} // Poly
} // Transform
} // QuICC
