/**
 * @file Transform.cpp
 * @brief Source of the implementation of the Finite Differences Sphere transform
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/FiniteDiff/Sphere/Transform.hpp"
#include "Types/Internal/Math.hpp"
#include "FiniteDiff/Sphere/UniformRadialGrid.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

   Array Transform::generateGrid(const int size)
   {
      // Initialise grid storage
      Internal::Array igrid(size);

      ::QuICC::FiniteDiff::Sphere::UniformRadialGrid quad;
      quad.computeGrid(igrid, size);

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
      if(this->mIGrid.size() == 0)
      {
         throw std::logic_error("FD Sphere transform has not been initialised!");
      }

      return this->mIGrid.cast<MHDFloat>();
   }

   void Transform::initQuadrature()
   {
      if(this->mspSetup->purpose() == GridPurpose::SIMULATION)
      {
         // Set the grid and weights
         ::QuICC::FiniteDiff::Sphere::UniformRadialGrid quad;
         quad.computeGrid(this->mIGrid, this->mspSetup->fwdSize());

      } else if(this->mspSetup->purpose() == GridPurpose::VISUALIZATION)
      {
         // Set the grid and weights
         ::QuICC::FiniteDiff::Sphere::UniformRadialGrid quad;
         quad.computeGrid(this->mIGrid, this->mspSetup->fwdSize());
      }
   }

   void Transform::transform(MatrixZ& rOut, const MatrixZ& in, const IOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup, this->mIGrid);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(Matrix& rOut, const MatrixZ& in, const IOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup, this->mIGrid);
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
         throw std::logic_error("Requested Finite Differences transform operator is not avaible");
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
         throw std::logic_error("Requested Finite Differences transform operator is not avaible");
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

} // Sphere
} // FiniteDiff
} // Transform
} // QuICC
