/**
 * @file Transform.cpp
 * @brief Source of the implementation of the mixed Fourier transform
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Transform.hpp"

// Project includes
//
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

   Array Transform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create equispaced FFT grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = 2.0*Math::PI*static_cast<MHDFloat>(k)/static_cast<MHDFloat>(size);
      }

      return grid;
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
      return Transform::generateGrid(this->mspSetup->fwdSize());
   }

   void Transform::transform(MatrixZ& rOut, const Matrix& in, const Mixed::IMixedOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(Matrix& rOut, const MatrixZ& in, const Mixed::IMixedOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(MatrixZ& rOut, const Matrix& in, const std::size_t id)
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
         throw std::logic_error("Requested Mixed Fourier transform operator is not avaible");
      }
   }

   MHDFloat Transform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);

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
}
