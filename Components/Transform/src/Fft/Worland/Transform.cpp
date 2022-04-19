/**
 * @file Transform.cpp
 * @brief Source of the implementation of the FFT based Worland transform
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Transform.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tools.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

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
      internal::Array grid;
      Tools::computeGrid(grid, this->mspSetup->fwdSize());
      return grid.cast<MHDFloat>();
   }

   void Transform::transform(MatrixZ& rOut, const MatrixZ& in, const IWorlandOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(Matrix& rOut, const MatrixZ& in, const IWorlandOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup);
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
