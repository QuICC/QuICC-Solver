/**
 * @file Transform.cpp
 * @brief Source of the implementation of the Cartesian Chebyshev transform
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Transform.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/NonDimensional/Lower2D.hpp"
#include "QuICC/NonDimensional/Upper2D.hpp"
#include "QuICC/NonDimensional/Lower3D.hpp"
#include "QuICC/NonDimensional/Upper3D.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

   Array Transform::generateGrid(const int size, const MHDFloat xi, const MHDFloat xo)
   {
      if(xo > xi)
      {
         // Initialise grid storage
         Array grid(size);

         // Compute linear map y = ax + b
         MHDFloat b = (xo + xi)/2.0;
         MHDFloat a = (xo - xi)/2.0;

         // Create Chebyshev grid
         for(int k = 0; k < size; k++)
         {
            grid(k) = std::cos((Math::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));

            grid(k) = a*grid(k) + b;
         }

         return grid;
      } else
      {
         throw std::logic_error("generateGrid called with incompatible gap bounds: xi = " + std::to_string(xi) + ", xo = " + std::to_string(xo));
      }
   }

   Transform::Transform()
      : mLower(4242), mUpper(-4242)
   {
   }

   Transform::~Transform()
   {
   }

   void Transform::init(Transform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;
   }

   void Transform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      if(dimId == Dimensions::Transform::TRA1D)
      {
         list.insert(NonDimensional::Lower1D::id());
         list.insert(NonDimensional::Upper1D::id());

      } else if(dimId == Dimensions::Transform::TRA2D)
      {
         list.insert(NonDimensional::Lower2D::id());
         list.insert(NonDimensional::Upper2D::id());

      } else if(dimId == Dimensions::Transform::TRA3D)
      {
         list.insert(NonDimensional::Lower3D::id());
         list.insert(NonDimensional::Upper3D::id());
      }
   }

   void Transform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      if(dimId == Dimensions::Transform::TRA1D)
      {
         assert(options.count(NonDimensional::Lower1D::id()) == 1);
         this->mLower = options.find(NonDimensional::Lower1D::id())->second->value();
         assert(options.count(NonDimensional::Upper1D::id()) == 1);
         this->mUpper = options.find(NonDimensional::Upper1D::id())->second->value();

      } else if(dimId == Dimensions::Transform::TRA2D)
      {
         assert(options.count(NonDimensional::Lower2D::id()) == 1);
         this->mLower = options.find(NonDimensional::Lower2D::id())->second->value();
         assert(options.count(NonDimensional::Upper2D::id()) == 1);
         this->mUpper = options.find(NonDimensional::Upper2D::id())->second->value();

      } else if(dimId == Dimensions::Transform::TRA3D)
      {
         assert(options.count(NonDimensional::Lower3D::id()) == 1);
         this->mLower = options.find(NonDimensional::Lower3D::id())->second->value();
         assert(options.count(NonDimensional::Upper3D::id()) == 1);
         this->mUpper = options.find(NonDimensional::Upper3D::id())->second->value();
      }

      // Set bounds and lock setup
      assert(this->mLower < this->mUpper);
      this->mspSetup->setBounds(this->mLower, this->mUpper);
      this->mspSetup->lock();
   }

   Array Transform::meshGrid() const
   {
      return Transform::generateGrid(this->mspSetup->fwdSize(), this->mLower, this->mUpper);
   }

   void Transform::transform(MatrixZ& rOut, const MatrixZ& in, const IChebyshevOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(Matrix& rOut, const MatrixZ& in, const IChebyshevOperator& op)
   {
      if(!op.isInitialized())
      {
         op.init(this->mspSetup);
      }

      op.transform(rOut, in);
   }

   void Transform::transform(Matrix& rOut, const Matrix& in, const IChebyshevOperator& op)
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
         throw std::logic_error("Requested Chebyshev LinearMap transform operator is not avaible");
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
         throw std::logic_error("Requested Chebyshev LinearMap transform operator is not avaible");
      }
   }

   void Transform::transform(Matrix& rOut, const Matrix& in, const std::size_t id)
   {
      auto it = this->mOps.find(id);

      if(it != this->mOps.end())
      {
         this->transform(rOut, in, *(it->second));
      } else
      {
         throw std::logic_error("Requested Chebyshev LinearMap transform operator is not avaible");
      }
   }

   MHDFloat Transform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
      mem += static_cast<MHDFloat>(2*Debug::MemorySize<MHDFloat>::BYTES);

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
