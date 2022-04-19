/**
 * @file Transform.hpp
 * @brief Implementation of the Cartesian Chebyshev transform
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_TRANSFORM_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_TRANSFORM_HPP

// Debug includes
//

// Configuration includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <set>
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the Cartsian Chebyshev transform
    */
   class Transform
   {
      public:
         /// Typedef for the configuration class
         typedef IChebyshevOperator::SetupType SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef IChebyshevOperator::SharedSetupType SharedSetupType;

         /**
          * @brief Constructor
          */
         Transform();

         /**
          * @brief Destructor
          */
         virtual ~Transform();

         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Map transform operator to ID
          */
         template <typename TOp> void addOperator(const std::size_t id = 0);

         /**
          * @brief set list of required options
          */
         void requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const;

         /**
          * @brief Set the required options
          */
         void setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId);

         /**
          * @brief Get the physical grid
          */
         Array meshGrid() const;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(MatrixZ& rOut, const MatrixZ& in, const IChebyshevOperator& op);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const MatrixZ& in, const IChebyshevOperator& op);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const Matrix& in, const IChebyshevOperator& op);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(MatrixZ& rOut, const MatrixZ& in, const std::size_t id);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const MatrixZ& in, const std::size_t id);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const Matrix& in, const std::size_t id);

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size, const MHDFloat xi, const MHDFloat xo);

         /**
          * @brief Polynomial setup object providing the sizes
          */
         SharedSetupType    mspSetup;

         /**
          * @brief Lower bound of Cartesian domain
          */
         MHDFloat mLower;

         /**
          * @brief Upper bound of Cartesian domain
          */
         MHDFloat mUpper;

         /**
          * @brief Store transform operator to ID mapping
          */
         std::map<std::size_t,std::shared_ptr<IChebyshevOperator> > mOps;
   };

   template <typename TOp> void Transform::addOperator(const std::size_t id)
   {
      std::size_t lid;

      if(id == 0)
      {
         lid = this->mOps.size();
      } else
      {
         lid = id;
      }

      this->mOps.insert(std::make_pair(lid, std::make_shared<TOp>()));
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_TRANSFORM_HPP
