/**
 * @file Transform.hpp
 * @brief Implementation of the mixed Fourier transform
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_TRANSFORM_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_TRANSFORM_HPP

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
#include "QuICC/Transform/Fft/Fourier/Mixed/IMixedOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

   /**
    * @brief Implementation of the mixed Fourier transform
    */
   class Transform
   {
      public:
         /// Typedef for the configuration class
         typedef IMixedOperator::SetupType SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef IMixedOperator::SharedSetupType SharedSetupType;

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
         void transform(MatrixZ& rOut, const Matrix& in, const IMixedOperator& op);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const MatrixZ& in, const IMixedOperator& op);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(MatrixZ& rOut, const Matrix& in, const std::size_t id);

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const MatrixZ& in, const std::size_t id);

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size);

         /**
          * @brief Polynomial setup object providing the sizes
          */
         SharedSetupType    mspSetup;

         /**
          * @brief Store transform operator to ID mapping
          */
         std::map<std::size_t,std::shared_ptr<IMixedOperator> > mOps;
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

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_TRANSFORM_HPP
