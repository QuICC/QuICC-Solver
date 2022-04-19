/**
 * @file ITransform.hpp
 * @brief Generic interface to a transform
 */

#ifndef QUICC_TRANSFORM_ITRANSFORM_HPP
#define QUICC_TRANSFORM_ITRANSFORM_HPP

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

namespace QuICC {

namespace Transform {

   /**
    * @brief Generic interface to a transform
    */
   class ITransform
   {
      public:
         /**
          * @brief Constructor
          */
         ITransform() = default;

         /**
          * @brief Destructor
          */
         virtual ~ITransform() = default;

         /**
          * @brief set list of required options
          */
         virtual void requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const = 0;

         /**
          * @brief Set the required options
          */
         virtual void setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId) = 0;

         /**
          * @brief Get the physical grid
          */
         virtual Array meshGrid() const = 0;

         /**
          * @brief Compute forward transform
          *
          * @param rOut Complex output values
          * @param in   Complex input values
          */
         virtual void forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id) = 0;

         /**
          * @brief Compute forward transform
          *
          * @param rOut Complex output values
          * @param in   Real input values
          */
         virtual void forward(MatrixZ& rOut, const Matrix& in, const std::size_t id) = 0;

         /**
          * @brief Compute forward transform
          *
          * @param rOut Real output values
          * @param in   Complex input values
          */
         virtual void forward(Matrix& rOut, const MatrixZ& in, const std::size_t id) = 0;

         /**
          * @brief Compute forward transform
          *
          * @param rOut Real output values
          * @param in   Real input values
          */
         virtual void forward(Matrix& rOut, const Matrix& in, const std::size_t id) = 0;

         /**
          * @brief Compute backward transform
          *
          * @param rOut Complex output values
          * @param in   Complex input values
          */
         virtual void backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id) = 0;

         /**
          * @brief Compute backward transform
          *
          * @param rOut Complex output values
          * @param in   Real input values
          */
         virtual void backward(MatrixZ& rOut, const Matrix& in, const std::size_t id) = 0;

         /**
          * @brief Compute backward transform
          *
          * @param rOut Real output values
          * @param in   Complex input values
          */
         virtual void backward(Matrix& rOut, const MatrixZ& in, const std::size_t id) = 0;

         /**
          * @brief Compute backward transform
          *
          * @param rOut Real output values
          * @param in   Real input values
          */
         virtual void backward(Matrix& rOut, const Matrix& in, const std::size_t id) = 0;

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Complex output values
          * @param in   Complex input values
          */
         virtual void reduce(MatrixZ& rOut, const MatrixZ& in, const std::size_t id) = 0;

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Complex output values
          * @param in   Real input values
          */
         virtual void reduce(MatrixZ& rOut, const Matrix& in, const std::size_t id) = 0;

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Real output values
          * @param in   Complex input values
          */
         virtual void reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id) = 0;

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Real output values
          * @param in   Real input values
          */
         virtual void reduce(Matrix& rOut, const Matrix& in, const std::size_t id) = 0;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const = 0;

         /**
          * @brief Profile the memory requirements
          */
         virtual void profileStorage() const = 0;

      protected:
         /**
          * @brief Generic function for unimplemented transform
          */
         void unimplemented();

      private:
   };

}
}

#endif // QUICC_TRANSFORM_ITRANSFORM_HPP
