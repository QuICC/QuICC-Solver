/**
 * @file ITransform.hpp
 * @brief Generic interface to a transform
 */

#ifndef QUICC_TRANSFORM_ITRANSFORM_HPP
#define QUICC_TRANSFORM_ITRANSFORM_HPP

// System includes
//
#include <set>
#include <map>

// Project includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "DenseSM/IGenericProfile.hpp"

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
         virtual void forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id);

         /**
          * @brief Compute forward transform
          *
          * @param rOut Complex output values
          * @param in   Real input values
          */
         virtual void forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id);

         /**
          * @brief Compute forward transform
          *
          * @param rOut Real output values
          * @param in   Complex input values
          */
         virtual void forward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id);

         /**
          * @brief Compute forward transform
          *
          * @param rOut Real output values
          * @param in   Real input values
          */
         virtual void forward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id);

         /**
          * @brief Compute backward transform
          *
          * @param rOut Complex output values
          * @param in   Complex input values
          */
         virtual void backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id);

         /**
          * @brief Compute backward transform
          *
          * @param rOut Complex output values
          * @param in   Real input values
          */
         virtual void backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id);

         /**
          * @brief Compute backward transform
          *
          * @param rOut Real output values
          * @param in   Complex input values
          */
         virtual void backward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id);

         /**
          * @brief Compute backward transform
          *
          * @param rOut Real output values
          * @param in   Real input values
          */
         virtual void backward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id);

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Complex output values
          * @param in   Complex input values
          */
         virtual void reduce(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id);

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Complex output values
          * @param in   Real input values
          */
         virtual void reduce(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id);

         /**
          * @brief Compute reduction transform
          * 
          * Overloaded to accept a shared pointer to a (e.g. radial) profile, pF
          *
          * @param rOut Real output values
          * @param in   Complex input values
          * @param pF   Shared pointer to the profile
          */
         virtual void reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id, std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF);

         /**
          * @brief Compute reduction transform
          * 
          * Overloaded to accept a shared pointer to a (e.g. radial) profile, pF
          *
          * @param rOut Real output values
          * @param in   Complex input values
          * @param pF   Shared pointer to the profile
          */
         virtual void reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id, std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF);

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Real output values
          * @param in   Complex input values
          */
         virtual void reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id);

         /**
          * @brief Compute reduction transform
          *
          * @param rOut Real output values
          * @param in   Real input values
          */
         virtual void reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id);

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
