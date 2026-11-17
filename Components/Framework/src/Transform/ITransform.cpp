/**
 * @file ITransform.cpp
 * @brief Source of the generic transform interface
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/ITransform.hpp"

namespace QuICC {

namespace Transform {

   void ITransform::unimplemented()
   {
      throw std::logic_error("Transform does not implemnent this option");
   }

   void ITransform::forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->unimplemented();
   }

   void ITransform::forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id)
   {
      this->unimplemented();
   }

   void ITransform::forward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->unimplemented();
   }

   void ITransform::forward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id)
   {
      this->unimplemented();
   }

   void ITransform::backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->unimplemented();
   }

    void ITransform::backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id)
    {
      this->unimplemented();
    }

    void ITransform::backward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
    {
      this->unimplemented();
    }

    void ITransform::backward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id)
    {
      this->unimplemented();
    }

    void ITransform::reduce(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
    {
      this->unimplemented();
    }

    void ITransform::reduce(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id)
    {
      this->unimplemented();
    }

    void ITransform::reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
    {
      this->unimplemented();
    }

    void ITransform::reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id)
    {
      this->unimplemented();
    }

    void ITransform::reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id, std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF)
    {
       this->reduce(rOut, in, id);
    }

    void ITransform::reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id, std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF)
    {
       this->reduce(rOut, in, id);
    }
}
}
