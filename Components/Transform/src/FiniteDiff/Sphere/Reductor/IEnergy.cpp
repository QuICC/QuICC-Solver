/**
 * @file IEnergy.cpp
 * @brief Source of the interface to a Finite Differences based projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/IEnergy.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

namespace Reductor {

   IEnergy::IEnergy()
   {
   }

   void IEnergy::initOperators(const Internal::Array& icompgrid) const
   {
#if 0
      // Energy calculation requires a different quadrature
      Internal::Array igrid, iweights;
      this->computePowerQuadrature(igrid, iweights, icompgrid.size());

      // Reserve storage for the operators
      this->mOps.reserve(this->mspSetup->slowSize());
      this->mEOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic degrees
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(igrid.size(), this->mspSetup->fastSize(i)));
         this->mEOps.push_back(Matrix(igrid.size(), 1));
         this->makeOperator(this->mOps.back(), this->mEOps.back(), igrid, iweights, i);
      }
#endif
   }

   void IEnergy::defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
#if 0
      rOut = (this->mEOps.at(i).leftCols(rOut.rows()).transpose()*this->mOps.at(i)*in).array().abs2();
#endif
   }

   void IEnergy::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Unused interface");
   }

   int IEnergy::outRows() const
   {
      return this->mspSetup->blockSize();
   }

   int IEnergy::outCols() const
   {
      return 1;
   }

}
}
}
}
}
