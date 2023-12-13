/**
 * @file EnergyReductor.hpp
 * @brief Interface for a Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYREDUCTOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYREDUCTOR_HPP

// System includes
//

// Project includes
//
#include "Profiler/Interface.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandReductor.hpp"
#include "Types/Typedefs.hpp"


namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

/**
 * @brief Interface for a Worland based energy operator
 */
template <typename T> class KokkosEnergyReductor : public T
{
public:
   using OpVectorI = typename KokkosIOperatorTypes::OpVectorI;
   using OpMatrixLZ = typename KokkosIOperatorTypes::OpMatrixLZ;
   using OpMatrixLZL = typename KokkosIOperatorTypes::OpMatrixLZL;
   using OpMatrixL = typename KokkosIOperatorTypes::OpMatrixL;

   /**
    * @brief Constructor
    */
   KokkosEnergyReductor() = default;

   /**
    * @brief Destructor
    */
   virtual ~KokkosEnergyReductor() = default;

   /**
    * @brief Rows of output data
    */
   virtual int outRows() const override;

   /**
    * @brief Columns of output data
    */
   virtual int outCols() const override;

protected:
   virtual void defaultApplyUnitOperator(const OpMatrixL& rOut,
      const OpMatrixLZL& in, const OpVectorI& scan,
      const int totalOpsCols) const;

private:
   /**
    * @brief Compute energy (integral of squared values)
    *
    * @param rOut Output physical values
    * @param in   Input spectral coefficients
    */
   virtual void applyOperators(Matrix& rOut, const MatrixZ& in) const override;

   /**
    * @brief Compute energy (integral of squared values)
    *
    * @param rOut Output physical values
    * @param in   Input spectral coefficients
    */
   virtual void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;
};

template <typename T>
void KokkosEnergyReductor<T>::applyOperators(Matrix& rOut,
   const MatrixZ& in) const
{
   Profiler::RegionFixture<3> fix("KokkosEnergyReductor::applyOperators");

   // assert right sizes for input  matrix
   assert(in.cols() == this->mspSetup->blockSize());
   // assert right sizes for output matrix
   assert(rOut.rows() == this->outRows());
   assert(rOut.cols() == this->outCols());

   auto slowSize = this->mspSetup->slowSize();
   OpVectorI scan("outRows Scan", slowSize + 1);
   auto hostScan = Kokkos::create_mirror_view(scan);

   Profiler::RegionStart<4>("KokkosEnergyReductor::hostScan");

   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      int cols = this->mspSetup->mult(i);
      hostScan[i + 1] = hostScan[i] + cols;
   }

   Profiler::RegionStop<4>("KokkosEnergyReductor::hostScan");

   Kokkos::deep_copy(scan, hostScan);

   auto total = hostScan(slowSize);
   OpMatrixLZL inView("inView", in.rows(), in.cols());
   DeepCopyEigen(inView, in);
   OpMatrixL rOutView("rOutView", rOut.rows(), rOut.cols());

   Profiler::RegionStart<4>("KokkosEnergyReductor::applyUnitOperator");
   this->applyUnitOperator(rOutView, inView, scan, total);
   Profiler::RegionStop<4>("KokkosEnergyReductor::applyUnitOperator");

   DeepCopyEigen(rOut, rOutView);
}

template <typename T>
void KokkosEnergyReductor<T>::applyOperators(MatrixZ& rOut,
   const MatrixZ& in) const
{
   throw std::logic_error("Unused interface");
}

template <typename T>
void KokkosEnergyReductor<T>::defaultApplyUnitOperator(
   const OpMatrixL& rOutView, const OpMatrixLZ& inView, const OpVectorI& scan,
   const int total) const
{
   OpMatrixLZ temp("Temp View", this->vmOps.extent(0), inView.extent(1));
   applyBlockOperator<4>(this->mspSetup, this->vmOps, temp, inView, scan,
      total);
   OpMatrixL tempRout("rOut Temp View", inView.extent(0), inView.extent(1));
   applyBlockOperator<3>(this->mspSetup, this->vmEOps, tempRout, temp, scan,
      total, Abs2Complex());
   colwiseSum(tempRout, rOutView);
}

template <typename T> int KokkosEnergyReductor<T>::outRows() const
{
   return this->mspSetup->blockSize();
}

template <typename T> int KokkosEnergyReductor<T>::outCols() const
{
   return 1;
}
} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYREDUCTOR_HPP
