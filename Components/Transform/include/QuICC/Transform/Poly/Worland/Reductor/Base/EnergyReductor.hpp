/**
 * @file EnergyReductor.hpp
 * @brief Interface for a Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYREDUCTOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYREDUCTOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandReductor.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Interface for a Worland based energy operator
    */
   template <typename T> class EnergyReductor: public T
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyReductor() = default;

         /**
          * @brief Destructor
          */
         virtual ~EnergyReductor() = default;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

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

   template <typename T> void EnergyReductor<T>::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag);

      // assert right sizes for input  matrix
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.rows() == this->outRows());
      assert(rOut.cols() == this->outCols());

      int start = 0;
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int inRows = this->mspSetup->fastSize(i);

         this->applyOperator(rOut.block(start,0, cols, 1), i, in.block(0,start, inRows, cols));

         start += cols;
      }
   }

   template <typename T> void EnergyReductor<T>::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Unused interface");	   
   }


   template <typename T> void EnergyReductor<T>::defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut.transpose() = (this->mEOps.at(i).transpose()*this->mOps.at(i)*in).array().abs2().colwise().sum();
   }

   template <typename T> int EnergyReductor<T>::outRows() const
   {
      return this->mspSetup->blockSize();
   }

   template <typename T> int EnergyReductor<T>::outCols() const
   {
      return 1;
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYREDUCTOR_HPP
