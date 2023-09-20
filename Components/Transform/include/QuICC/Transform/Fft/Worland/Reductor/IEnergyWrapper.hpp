/**
 * @file IEnergyWrapper.hpp
 * @brief Wrapper to call Poly implementation in FFT transform
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_IENERGYWRAPPER_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_IENERGYWRAPPER_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   /**
    * @brief Interface for a generic Worland FFT based energy operator
    */
   template <typename T> class IEnergyWrapper: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IEnergyWrapper();

         /**
          * @brief Destructor
          */
         virtual ~IEnergyWrapper();

         /**
          * @brief Compute transform R2R componentwise
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform R2R
          *
          * @param rOut Output values
          * @param in   Input values
          */
 	 virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

	 virtual void transform(Matrix& rOut, const Matrix& in) const override;
         virtual void transform(MatrixZ& rOut, const Matrix& in) const override;


         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief Poly setup
          */
         mutable std::shared_ptr<typename T::SetupType> mspPSetup;

         /**
          * @brief Poly operator
          */
         T mOp;

      private:
   };

   template <typename T> IEnergyWrapper<T>::IEnergyWrapper()
   {
   }

   template <typename T> IEnergyWrapper<T>::~IEnergyWrapper()
   {
   }

   template <typename T> void IEnergyWrapper<T>::initBackend() const
   {
      this->mspPSetup = std::make_shared<typename T::SetupType>(this->mspSetup->fwdSize(), this->mspSetup->specSize(),this->mspSetup->purpose());

      for(auto i = 0; i < this->mspSetup->slowSize(); i++)
      {
         this->mspPSetup->addIndex(this->mspSetup->slow(i), this->mspSetup->mult(i));
      }
      this->mspPSetup->lock();

      internal::Array igrid, iweights;
      igrid.resize(this->mspPSetup->fwdSize());
      this->mOp.init(this->mspPSetup, igrid, iweights);
   }

   template <typename T> void IEnergyWrapper<T>::transform(Matrix& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());

      this->mOp.transform(rOut, in);
   }

   template <typename T> void IEnergyWrapper<T>::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());

      this->mOp.transform(rOut, in);
   }

   template <typename T> void IEnergyWrapper<T>::transform(Matrix& rOut, const Matrix& in) const
   {
      throw std::logic_error("Unused interface");
   }

   template <typename T> void IEnergyWrapper<T>::transform(MatrixZ& rOut, const Matrix& in) const
   {
      throw std::logic_error("Unused interface");
   }

   template <typename T> int IEnergyWrapper<T>::outRows() const
   {
      return this->mOp.outRows();
   }

   template <typename T> int IEnergyWrapper<T>::outCols() const
   {
      return this->mOp.outCols();
   }

   template <typename T> MHDFloat IEnergyWrapper<T>::requiredStorage() const
   {
      return this->mOp.requiredStorage();
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_IENERGYWRAPPER_HPP
