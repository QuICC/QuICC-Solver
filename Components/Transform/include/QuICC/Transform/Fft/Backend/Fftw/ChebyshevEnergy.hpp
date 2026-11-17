/**
 * @file ChebyshevEnergy.hpp
 * @brief Interface for a generic Chebyshev FFTW based energy reductor
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVENERGY_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVENERGY_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/IChebyshevBackend.hpp"
#include "QuICC/Transform/Fft/Backend/StorageKind.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/DifferentialSolver.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic Chebyshev FFTW based energy reductor
    */
   // Inheritance: ChebyshevEnergy
   //                -> IChebyshevBackend 
   //                   -> IFftwBackend
   class ChebyshevEnergy: public IChebyshevBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ChebyshevEnergy() = default;

         /**
          * @brief Destructor
          */
         ~ChebyshevEnergy() = default;

         /**
          * @brief Initialise the FFTW transforms
          */
         void init(const SetupType& setup) const final;

         /**
          * @brief Set Scaler array
          */
         void setScaler(const Array& scaler) const;

         /**
          * @brief Initialise the FFTW transforms, anelastic overload
          */
         void init(const SetupType& setup, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const;

         /**
          * @brief set spectral operator
          */
         void setSpectralOperator(const SparseMatrix& mat) const;

         /**
          * @brief Square intermediate result
          */
         void square(Matrix& tmp, const Matrix& in, const bool isFirst) const;

         /**
          * @brief Set output
          */
         void output(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const;

         /**
          * @brief Set output on grid
          */
         void outputGrid(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const;

         /**
          * @brief Set output mutliplied by scalar operator
          */
         void outputSpectral(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const;

         /**
          * @brief Apply forward FFT
          */
         void applyFwdFft(Eigen::Ref<Matrix> mods, const Eigen::Ref<const Matrix>& phys) const;

         /**
          * @brief Add new linear solver
          *
          * @param extra Extra spectral modes
          */
         void addSolver(const int extra = 0) const;

         /**
          * @brief Get solver
          */
         DifferentialSolver& solver() const;

         /**
          * @brief Get solution from solver
          *
          */
         void getSolution(Matrix& tmp, const int zeroRows = 0, const int extraRows = 0) const;

         /**
          * @brief Get the temporary storage
          *
          * @param getOut return input or ouput storage
          */
         Matrix& getStorage(const StorageKind = StorageKind::in) const final;

         /**
          * @brief Apply padding
          */
         void applyPadding(Eigen::Ref<Matrix> rData, const int extraRows = 0) const final;

         /**
          * @brief Get the energy grid
          */
         Array& getEGrid() const;

         /**
          * @brief Get fft scaling
          */
         MHDFloat getFftScaling() const;
         
         void setExtraSize(int extraSize) const;

         int getExtraSize() const;

      protected:

      private:

         /**
          * @brief Compute energy weights
          */
         void computeEWeights(const int size, const MHDFloat lower, const MHDFloat upper) const;

         /**
          * @brief Compute chebyshev grid
          */
         void computeEGrid(const int size, const MHDFloat lower, const MHDFloat upper) const;

         /**
          * @brief Temporary data for mid operations
          */
         mutable Matrix  mTmpMid;

         /**
          * @brief Plan for the transform
          */
         mutable fftw_plan   mFwdPlan;

         /**
          * @brief Fwd size
          */
         mutable int mFwdSize;

         /**
          * @brief Padding size
          */
         mutable int mPadSize;

         /**
          * @brief Solver for differential operators
          */
         mutable std::shared_ptr<DifferentialSolver> mspSolver;

         /**
          * @brief Scaler array
          */
         mutable Array mScaler;

         /**
          * @brief FFT scaling factor
          */
         mutable MHDFloat mFftScaling;

         /**
          * @brief Energy weights
          */
         mutable Array mEWeights;

          /**
          * @brief Energy weights
          */
         mutable Array mEGrid;

         /**
          * @brief Spectral operator
          */
         mutable SparseMatrix mSpecOp;

         mutable int mExtraSize = 0; 

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVENERGY_HPP
