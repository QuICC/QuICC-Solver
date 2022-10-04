/** 
 * @file D.hpp
 * @brief Implementation of the Chebyshev based D projector, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_D_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_D_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/IChebyshevProjector.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   /**
    * @brief Implementation of the Chebyshev based D projector, with linear map y = ax + b
    *
    * @tparam DO Derivative order
    */
   template <int DO>
   class D: public IChebyshevProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D(){};

         /**
          * @brief Destructor
          */
         ~D(){};

      protected:
         /**
          * @brief Initialize storage
          */
         void initBackend() const final
         {
            // Call parent initializer
            IChebyshevProjector::initBackend();

            // Initialize the solver
            this->mBackend.addSolver();
         }

      private:
         // Nasty typedef to map the implementation of
         // the I sparse integration operator
         typedef std::conditional_t<DO == 1 ,
            ::QuICC::SparseSM::Chebyshev::LinearMap::I1,
            std::conditional_t<DO == 2 ,
               ::QuICC::SparseSM::Chebyshev::LinearMap::I2,
               std::conditional_t<DO == 3 ,
                  ::QuICC::SparseSM::Chebyshev::LinearMap::I3,
                  std::conditional_t<DO == 4 ,
                     ::QuICC::SparseSM::Chebyshev::LinearMap::I4,
                     void>
                  >
               >
            > III;
         /**
          * @brief Initialize solver operators
          */
         void initOperator() const final
         {
            static_assert(!std::is_same_v<III, void>, "missing sparse integrator.");
            III
               op(this->mspSetup->specSize()+DO,
                  this->mspSetup->specSize()+DO,
                  this->mspSetup->lower(),
                  this->mspSetup->upper());

            this->mBackend.solver().setOperator(op.mat());
         }

         /**
          * @brief Apply pre FFT operator
          *
          * @param tmp  Temporary padded modal values
          * @param in   Input values
          */
         void applyPreOperator(Matrix& tmp, const Matrix& in) const final
         {
            this->mBackend.input(tmp, in, DO);
            this->mBackend.getSolution(tmp, DO);
         }

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(Matrix& rOut) const final{};

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param tmp Temporary padded modal values, either real or imag part only
          * @param in Input values
          * @param useReal Real vs Imag flag
          */
         void applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const final
         {
            this->mBackend.input(tmp, in, DO, useReal);
            this->mBackend.getSolution(tmp, DO);
         }

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         void applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const final
         {
            this->mBackend.output(rOut, tmp, useReal);
         }
   };

   /**
    * @brief type alias for backward compatibility
    */
   using D1 = D<1>;
   /**
    * @brief type alias for backward compatibility
    */
   using D2 = D<2>;
   /**
    * @brief type alias for backward compatibility
    */
   using D3 = D<3>;
   /**
    * @brief type alias for backward compatibility
    */
   using D4 = D<4>;

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_D1_HPP
