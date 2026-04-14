/**
 * @file IWorlandOperator.hpp
 * @brief Interface for a generic Worland FFT based operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_IWORLANDOPERATOR_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_IWORLANDOPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"
#include "QuICC/Transform/Fft/IFftOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

   /**
    * @brief Interface for a generic Worland FFT based operator
    */
   class IWorlandOperator: public IFftOperator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IWorlandOperator();

         /**
          * @brief Destructor
          */
         virtual ~IWorlandOperator() = default;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup) const override;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in) const override;

      protected:
         /**
          * @brief Polynomial setup object providing the sizes
          */
         mutable SharedSetup    mspSetup;

         /**
          * @brief Compute transform R2C (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in) const override;
   };

} // namespace Worland
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_WORLAND_IWORLANDOPERATOR_HPP
