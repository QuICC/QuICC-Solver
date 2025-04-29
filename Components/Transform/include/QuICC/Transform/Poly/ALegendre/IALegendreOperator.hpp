/**
 * @file IALegendreOperator.hpp
 * @brief Interface for a associated Legendre based operator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_IALEGENDREOPERATOR_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_IALEGENDREOPERATOR_HPP

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
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"
#include "QuICC/Transform/ITransformOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

   /**
    * @brief Interface for a associated Legendre based operator
    */
   class IALegendreOperator: public ITransformOperator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IALegendreOperator();

         /**
          * @brief Destructor
          */
         virtual ~IALegendreOperator() = default;

         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         virtual void init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const override;

         void init(SharedTransformSetup spSetup) const override;

         /**
          * @brief Compute polynomial transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute polynomial transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const MatrixZ& in) const;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const = 0;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const = 0;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

      protected:
         /**
          * @brief Polynomial setup object providing the sizes
          */
         mutable SharedSetup    mspSetup;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const = 0;

         /**
          * @brief Apply operators
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyOperators(MatrixZ& rOut, const MatrixZ& in) const;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_IALEGENDREOPERATOR_HPP
