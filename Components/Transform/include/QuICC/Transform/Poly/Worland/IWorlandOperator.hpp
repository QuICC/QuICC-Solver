/**
 * @file IWorlandOperator.hpp
 * @brief Interface for a Worland based operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_IWORLANDOPERATOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_IWORLANDOPERATOR_HPP

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
#include "QuICC/Precision.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"
#include "QuICC/Transform/ITransformOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

   /**
    * @brief Interface for a Worland based operator
    */
   class IWorlandOperator: public ITransformOperator
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
         virtual ~IWorlandOperator();

         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup, const internal::Array& igrid, const internal::Array& iweights) const;

         /**
          * @brief Compute polynomial transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(MatrixZ& rOut, const MatrixZ& in) const;

         /**
          * @brief Compute polynomial transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const MatrixZ& in) const;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const = 0;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const = 0;

      protected:
         /**
          * @brief Check grid size can accommodate exact computation
          */
         void checkGridSize(const int n, const int l, const int nG) const;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         mutable SharedSetup    mspSetup;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const internal::Array& igrid, const internal::Array& iweights) const = 0;

         /**
          * @brief Apply operators
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyOperators(MatrixZ& rOut, const MatrixZ& in) const;

         /**
          * @brief Apply operators
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyOperators(Matrix& rOut, const MatrixZ& in) const;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_IWORLANDOPERATOR_HPP
