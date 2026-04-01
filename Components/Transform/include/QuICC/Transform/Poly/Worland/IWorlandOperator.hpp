/**
 * @file IWorlandOperator.hpp
 * @brief Interface for a Worland based operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_IWORLANDOPERATOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_IWORLANDOPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
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
         virtual ~IWorlandOperator() = default;

         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const override;
         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup) const override;

         /**
          * @brief Compute polynomial transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief Compute polynomial transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in) const;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

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
         virtual void initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const = 0;

         /**
          * @brief Apply operators
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyOperators(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in) const;

         /**
          * @brief Apply operators
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyOperators(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in) const;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_IWORLANDOPERATOR_HPP
