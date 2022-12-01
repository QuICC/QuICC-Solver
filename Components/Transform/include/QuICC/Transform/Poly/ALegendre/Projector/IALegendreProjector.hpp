/** 
 * @file IALegendreProjector.hpp
 * @brief Interface for a associated Legendre based projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_IALEGENDREPROJECTOR_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_IALEGENDREPROJECTOR_HPP

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
#include "QuICC/Transform/Poly/ALegendre/IALegendreOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   /**
    * @brief Interface for a associated Legendre based projector
    */ 
    template<typename OpTypes = IALegendreOperatorTypes>
   class IALegendreProjector: public IALegendreOperator<OpTypes>
   {
      public:
         using OpArray = typename OpTypes::OpArray;
         using OpMatrix = typename OpTypes::OpMatrix;
         using OpMatrixZ = typename OpTypes::OpMatrixZ;
         using OpMatrixR = typename OpTypes::OpMatrixR;
         using OpMatrixCR = typename OpTypes::OpMatrixCR;

         /**
          * @brief Constructor
          */
         IALegendreProjector();

         /**
          * @brief Destructor
          */
         virtual ~IALegendreProjector();

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
          * @brief Storage for the operators 
          */
         mutable std::vector<OpMatrix>  mOps;

        /**
         * @brief Storage for the quadrature grid
         */
        mutable OpArray mGrid;

        /**
         * @brief Storage for the quadrature weights
         */
        mutable OpArray mWeights;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const OpArray& igrid, const OpArray& iweights) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const = 0;

         /**
          * @brief Compute projection
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyOperators(OpMatrixZ& rOut, const OpMatrixZ& in) const override;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const = 0;

         /**
          * @brief Allow specific initialisation for operators
          */
         virtual void initSpecial() const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_IALEGENDREPROJECTOR_HPP
