/**
 * @file PD1.hpp
 * @brief Implementation of the associated Legendre based D projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_PD1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_PD1_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Projector/PIALegendreProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   /**
    * @brief Implementation of the associated Legendre based D projector
    */
   template<typename OpTypes = PIALegendreOperatorTypes>
   class PD1: public PIALegendreProjector<OpTypes>
   {
      public:
        using OpArray = typename OpTypes::OpArray;
        using OpMatrix = typename OpTypes::OpMatrix;
        using OpMatrixR = typename OpTypes::OpMatrixR;
        using OpMatrixCR = typename OpTypes::OpMatrixCR;

        using OpVectorI = typename OpTypes::OpVectorI;
        using OpMatrixI = typename OpTypes::OpMatrixI;
        using OpMatrixLZ = typename OpTypes::OpMatrixLZ;
        using OpMatrixL = typename OpTypes::OpMatrixL;

        using DataType = typename OpTypes::DataType;

        /**
         * @brief Constructor
         */
        PD1();

        /**
         * @brief Destructor
         */
        virtual ~PD1();

        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const;

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D1_HPP
