/**
 * @file PDivLlDivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) 1/Sin D_phi parallel integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_PDIVLLDIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_PDIVLLDIVS1DP_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template<typename OpTypes = PIALegendreOperatorTypes>
   class PDivLlDivS1Dp: public PDivS1<OpTypes>
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
        PDivLlDivS1Dp();

        /**
         * @brief Destructor
         */
        virtual ~PDivLlDivS1Dp();

      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const;

        /**
         * @brief l(l+1) scaling factors
         */
        virtual void initSpecial() const override;

        /**
         * @brief Storage for l(l+1) factors
         */
        mutable Array mDivLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP
