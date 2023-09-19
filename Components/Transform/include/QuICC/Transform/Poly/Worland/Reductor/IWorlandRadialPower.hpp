/**
 * @file IWorlandRadialPower.hpp
 * @brief Interface for a Worland based power operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDRADIALPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDRADIALPOWER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandReductor.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Interface for a Worland based power operator
    */
   class IWorlandRadialPower: public IWorlandReductor
   {
      public:
         /**
          * @brief Constructor
          */
         IWorlandRadialPower();

         /**
          * @brief Destructor
          */
         virtual ~IWorlandRadialPower() = default;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const = 0;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const internal::Array& igrid, const internal::Array& iweights) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const = 0;

         /**
          * @brief Compute power on grid (squared values)
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         virtual void applyOperators(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute power (integral of squared values)
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         virtual void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;
   };

} // Reductor
} // Worland
} // Poly
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDRADIALPOWER_HPP
