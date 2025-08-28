/**
 * @file IEnergy.hpp
 * @brief Interface for a Finite Difference based energy operator
 */

#ifndef QUICC_TRANSFORM_FINITEDIFF_SPHERE_REDUCTOR_IENERGY_HPP
#define QUICC_TRANSFORM_FINITEDIFF_SPHERE_REDUCTOR_IENERGY_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/IReductor.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

namespace Reductor {

   /**
    * @brief Interface for a Finite Differences based energy operator
    */
   class IEnergy: public IReductor
   {
      public:
         /**
          * @brief Constructor
          */
         IEnergy();

         /**
          * @brief Destructor
          */
         virtual ~IEnergy() = default;

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
         virtual void initOperators(const Internal::Array& igrid) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const int i) const = 0;

         /**
          * @brief Compute power (integral of squared values)
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         virtual void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;
   };

} // namespace Reductor
} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FINITEDIFF_SPHERE_REDUCTOR_IENERGY_HPP
