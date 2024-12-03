/**
 * @file IMatrixSMOperator.hpp
 * @brief Implementation of the generic interface for a dense spectral
 * matrix operator
 */

#ifndef QUICC_DENSESM_IMATRIXSMOPERATOR_HPP
#define QUICC_DENSESM_IMATRIXSMOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

/**
 * @brief Implementation of the generic interface to the a dense
 * spectral matrix operator
 */
class IMatrixSMOperator
{
public:
   /// Typedef for scalar
   typedef Internal::MHDFloat Scalar_t;

   /// Typedef for coefficient array
   typedef Internal::ACoeff ACoeff_t;

   /**
    * @brief Constructor
    *
    * @param rows Number of rows
    * @param cols Number of columns
    */
   IMatrixSMOperator(const int rows, const int cols);

   /**
    * @brief Destructor
    */
   virtual ~IMatrixSMOperator() = default;

   /**
    * @brief Get dense matrix
    */
   Matrix mat() const;

   /**
    * @brief Build matrix operator
    * @param output operator, might be banded or dense
    * @tparam T matrix type
    *
    * Backend has no MP, call directly
    */
   template <class T, typename TScalar = Scalar_t,
      typename std::enable_if_t<std::is_same_v<TScalar, MHDFloat>, bool> = true>
   void buildOp(T& mat) const
   {
      this->buildOpImpl(mat, this->rows(), this->cols());
   }

   /**
    * @brief Build dense matrix operator
    * @param output operator
    * @tparam T matrix type
    *
    * Backend has MP, needs casting before returning the operator
    */
   template <class T, typename TScalar = Scalar_t,
      typename std::enable_if_t<!std::is_same_v<TScalar, MHDFloat> &&
                                   std::is_same_v<T, Matrix>,
         bool> = true>
   void buildOp(T& mat) const
   {
      Internal::Matrix imat;
      this->buildOpImpl(imat, this->rows(), this->cols());
      mat = imat.cast<MHDFloat>();
   }

   /**
    * @brief Number of rows
    */
   int rows() const;

   /**
    * @brief Number of columns
    */
   int cols() const;

protected:
   /**
    * @brief Implementation of build dense matrix operator
    * @param output operator
    */
   virtual void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const = 0;

private:
   /**
    * @brief Number of rows
    */
   int mRows;

   /**
    * @brief Number of columns
    */
   int mCols;
};

} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_IMATRIXSMOPERATOR_HPP
