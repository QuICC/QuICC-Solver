/**
 * @file ISparseSMOperator.hpp
 * @brief Implementation of the generic interface sparse operator
 */

#ifndef QUICC_SPARSESM_ISPARSESMPERATOR_HPP
#define QUICC_SPARSESM_ISPARSESMPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace SparseSM {

   /**
    * @brief Implementation of the generic interface to the sparse operator
    */
   class ISparseSMOperator
   {
      public:
         /// Typedef for scalar
         typedef internal::MHDFloat Scalar_t;

         /// Typedef for coefficient array
         typedef internal::ACoeff ACoeff_t;

         /**
          * @brief Typedef for a triplets
          *
          * (i,j,value) triplet needed to build a sparse matrix
          *
          */
         typedef Eigen::Triplet<Scalar_t>  Triplet_t;

         /// Typedef for a list of triplets
         typedef std::vector<Triplet_t>  TripletList_t;

         /**
          * @brief Constructor
          *
          * @param rows Number of rows
          * @param cols Number of columns
          */
         ISparseSMOperator(const int rows, const int cols);

         /**
          * @brief Destructor
          */
         virtual ~ISparseSMOperator() = default;

         /**
          * @brief Get sparse matrix
          */
         SparseMatrix mat() const;

         /**
          * @brief Get BLAS banded storage matrix and KL, KU
          */
         Matrix banded(unsigned int& kL, unsigned int& kU) const;

         /**
          * @brief Build matrix operator
          * @param output operator, might be banded or sparse
          * @tparam T matrix type
          *
          * Backend has no MP, call directly
          */
         template<class T, typename std::enable_if_t<std::is_same_v<Scalar_t, MHDFloat>, bool> = true>
         void buildOp(T& mat) const
         {
            this->buildOpImpl(mat);
         }

         /**
          * @brief Build sparse matrix operator
          * @param output operator
          * @tparam T matrix type
          *
          * Backend has MP, needs casting before returning the operator
          */
         template<class T, typename std::enable_if_t<!std::is_same_v<Scalar_t, MHDFloat> &&
            std::is_same_v<T, SparseMatrix>, bool> = true>
         void buildOp(T& mat) const
         {
            internal::SparseMatrix imat;
            this->buildOpImpl(imat);
            mat = imat.cast<MHDFloat>();
         }

         /**
          * @brief Build banded matrix operator
          * @param output operator
          * @tparam T matrix type
          *
          * Backend has MP, needs casting before returning the operator
          */
         template<class T, typename std::enable_if_t<!std::is_same_v<Scalar_t, MHDFloat> &&
            std::is_same_v<T, Matrix>, bool> = true>
         void buildOp(T& mat) const
         {
            internal::Matrix imat;
            this->buildOpImpl(imat);
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
          * @brief Implementation of build sparse matrix operator
          * @param output operator
          */
         void buildOpImpl(internal::SparseMatrix& mat) const;

         /**
          * @brief Implementation of build banded matrix operator
          * @param output operator
          */
         void buildOpImpl(internal::Matrix& mat) const;

         /**
          * @brief Convert diagonal to triplets
          *
          * @param list List of triplets to expand
          * @param d    Index of diagonal
          * @param rows Index of rows
          * @param diag Coefficients on diagonal
          */
         void convertToTriplets(TripletList_t& list, const int d, const ACoeffI& row, const ACoeff_t& diag) const;

      private:
         /**
          * @brief Create triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const = 0;

         /**
          * @brief Create triplet representation of matrix
          */
         virtual void buildBanded(internal::Matrix& bd, unsigned int& kL, unsigned int& kU) const;

         /**
          * @brief Handle entries to the left of the matrix (negative column index)
          *
          * Default implementation drops element
          */
         virtual void leftOutOfMatrix(TripletList_t& list, const int row, const int col, const Scalar_t value) const;

         /**
          * @brief Handle entries to the right of the matrix (column index > size of matrix)
          *
          * Default implementation drops element
          */
         virtual void rightOutOfMatrix(TripletList_t& list, const int row, const int col, const Scalar_t value) const;

         /**
          * @brief Number of rows
          */
         int mRows;

         /**
          * @brief Number of columns
          */
         int mCols;
   };

}
}

#endif // QUICC_SPARSESM_ISPARSESMPERATOR_HPP
