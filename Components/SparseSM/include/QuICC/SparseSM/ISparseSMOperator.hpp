/**
 * @file ISparseSMOperator.hpp
 * @brief Implementation of the generic interface sparse operator
 */

#ifndef QUICC_SPARSESM_ISPARSESMPERATOR_HPP
#define QUICC_SPARSESM_ISPARSESMPERATOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

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
          */
         ISparseSMOperator(const int rows, const int cols);

         /**
          * @brief Destructor
          */
         virtual ~ISparseSMOperator();

         /**
          * @brief Get sparse matrix
          */
         SparseMatrix mat() const;

         /**
          * @brief Get BLAS banded storage matrix and KL, KU
          */
         Matrix banded(unsigned int& kL, unsigned int& kU) const;

         /**
          * @brief Build sparse matrix operator
          */
         void buildOp(internal::SparseMatrix& mat) const;

         /**
          * @brief Build banded matrix operator
          */
         void buildOp(internal::Matrix& mat) const;

         /**
          * @brief Build sparse matrix operator
          */
         template<class T>
            void buildOp(T& mat, typename std::enable_if_t<!std::is_same_v<Scalar_t,MHDFloat>,bool> = true) const
            {
               if constexpr(std::is_same_v<T, SparseMatrix>)
               {
                  internal::SparseMatrix imat;
                  this->buildOp(imat);
                  mat = imat.cast<MHDFloat>();
               }
               else if constexpr(std::is_same_v<T, Matrix>)
               {
                  internal::Matrix imat;
                  this->buildOp(imat);
                  mat = imat.cast<MHDFloat>();
               }
               else
               {
                  static_assert(true, "Unknown type");
               }
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
