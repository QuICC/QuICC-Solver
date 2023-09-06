/** 
 * @file Operator.hpp
 * @brief Implementation of the bounary operator for Chebyshev linear map polynomials
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_OPERATOR_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_OPERATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   /**
    * @brief Implementation of the boundary operator for Chebyshev linear map polynomial
    */ 
   class Operator: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param rows    number of rows of operator
          * @param cols    number of columns of operator
          * @param lower   Lower bound of y
          * @param upper   Upper bound of y
          * @param atTop   Tau lines at top of operator
          */
         Operator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const bool atTop = true);

         /**
          * @brief Destructor
          */
         ~Operator() = default;

         /**
          * @brief Add boundary tau row
          *
          * @param pos  Position of boundary condition
          */
         template <typename TCondition> void addRow(const ICondition::Position pos);

         /**
          * @brief Add boundary tau row which depends on harmonic degree
          *
          * @param pos  Boundary position
          * @param l    harmonic degree
          */
         template <typename TCondition> void addRow(const ICondition::Position pos, const int l);

      private:
         /**
          * @brief Compute list of boundary values
          *
          * @param list List of triplets (row, col, value)
          */
         virtual void buildTriplets(TripletList_t& list) const;

         /**
          * @brief Tau lines at top?
          */
         bool mAtTop;

         /**
          * @brief Lower bound
          */
         Scalar_t mLower;

         /**
          * @brief Upper bound
          */
         Scalar_t mUpper;

         /**
          * @brief Boundary tau rows
          */
         std::vector<ACoeff_t> mBcs;
   };

   template <typename TCondition> void Operator::addRow(const ICondition::Position pos)
   {
      TCondition bc(this->mLower, this->mUpper, pos);

      auto val = bc.compute(this->cols() - 1);
      this->mBcs.emplace_back(val);
   }

   template <typename TCondition> void Operator::addRow(const ICondition::Position pos, const int l)
   {
      TCondition bc(this->mLower, this->mUpper, pos, l);

      auto val = bc.compute(this->cols() - 1);
      this->mBcs.emplace_back(val);
   }

} // Boundary
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_OPERATOR_HPP
