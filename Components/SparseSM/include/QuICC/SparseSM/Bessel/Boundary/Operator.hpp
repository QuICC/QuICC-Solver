/** 
 * @file Operator.hpp
 * @brief Implementation of the bounary operator for Bessel polynomials
 */

#ifndef QUICC_SPARSESM_BESSEL_BOUNDARY_OPERATOR_HPP
#define QUICC_SPARSESM_BESSEL_BOUNDARY_OPERATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/IBesselOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

   /**
    * @brief Implementation of the boundary operator for Bessel polynomial
    */ 
   class Operator: public IBesselOperator
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param rows    number of rows of operator
          * @param cols    number of columns of operator
          * @param l       harmonic degree l
          * @param atTop   Tau lines at top of operator
          */
         Operator(const int rows, const int cols, const BesselKind type, const int l, const bool atTop = true);

         /**
          * @brief Destructor
          */
         ~Operator() = default;

         /**
          * @brief Add boundary tau row
          */
         template <typename TCondition> void addRow();

         /**
          * @brief Embed boundary tau row in bigger size
          *
          * @param n Truncation of BC row
          */
         template <typename TCondition> void embedRow(const int n);

      private:
         /**
          * @brief Compute list of boundary values
          *
          * @param list List of triplets (row, col, value)
          */
         virtual void buildTriplets(TripletList_t& list) const;

         /**
          * @brief Harmonic degree l
          */
         int mL;

         /**
          * @brief Tau lines at top?
          */
         bool mAtTop;

         /**
          * @brief Boundary tau rows
          */
         std::vector<ACoeff_t> mBcs;
   };

   template <typename TCondition> void Operator::addRow()
   {
      TCondition bc(this->type(), this->mL);

      auto val = bc.compute(this->cols() - 1);
      this->mBcs.emplace_back(val);
   }

   template <typename TCondition> void Operator::embedRow(const int n)
   {
      TCondition bc(this->type(), this->mL);

      auto val = bc.compute(this->cols() - 1);
      if(this->cols() - n > 0)
      {
         val.bottomRows(this->cols()-n).setZero();
      }
      this->mBcs.emplace_back(val);
   }

} // Boundary
} // Bessel
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_BESSEL_BOUNDARY_OPERATOR_HPP
