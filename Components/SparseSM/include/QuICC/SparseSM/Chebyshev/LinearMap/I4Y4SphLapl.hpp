/** 
 * @file I4Y4SphLapl.hpp
 * @brief Implementation of the I^4 Y^4 spherical laplacian sparse operator, with y = ax + b
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4Y4SPHLAPL_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4Y4SPHLAPL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/ISphericalOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the I^4 Y^4 spherical laplacian sparse operator, with y = ax + b
    */ 
   class I4Y4SphLapl: public ISphericalOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of cols
          * @param lower   Lower bound
          * @param upper   Upper bound
          * @param l       Harmonic degree l
          */
         I4Y4SphLapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l);

         /**
          * @brief Destructor
          */
         virtual ~I4Y4SphLapl() = default;
         
      protected:

      private:
         /**
          * @brief 6th subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_6(const ACoeff_t& n) const; 

         /**
          * @brief 5th subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_5(const ACoeff_t& n) const; 

         /**
          * @brief 4th subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_4(const ACoeff_t& n) const; 

         /**
          * @brief 3rd subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_3(const ACoeff_t& n) const; 

         /**
          * @brief 2nd subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_2(const ACoeff_t& n) const; 

         /**
          * @brief 1st subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_1(const ACoeff_t& n) const; 

         /**
          * @brief diagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const; 

         /**
          * @brief 1st superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d1(const ACoeff_t& n) const; 

         /**
          * @brief 2nd superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d2(const ACoeff_t& n) const; 

         /**
          * @brief 3rd superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d3(const ACoeff_t& n) const; 

         /**
          * @brief 4th superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d4(const ACoeff_t& n) const; 

         /**
          * @brief 5th superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d5(const ACoeff_t& n) const; 

         /**
          * @brief 6th superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d6(const ACoeff_t& n) const; 

         /**
          * @brief Build triplet representation of matrix
          *
          * @param[out] list containing triplets
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4Y4SPHLAPL_HPP
