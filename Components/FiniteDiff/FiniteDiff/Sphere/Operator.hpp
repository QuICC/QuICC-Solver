/**
 * @file Operator.hpp
 * @brief Implementation of generic operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_OPERATOR_HPP
#define QUICC_FINITEDIFF_SPHERE_OPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   /**
    * @brief Implementation of generic operator
    */
   class Operator
   {
      public:
         /**
          * @brief Constructor
          */
         Operator(const size_t order, const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          */
         Operator(const size_t order);

         /**
          * @brief Destructor
          */
         ~Operator() = default;

      protected:
         /**
          * @brief Convert order to stencil size
          *
          * @param m derivative order
          * @param n accuracy order
          * @param isCentral  is central differences?
          */
         std::size_t order2Stencil(const std::size_t m, const std::size_t n, const bool isCentral) const;

         /**
          * @brief Compute finite differences weights following Fornberg, 1998
          */
         void fdWeights(std::vector<std::vector<Internal::MHDFloat> >& w, const Internal::MHDFloat z, const std::vector<Internal::MHDFloat>& x, const std::size_t m) const;

         /**
          * @brief Compute finite differences matrices
          */
         void fdMatrices(std::vector<Internal::SparseMatrix>& w, const Internal::Array& grid, const std::size_t s, const std::size_t m) const;

         /**
          * @brief Get operator to zero top and bottom rows
          */
         Internal::SparseMatrix zeroTopBottom(const int nR, const int zTop, const int zBot) const;

         /**
          * @brief Get operator to zero top and bottom rows
          */
         Internal::SparseMatrix zeroTopBottom(const int nR) const;

         /**
          * @brief Scheme order
          */
         std::size_t mOrder;

         /**
          * @brief Zero rows at top
          */
         std::size_t mZtop;

         /**
          * @brief Zero rows at bottom
          */
         std::size_t mZbot;
   };

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_OPERATOR_HPP
