/**
 * @file ChebyshevRule.hpp
 * @brief Implementation of the Chebyshev quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_CHEBYSHEVRULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_CHEBYSHEVRULE_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "QuICC/Polynomial/Quadrature/traits.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

/**
 * @brief Implementation of the Chebyshev quadrature rule
 */
template <typename quadType = gauss_t>
class ChebyshevRule
{
   public:
      /**
       * @brief Compute the ith zero for Gauss quadrature
       */
      template <class TAG,
         typename std::enable_if<std::is_same<gauss_t, TAG>::value, bool>::type = 0
      >
      static inline Internal::MHDFloat computeZero(const std::uint32_t k, const std::uint32_t size)
      {
         assert(k >= 0);
         assert(size > 1);
         return -Internal::Math::cos(static_cast<Internal::MHDLong>(Internal::Math::PI_long
            *(Internal::MHDFloat(k)+ MHD_MP(0.5))/Internal::MHDFloat(size)));
      }

      /**
       * @brief Compute the ith zero for Gauss-Lobatto quadrature
       */
      template <class TAG,
         typename std::enable_if<std::is_same<gaussLobatto_t, TAG>::value, bool>::type = 0
      >
      static inline Internal::MHDFloat computeZero(const std::uint32_t k, const std::uint32_t size)
      {
         assert(k >= 0);
         assert(size >= 2);
         // use sin to enforce simmetry
         Internal::MHDFloat Nm1 = Internal::MHDFloat(size-1);
         Internal::MHDFloat twoNm1 = Internal::MHDFloat(2.0)*Nm1;
         Internal::MHDFloat pos = Internal::MHDFloat(2.0)*Internal::MHDFloat(k) - Nm1;
         return Internal::Math::sin(static_cast<Internal::MHDLong>(Internal::Math::PI_long
            *pos/twoNm1));
      }

      /**
       * @brief Compute the quadrature
       */
      void computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const std::uint32_t size)
      {
         // Internal grid and weights arrays
         igrid.resize(size);
         iweights.resize(size);

         for(std::uint32_t k = 0; k < size; k++)
         {
            igrid(k) = -computeZero<quadTypeTag_t>(k, size);
         }

         // Set exact zero
         if(size % 2 == 1)
         {
            igrid(size/2) = 0.0;
         }

         iweights.setConstant(Internal::Math::PI/Internal::MHDFloat(size));

         if constexpr (std::is_same_v<quadTypeTag_t, gaussLobatto_t>)
         {
            iweights(0) *= 0.5;
            iweights(size-1) *= 0.5;
         }
      }

      /**
       * @brief Compute the quadrature with lower and upper bound
       */
      void computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const std::uint32_t size, const MHDFloat lower, const MHDFloat upper)
      {
         computeQuadrature(igrid, iweights, size);

         MHDFloat a = (upper - lower)/2.0;
         MHDFloat b =  (upper + lower)/2.0;
         igrid.array() = a*igrid.array() + b;
      }


      /**
       * @brief Compute the collocation differentiation matrix (phys -> phys)
       *
       * note: need more tricks for high n
       *
       * ref
       *    Spectral differencing with a twist
       *    Baltensperger and Trummer 2003
       *
       */
      void computeDiffMat(Internal::Matrix& D)
      {
         if constexpr (std::is_same_v<quadTypeTag_t, gaussLobatto_t>)
         {
            auto nPts = D.rows();
            auto nPtsM1 = nPts - 1;
            // zero everything
            D = Internal::Matrix::Zero(nPts,nPts);

            // this could be stored
            std::vector<Internal::MHDFloat> grid(nPts);
            for(int k = 0; k < nPts; k++)
            {
               grid[k] = computeZero<quadTypeTag_t>(k, nPts);
            }

            // Set exact zero
            if(nPts % 2 == 1)
            {
               grid[nPts/2] = 0.0;
            }

            auto fi = 0.0;
            auto fj = 0.0;
            std::int32_t ci{}, cj{};
            Internal::Matrix::Index i{};
            for (i = 0, ci = 1; i < nPts; ++i, ci = -ci)
            {
               if (i == 0 || i == nPts-1)
               {
                  fi = 2.0;
               }
               else
               {
                  fi = 1.0;
               }
               Internal::Matrix::Index j{};
               for (j = 0, cj = 1; j < nPts; ++j, cj = -cj)
               {
                  if (j == 0 || j == nPts-1)
                  {
                     fj = 0.5;
                  }
                  else
                  {
                     fj = 1.0;
                  }
                  D(i,j) = fi*fj*cj*ci/(grid[i]-grid[j]);
               }
            }
            for (Internal::Matrix::Index i = 1; i < nPtsM1; ++i)
            {
               D(i,i) = 0.0;
               Internal::MHDFloat tmp= 0.0;
               for (Internal::Matrix::Index j = 0; j < nPts; ++j)
               {
                  tmp -= D(i,j);
               }
               D(i,i) = tmp;
            }
            D(0,0) = -(2*nPtsM1*nPtsM1+1)/6.0;
            D(nPtsM1, nPtsM1) = -D(0,0);
         }
         else
         {
            throw std::logic_error("Chebyshev Gauss diff mat not implemented");
         }
      }

      /**
       * @brief Compute the integration (operational) matrix (modal -> modal)
       *
       * ref
       *    Chebyshev series approach to system identification, analysis and optimal control
       *    N.Paraskevopoulos 1983
       *
       */
      void computeIntMat(Internal::Matrix& I)
      {
         if constexpr (std::is_same_v<quadTypeTag_t, gaussLobatto_t>)
         {
            auto nPts = I.rows();
            // zero everything
            I = Internal::Matrix::Zero(nPts,nPts);

            // alpha_0
            I(0,0) = 1.0;
            // beta_0
            I(0,1) = 1.0;
            // alpha_1
            I(1,0) = -0.25;
            // gamma_1
            I(1,2) = 0.25;
            auto sgn = -1.0;
            for (Internal::Matrix::Index i = 2; i < nPts; ++i)
            {
               // alpha_n
               I(i,0) = sgn/(i*i-1.0);
               sgn = -sgn;
               // gamma_n
               if (i < nPts -1)
               {
                  I(i,i+1) = 1.0/(2.0*(i+1.0));
               }
               // beta_n
               I(i,i-1) = -1.0/(2.0*(i-1.0)) ;
            }
            I.transposeInPlace();
         }
         else
         {
            throw std::logic_error("Chebyshev Gauss int mat not implemented");
         }
      }

      /**
       * @brief Forward/Backward Vandermonde matrix
       *  FWD = 1       -> (modal -> phys)
       *  BWD (FWD = 0) -> (phys -> modal)
       *
       * note: need more tricks for high n
       *
       * ref
       *    Spectral differencing with a twist
       *    Baltensperger and Trummer 2003
       *
       */
      template <bool FWD>
      void computeTMat(Internal::Matrix& T)
      {
         if constexpr (std::is_same_v<quadTypeTag_t, gaussLobatto_t>)
         {
            auto nPts = T.rows();
            Internal::Array grid(nPts), weights(nPts);
            Quadrature::ChebyshevRule<Quadrature::gaussLobatto_t> quad;
            quad.computeQuadrature(grid, weights, nPts, 1.0, -1.0);

            Internal::MHDFloat one_n = 1.0 / (nPts-1);
            Internal::MHDFloat c = 1.0;
            for (Internal::Matrix::Index i = 0; i < nPts; ++i)
            {
               if(!FWD)
               {
                  c = 2.0*one_n;
                  if (i == 0 || i == nPts-1)
                  {
                     c *= 0.5;
                  }
               }
               for (Internal::Matrix::Index j = 0; j < nPts; ++j)
               {
                  Internal::MHDFloat c2 = 1.0;
                  if(!FWD && (j == 0 || j == nPts-1))
                  {
                     c2 *= 0.5;
                  }
                  T(i,j) = c*c2*Internal::Math::cos(i*Internal::Math::acos(grid(j)));
               }
            }
            if(FWD)
            {
               T.transposeInPlace();
            }
         }
         else
         {
            throw std::logic_error("Chebyshev Gauss Vandermonde mat not implemented");
         }
      }





   protected:

   private:
      using quadTypeTag_t = quadType;
};

} // namespace Quadrature
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_QUADRATURE_CHEBYSHEVRULE_HPP
