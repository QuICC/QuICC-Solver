/**
 * @file Operator.cpp
 * @brief Implementation of generic operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Operator.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   Operator::Operator(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : mOrder(order), mZtop(zTop), mZbot(zBot)
   {
   }

   Operator::Operator(const std::size_t order)
      : Operator(order, 0, 0)
   {
   }

   std::size_t Operator::order2Stencil(const std::size_t m, const std::size_t n, const bool isCentral) const
   {
      std::size_t s;
      if(isCentral)
      {
         s = 2*((m + 1)/2) - 1 + n;
      }
      else
      {
         s = m + n;
      }

      return s;
   }

   void Operator::fdWeights(std::vector<std::vector<Internal::MHDFloat> >& w, const Internal::MHDFloat z, const std::vector<Internal::MHDFloat>& x, const std::size_t m) const
   {
      using namespace Internal::Literals;

      assert(w.size() == 0);

      // Initialize storage
      for(std::size_t i = 0; i <= m; i++)
      {
         w.emplace_back(x.size(), 0_mp);
      }

      Internal::MHDFloat c1,c2,c3,c4,c5;

      c1 = 1.0_mp;
      c4 = x.at(0) - z;
      w[0][0] = 1.0_mp;
      for(std::size_t i = 1; i < x.size(); i++)
      {
         std::size_t mn = std::min(i, m);
         c2 = 1.0_mp;
         c5 = c4;
         c4 = x[i] - z;
         for(std::size_t j = 0; j < i; j++)
         {
            c3 = x[i] - x[j];
            c2 *= c3;
            if(j == i-1)
            {
               for(int k = mn; k > 0; k--)
               {
                  w[k][i] = c1*(k*w[k-1][i-1] - c5*w[k][i-1])/c2;
               }
               w[0][i] = -c1*c5*w[0][i-1]/c2;
            }
            for(int k = mn; k > 0; k--)
            {
               w[k][j] = (c4*w[k][j] - k*w[k-1][j])/c3;
            }
            w[0][j] = c4*w[0][j]/c3;
         }
         c1 = c2;
      }
   }

   void Operator::fdMatrices(std::vector<Internal::SparseMatrix>& wMat, const Internal::Array& grid, const std::size_t n, const std::size_t m) const
   {

      assert(wMat.size() == 0);
      assert(grid.size() >= static_cast<int>(this->order2Stencil(m, n, false)));

      std::size_t nR = static_cast<std::size_t>(grid.size());

      // Initialize storage
      for(std::size_t p = 0; p <= m; p++)
      {
         wMat.emplace_back(nR, nR);
      }

      for(std::size_t p = 0; p < nR; p++)
      {
         std::size_t s_ = 0;
         const auto& z = grid(p);
         std::vector<Internal::MHDFloat> x;
         std::vector<std::vector<Internal::MHDFloat>> w;
         int i0;
         for(std::size_t j = m; j > 0; j--)
         {
            std::size_t s = this->order2Stencil(j, n, true);
            if(p < s/2)
            {
               s = this->order2Stencil(j, n, false);
               i0 = 0;
            }
            else if((p + s/2) > (nR - 1))
            {
               s = this->order2Stencil(j, n, false);
               i0 = nR-s;
            }
            else
            {
               assert(s % 2 == 1);
               i0 = std::max<int>(0, (p - s/2));
            }
            
            if(s != s_)
            {
               x.clear();
               w.clear();
               assert(s > 1);
               for(std::size_t i = 0; i < s; i++)
               {
                  x.push_back(grid[i0+i]);
               }

               this->fdWeights(w, z, x, m);
            }

            Internal::SparseMatrix row(nR, nR);
            std::vector<Eigen::Triplet<Internal::MHDFloat>> coeffs;
            for(std::size_t i = 0; i < s; i++)
            {
               coeffs.emplace_back(p, i+i0, w[j][i]);
            }

            row.setFromTriplets(coeffs.begin(), coeffs.end());
            row.makeCompressed();
            wMat.at(j).makeCompressed();
            wMat.at(j) += row;

            // Store previous stencil size
            s_ = s;
         }

         // identity for 0th order derivative
         Internal::SparseMatrix row(nR, nR);
         std::vector<Eigen::Triplet<Internal::MHDFloat>> coeffs;
         coeffs.emplace_back(p, p, 1);

         row.setFromTriplets(coeffs.begin(), coeffs.end());
         wMat.at(0) += row;
      }
         
      for(std::size_t j = 0; j <= m; j++)
      {
         wMat.at(j).makeCompressed();
      }
   }

   Internal::SparseMatrix Operator::zeroTopBottom(const int nR, const int zTop, const int zBot) const
   {
      using namespace Internal::Literals;

      Internal::SparseMatrix op(nR, nR);

      std::vector<Eigen::Triplet<Internal::MHDFloat>> coeffs;
      for(int i = zTop; i < nR - zBot; i++)
      {
         coeffs.emplace_back(i,i, 1_mp);
      }

      op.setFromTriplets(coeffs.begin(), coeffs.end());

      return op;
   }

   Internal::SparseMatrix Operator::zeroTopBottom(const int nR) const
   {
      return this->zeroTopBottom(nR, static_cast<int>(this->mZtop), static_cast<int>(this->mZbot));
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
