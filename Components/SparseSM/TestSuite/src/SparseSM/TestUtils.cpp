/**
 * @file TestUtils.cpp
 * @brief Source of test utils
 */

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/SparseSM/TestUtils.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

/**
 * @brief computing ULP
 */
std::pair<bool, MHDFloat> computeUlp(const ::QuICC::MHDFloat data,
   const ::QuICC::MHDFloat ref, const ::QuICC::MHDFloat refMod,
   const ::QuICC::MHDFloat maxUlp, const ::QuICC::MHDFloat epsilon)
{
   bool isEqual = false;

   auto diff = std::abs(data - ref);
   auto tol = maxUlp * epsilon;

   if (diff < tol)
   {
      isEqual = diff < (tol * refMod);
   }
   else
   {
      isEqual = (diff / refMod) < tol;
   }

   auto ulp = diff / (refMod * epsilon);

   return std::make_pair(isEqual, ulp);
};

/**
 * @brief Check stencil
 */
void checkStencil(const Matrix& bc, const Matrix& stencil,
   const MHDFloat maxUlp)
{
   ::QuICC::MHDFloat eps = std::numeric_limits<::QuICC::MHDFloat>::epsilon();

   for (int j = 0; j < bc.rows(); j++)
   {
      INFO("BC row: " << j);

      Array b = bc.row(j).transpose();
      for (int i = 0; i < stencil.cols(); i++)
      {
         INFO("Stencil col: " << i);

         Array bs = b.array() * stencil.col(i).array();
         MHDFloat refMod = bs.array().abs().maxCoeff();
         MHDFloat err = std::abs(bs.sum());
         auto e = computeUlp(err, 0.0, refMod, maxUlp, eps);

         INFO("Measured ulp: " << e.second);
         CHECK(e.first);
      }
   }
}

} // namespace SparseSM
} // namespace TestSuite
} // namespace QuICC
