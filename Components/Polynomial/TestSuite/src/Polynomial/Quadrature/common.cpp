/** 
 * @file common.cpp
 * @brief Source of the useful functions for transform tests
 */


// Debug includes
//

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <iostream>

// External includes
//

// Project includes
//
#include "QuICC/TestSuite/Polynomial/Quadrature/common.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Quadrature {

   std::string resname(const int physN, const std::map<std::string,MHDFloat>& params)
   {
      std::string s = "";

      for(const auto& p: params)
      {
         s += "_" + p.first + std::to_string(p.second);
      }

      s += "_g" + std::to_string(physN);

      return s;
   }

   ErrorDataType computeError(const int physN, const std::map<std::string,MHDFloat>& params, const Matrix& outData, const std::string& fname)
   {
      // Read reference data
      Matrix refData(outData.rows(), outData.cols());
      std::string ext = resname(physN, params);
      std::string filename = fname;
      std::string fullname = "_refdata/Polynomial/Quadrature/" + filename.insert(fname.length()-4, ext);
      readData(refData, fullname);

      INFO(fullname);

      // Compute relative error. Subnormal values are forced to optimal precision (epsilon) 
      // as relative errors doesn't make sense in this case
      ErrorDataType error;

      constexpr auto epsilon = std::numeric_limits<MHDFloat>::epsilon();
      constexpr uint ulp = 11;
      constexpr auto tol = ulp * epsilon;

      // j = 0: grid
      // j = 1: weight
      for(int j = 0; j < refData.cols(); j++)
      {
         std::pair<MHDFloat,int> err = {0.0,-1};
         for(int i = 0; i < refData.rows(); i++)
         {
            auto diff = std::abs(refData(i,j)-outData(i,j));
            auto ref = std::abs(refData(i,j));

            // when comparing to zero there is no reference scale
            // use closest non-zero point instead
            if(ref == 0.0)
            {
               int shift = -1;
               if (i == 0)
               {
                  shift = 1;
               }
               ref = std::abs(refData(i+shift,j));
            }

            // check for that is a normalized value
            REQUIRE(ref >= std::numeric_limits<MHDFloat>::min());


            bool isEqual = false;
            if(diff < tol)
            {
               isEqual = diff < (tol * ref);
            }
            else
            {
               isEqual = (diff / ref ) < tol;
            }

            INFO( "grid(0) or weight(1): " << j );
            INFO( "position: " << i );
            INFO( "refData: " << std::scientific << std::setprecision(16) << refData(i,0) << '\t' << refData(i,1) );
            INFO( "outData: " << std::scientific << std::setprecision(16) << outData(i,0) << '\t' << outData(i,1) );
            INFO( "max ulp: " << ulp);
            INFO( "measured ulp: " << diff / (ref * epsilon));
            CHECK(isEqual);

            MHDFloat e = std::abs(diff/refData(i,j));


            if(std::abs(refData(i,j)) < std::numeric_limits<MHDFloat>::min())
            {
               e = std::numeric_limits<MHDFloat>::epsilon();
            }
            if(e > err.first)
            {
               err = std::make_pair(e, i);
            }
         }
         error.push_back(err);
      }

      return error;
   }

}
}
}
}
