/**
 * @file common.hpp
 * @brief Functions common to all quadrature tests
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_COMMON_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_COMMON_HPP

// Configuration includes
//

// System includes
//
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdexcept>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Casts.hpp"
#include "QuICC/TestSuite/Polynomial/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Quadrature {

   /// Typedef for error output
   typedef std::vector<std::pair<MHDFloat,int> > ErrorDataType;

   /**
    * @brief Build filenmae extension with resolution information
    */
   std::string resname(const int physN, const std::map<std::string,MHDFloat>& params);

   ///
   template <typename TOp, int TN> ErrorDataType maxOperatorError(const int size, const std::map<std::string,MHDFloat>& params, const std::string& fname, const bool keepData);

   template <typename TQuad> Matrix computeRule(const int size, const std::string& fname);

   template <typename TQuad> Matrix computeRule(const int size, const std::string& fname, const MHDFloat a, const MHDFloat b);

   /// Compute error with respect to reference data
   ErrorDataType computeError(const int physN, const std::map<std::string,MHDFloat>& params, const Matrix& outData, const std::string& fname);

   template <typename TQuad, int TN> ErrorDataType maxOperatorError(const int physN, const std::map<std::string,MHDFloat>& params, const std::string& fname, const bool keepData)
   {
      Matrix outData;
      std::string foutname = "";

      // Output data should be written to file
      if(keepData)
      {
         foutname = fname;
      }

      if constexpr(TN == 0)
      {
         assert(params.size() == 0);
         outData = computeRule<TQuad>(physN, foutname);
      } else if(TN == 2)
      {
         assert(params.size() == 2);
         outData = computeRule<TQuad>(physN, foutname, params.find("a")->second, params.find("b")->second);
      } else
      {
         throw std::logic_error("Number of parameters doesn't match existing quadrature rule");
      }

      // Compare data to reference
      ErrorDataType error = computeError(physN, params, outData, fname);

      return error;
   }

   template <typename TQuad> Matrix computeRule(const int size, const std::string& fname)
   {
      Internal::Array igrid;
      Internal::Array iweights;
      TQuad quad;
      quad.computeQuadrature(igrid, iweights, size);

      Matrix outData(igrid.size(),2);
      outData.col(0) = igrid.cast<MHDFloat>();
      outData.col(1) = iweights.cast<MHDFloat>();

      if(fname != "")
      {
         std::map<std::string,MHDFloat> mapped;
         std::string ext = resname(size, mapped);
         std::string filename = fname;
         writeData("_data/Polynomial/Quadrature/" + filename.insert(fname.length()-4, ext), outData);
      }

      return outData;
   }

   template <typename TQuad> Matrix computeRule(const int size, const std::string& fname, const MHDFloat a, const MHDFloat b)
   {
      Internal::Array igrid;
      Internal::Array iweights;
      TQuad quad(a, b);
      quad.computeQuadrature(igrid, iweights, size);

      Matrix outData(igrid.size(),2);
      outData.col(0) = igrid.cast<MHDFloat>();
      outData.col(1) = iweights.cast<MHDFloat>();


      if(fname != "")
      {
         std::map<std::string,MHDFloat> mapped;
         mapped.insert(std::make_pair("a",a));
         mapped.insert(std::make_pair("b",b));
         std::string ext = resname(size, mapped);
         std::string filename = fname;
         writeData("_data/Polynomial/Quadrature/" + filename.insert(fname.length()-4, ext), outData);
      }

      return outData;
   }

}
}
}
}

#endif //QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_COMMON_HPP
