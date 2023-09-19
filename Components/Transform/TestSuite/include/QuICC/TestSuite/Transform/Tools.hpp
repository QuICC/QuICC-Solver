/**
 * @file Tools.hpp
 * @brief Useful utility functions for tests
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_TOOLS_HPP
#define QUICC_TESTSUITE_TRANSFORM_TOOLS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <vector>
#include <tuple>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Precision.hpp"

namespace QuICC {

namespace TestSuite {

namespace Transform {

   /// Typedef for error output
   typedef std::vector<std::pair<MHDFloat,int> > ErrorDataType;

   /// Typedef for filename descritpion
   typedef std::vector<std::tuple<std::string,MHDFloat,bool> > FilenameExtDescr;

   /**
    * @brief Build filename extension
    */
   std::string filenameExt(const FilenameExtDescr& descr);

   template <typename TData> ErrorDataType compareData(const TData& outData, const TData& refData, const MHDFloat relThreshold);

   /// Compute error
   template <typename TData> ErrorDataType compareData(const TData& outData, const TData& refData, const MHDFloat relThreshold)
   {
      // Compute relative error. Subnormal values are forced to optimal precision (epsilon)
      // as relative errors doesn't make sense in this case
      ErrorDataType error;
      for(int j = 0; j < refData.cols(); j++)
      {
         std::pair<MHDFloat,int> err = {std::numeric_limits<MHDFloat>::epsilon(),-1};
         for(int i = 0; i < refData.rows(); i++)
         {
            MHDFloat e;
            if(std::abs(refData(i,j)) == 0.0 || (std::abs(refData(i,j)) < relThreshold))
            {
               e = std::abs((refData(i,j)-outData(i,j)));
            } else
            {
               e = std::abs((refData(i,j)-outData(i,j))/refData(i,j));
            }
            if(std::abs(refData(i,j)) < std::numeric_limits<MHDFloat>::min())
            {
               e = std::numeric_limits<MHDFloat>::epsilon();
            }
            if(e >= err.first)
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

#endif // QUICC_TESTSUITE_TRANSFORM_TOOLS_HPP
