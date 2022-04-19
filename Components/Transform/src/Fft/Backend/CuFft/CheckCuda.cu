/**
 * @file CheckCuda.cu
 * @brief Defines some helper functions for checking CUDA errors
 */

// System includes
//
#include <iostream>
#include <string>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/CheckCuda.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

void CheckCuda(cudaError_t status, const int line, const bool useException)
{
   if (status != cudaSuccess)
   {
      std::string msg;
      msg  = "CUDA API failed at line " + std::to_string(line);
      msg += " with error: " + std::string(cudaGetErrorString(status));
      msg += " (" + std::to_string(status) + ")";

      if(useException)
      {
        throw std::logic_error(msg);
      } else
      {
        std::cerr << msg << std::endl;
      }
   }
}

void CheckCuda(cublasStatus_t status, const int line, const bool useException)
{
   if (status != CUBLAS_STATUS_SUCCESS)
   {
      std::string msg;
      msg  = "CUBLAS API failed at line " + std::to_string(line);
      msg += " with error:";
      msg += " (" + std::to_string(status) + ")";

      if(useException)
      {
        throw std::logic_error(msg);
      } else
      {
        std::cerr << msg << std::endl;
      }
   }
}

void CheckCuda(cusparseStatus_t status, const int line, const bool useException)
{
   if (status != CUSPARSE_STATUS_SUCCESS)
   {
      std::string msg;
      msg  = "CUSPARSE API failed at line " + std::to_string(line);
      msg += " with error: " + std::string(cusparseGetErrorString(status));
      msg += " (" + std::to_string(status) + ")";

      if(useException)
      {
        throw std::logic_error(msg);
      } else
      {
        std::cerr << msg << std::endl;
      }
   }
}


}
}
}
}
}
