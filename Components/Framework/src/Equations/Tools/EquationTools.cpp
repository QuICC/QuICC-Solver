/**
 * @file EquationTools.cpp
 * @brief Source of the requirement tools to work with variables and equations
 */

// System includes
//
#include <cassert>
#include <algorithm>

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/EquationTools.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Tools {

   SparseMatrix makeBlockMatrix(const int nBlocks, const int row, const int col)
   {
      SparseMatrix  mat(nBlocks, nBlocks);

      mat.insert(row, col) = 1;
      mat.makeCompressed();

      return mat;
   }

   void setupPhysicalKernels(std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& rKernels, const std::vector<Array>& mesh)
   {
      auto spMesh = std::make_shared<std::vector<Array> >(mesh);

      for(auto &k: rKernels)
      {
         k.second->setMesh(spMesh);
      }
   }
}
}
}
