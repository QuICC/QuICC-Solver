/** 
 * @file Triangular.hpp
 * @brief Implementation of the tools for the triangular triangular schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_TRIANGULAR_HPP
#define QUICC_SPATIALSCHEME_TOOLS_TRIANGULAR_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the triangular spatial schemes
    */
   class Triangular
   {
      public:
         /**
          * @brief Compute triangular truncation
          */
         static int truncation(const int nN, const int l);

         /**
          * @brief Check if chosen resolution is optimal
          */
         static bool isOptimal(const int nN, const int maxL);

         /**
          * @brief Build map of indexes for a triangular spatial schemes
          */
         static void buildMap(std::multimap<int,int>& modes, const int i0, const int iN, const ArrayI& j0, const ArrayI& jN, const int c0, const int cN);

         /**
          * @brief Fill the indexes for 2D and 3D for triangular spatial schemes
          */
         static void fillIndexes2D3D(std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::multimap<int,int>& modes);

         /**
          * @brief Fill the indexes for 1D for triangular spatial schemes
          */
         static void fillIndexes1D(std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, const ArrayI& idx3D, const int nF1D, const int nB1D);

   };
}
}
}

#endif // QUICC_SPATIALSCHEME_TOOLS_TRIANGULAR_HPP
