/** 
 * @file Parity.hpp
 * @brief Implementation of the tools for parity splitting
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_PARITY_HPP
#define QUICC_SPATIALSCHEME_TOOLS_PARITY_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for parity splitting
    */
   class Parity
   {
      public:
         /**
          * @brief Build map of even and odd harmonics for L ordering
          */
         static void splitParityL(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& blockSize, MatrixI& evenBlocks, MatrixI& oddBlocks);

         /**
          * @brief Build map of even and odd harmonics for M ordering
          */
         static void splitParityM(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& blockSize, MatrixI& evenBlocks, MatrixI& oddBlocks);

   };
} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_PARITY_HPP
