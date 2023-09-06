/** 
 * @file UniformSH.hpp
 * @brief Implementation of the tools for the uniform + spherical harmonics schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_UNIFORMSH_HPP
#define QUICC_SPATIALSCHEME_TOOLS_UNIFORMSH_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBaseSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the uniform + spherical harmonic spatial schemes
    */
   class UniformSH: public IBaseSH
   {
      public:
         /**
          * @brief Default ctor
          */
         UniformSH() = default;

         /**
          * @brief Default dtor
          */
         ~UniformSH() = default;

         /**
          * @brief Compute uniform forward truncation
          */
         int truncationFwd(const int nN, const int l) final;

         /**
          * @brief Compute uniform backward truncation
          */
         int truncationBwd(const int nN, const int l) final;

         /**
          * @brief Compute index
          */
         int index(const int nN, const int k) final;

         /**
          * @brief Check if chosen resolution is optimal
          */
         bool isOptimal(const int nN, const int maxL) final;

      protected:
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_UNIFORMSH_HPP
