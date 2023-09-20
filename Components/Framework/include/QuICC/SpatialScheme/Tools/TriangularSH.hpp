/**
 * @file TriangularSH.hpp
 * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_TRIANGULARSH_HPP
#define QUICC_SPATIALSCHEME_TOOLS_TRIANGULARSH_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBaseSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes
    */
   class TriangularSH: public IBaseSH
   {
      public:
         /**
          * @brief Default ctor
          */
         TriangularSH() = default;

         /**
          * @brief Default dtor
          */
         ~TriangularSH() = default;

         /**
          * @brief Compute triangular forward truncation
          */
         int truncationFwd(const int nN, const int j, const int k) final;

         /**
          * @brief Compute triangular backward truncation
          */
         int truncationBwd(const int nN, const int j, const int k) final;

         /**
          * @brief Compute index
          */
         int index(const int nN, const int j, const int k) final;

         /**
          * @brief Check if chosen resolution is optimal
          */
         bool isOptimal(const int nN, const int maxL) final;

      private:
         /**
          * @brief Minimal truncation for highest modes
          */
         static constexpr const int MIN_TRUNCATION = 3;
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_TRIANGULARSH_HPP
