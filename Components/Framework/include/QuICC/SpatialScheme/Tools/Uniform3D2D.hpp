/** 
 * @file Uniform3D2D.hpp
 * @brief Implementation of the tools for the uniform + uniform + uniform schemes (3D), 2D transform stage
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_UNIFORM3D2D_HPP
#define QUICC_SPATIALSCHEME_TOOLS_UNIFORM3D2D_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the uniform + uniform + uniform spatial schemes (3D), 2D transform stage
    */
   class Uniform3D2D: public IBase
   {
      public:
         /**
          * @brief Default ctor
          */
         Uniform3D2D() = default;

         /**
          * @brief Default dtor
          */
         ~Uniform3D2D() = default;

         /**
          * @brief Compute forward truncation
          */
         int truncationFwd(const int nN, const int l) final;

         /**
          * @brief Compute backward truncation
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

         /**
          * @brief Compute total number of modes
          */
         int totalModes(const int n2D, const int n3D, const Splitting::Locations::Id split) final;

         /**
          * @brief Build map of balanced indexes for a generic spatial schemes
          */
         void buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>& n0, const std::vector<int>& nN) final;
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_UNIFORM3D2D_HPP
