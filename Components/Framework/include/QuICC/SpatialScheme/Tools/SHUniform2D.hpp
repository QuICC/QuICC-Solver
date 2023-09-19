/**
 * @file SHUniform2D.hpp
 * @brief Implementation of the tools for the spherical harmonic + uniform + uniform based schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_SHUNIFORM2D_HPP
#define QUICC_SPATIALSCHEME_TOOLS_SHUNIFORM2D_HPP

// System includes
//
#include <vector>
#include <deque>
#include <queue>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the spherical harmonic + uniform + uniform based schemes
    */
   class SHUniform2D: public IBase
   {
      public:
         /**
          * @brief Default ctor
          *
          * @param nL   Number of harmonic degrees
          */
         SHUniform2D(const int nL);

         /**
          * @brief Default dtor
          */
         ~SHUniform2D() = default;

         /**
          * @brief Compute forward truncation
          *
          * @param nN   Reference truncation
          * @param j    second dimension
          * @param k    third dimension
          */
         int truncationFwd(const int nN, const int j, const int k) final;

         /**
          * @brief Compute backward truncation
          *
          * @param nN   Reference truncation
          * @param j    second dimension
          * @param k    third dimension
          */
         int truncationBwd(const int nN, const int j, const int k) final;

         /**
          * @brief Compute index
          *
          * @param nN   Reference truncation
          * @param j    second dimension
          * @param k    third dimension
          */
         int index(const int nN, const int j, const int k) final;

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

      protected:
         /**
          * @brief Number of harmonic degrees
          */
         int mNl;
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_SHUNIFORM2D_HPP
