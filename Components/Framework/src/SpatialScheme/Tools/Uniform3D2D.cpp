/** 
 * @file Uniform3D2D.cpp
 * @brief Source of the tools for uniform + uniform + uniform spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D2D.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int Uniform3D2D::truncationFwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int Uniform3D2D::truncationBwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int Uniform3D2D::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool Uniform3D2D::isOptimal(const int nN, const int maxL)
   {
      return true;
   }

   int Uniform3D2D::totalModes(const int n2D, const int n3D, Splitting::Locations::Id split)
   {
      int nModes = 0;
      if(split == Splitting::Locations::FIRST)
      {
         nModes = n2D;
      }
      else if(split == Splitting::Locations::SECOND)
      {
         nModes = n3D;
      }
      else if(split == Splitting::Locations::BOTH)
      {
         nModes = n3D;
      }
      else
      {
         throw std::logic_error("Splitting location is not implemented");
      }

      return nModes;
   }

   void Uniform3D2D::buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>&, const std::vector<int>&)
   {
      // Create start indexes and length
      std::vector<int> i0(2, 0);
      std::vector<int> iN(2, 0);

      // Get size of the splittable dimension(s)
      int tot = this->totalModes(n2D, n3D, Splitting::Locations::BOTH);

      // Build a simple balanced split
      this->balancedSplit(i0.at(0), iN.at(0), tot, bins.at(1), id.at(1), false);

      // Get size of the splittable dimension(s)
      tot = this->totalModes(n2D, n3D, Splitting::Locations::FIRST);

      // Compute a balanced splitting
      this->balancedSplit(i0.at(1), iN.at(1), tot, bins.at(0), id.at(0), false);

      int k0 = i0.front();
      int kN = iN.front();
      std::vector j0(kN, i0.at(1));
      std::vector jN(kN, iN.at(1));
      int c0 = 0;
      int cN = n2D*n3D;

      this->buildMap(modes, k0, kN, j0, jN, c0, cN);
   }

} // Tools
} // SpatialScheme
} // QuICC
