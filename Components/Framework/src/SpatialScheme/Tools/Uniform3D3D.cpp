/** 
 * @file Uniform3D3D.cpp
 * @brief Source of the tools for uniform + uniform + uniform spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D3D.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int Uniform3D3D::truncationFwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int Uniform3D3D::truncationBwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int Uniform3D3D::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool Uniform3D3D::isOptimal(const int nN, const int maxL)
   {
      return true;
   }

   int Uniform3D3D::totalModes(const int n2D, const int n3D, Splitting::Locations::Id split)
   {
      int nModes = 0;
      if(split == Splitting::Locations::FIRST)
      {
         nModes = n3D;
      }
      else if(split == Splitting::Locations::SECOND)
      {
         nModes = n2D*n3D;
      }
      else if(split == Splitting::Locations::BOTH)
      {
         nModes = n2D;
      }
      else
      {
         throw std::logic_error("Splitting location is not implemented");
      }

      return nModes;
   }

   void Uniform3D3D::buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>&, const std::vector<int>&)
   {
      // Create start indexes and length
      std::vector<int> i0(2, 0);
      std::vector<int> iN(2, 0);

      // Get size of the splittable dimension(s)
      int tot = this->totalModes(n2D, n3D, Splitting::Locations::FIRST);

      // Build a simple balanced split
      this->balancedSplit(i0.at(0), iN.at(0), tot, bins.at(0), id.at(0), false);

      // Get size of the splittable dimension(s)
      tot = this->totalModes(n2D, n3D, Splitting::Locations::BOTH);
      tot *= iN.at(0);

      // Compute a balanced splitting
      this->balancedSplit(i0.at(1), iN.at(1), tot, bins.at(1), id.at(1), false);

      int& p0 = i0.at(0);
      int& pN = iN.at(0);
      int& t0 = i0.at(1);
      int& tN = iN.at(1);

      // Create start indexes and length
      std::vector<int> n0(pN+1,0);
      std::vector<int> nN(pN+1,0);
      n0.at(0) = p0;
      nN.at(0) = pN;

      // Compute starting point of grid points
      for(int i = 0; i < t0; ++i)
      {
         n0.at((i % nN.at(0)) + 1) += 1;
      }

      // Get smallest th0 to use as offset
      t0 = 0;
      int n0Min = n0.at(1);
      for(int r = 1; r < nN.at(0); ++r)
      {
         if(n0Min > n0.at(r+1))
         {
            t0 = r;
            n0Min = n0.at(r+1);
         }
      }

      // Compute optimal number of grid points
      for(int i = 0; i < tN; ++i)
      {
         nN.at(((i + t0) % nN.at(0)) + 1) += 1;
      }

      int k0 = n0.front();
      int kN = nN.front();
      std::vector j0(n0.end()-kN, n0.end());
      std::vector jN(nN.end()-kN, nN.end());
      int c0 = 0;
      int cN = n2D*n3D;

      this->buildMap(modes, k0, kN, j0, jN, c0, cN);
   }

} // Tools
} // SpatialScheme
} // QuICC
