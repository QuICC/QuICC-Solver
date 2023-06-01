/** 
 * @file WFTMesher.cpp
 * @brief Source of the WFT spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/WFTMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   WFTMesher::WFTMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNr(-1), mNt(-1), mNz(-1)
   {
   }

   void WFTMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& N = this->mDims.at(0);
      int& M = this->mDims.at(1);
      int& Z = this->mDims.at(2);

      // Get dealiased Worland size
      this->mNr = Transform::Poly::Tools::dealias(N + M/2 + 8);

      // Get mixed dealiased FFT size
      this->mNt = Transform::Fft::Tools::dealiasMixedFft(M+1);
      // Check for optimised FFT sizes
      this->mNt = Transform::Fft::Tools::optimizeFft(this->mNt);

      // Get standard dealiased FFT size
      this->mNz = Transform::Fft::Tools::dealiasCosFft(Z+1);
      // Check for optimised FFT sizes
      this->mNz = Transform::Fft::Tools::optimizeFft(this->mNz);
   }

   int WFTMesher::nPhys1D() const
   {
      return this->mNr;
   }

   int WFTMesher::nPhys2D() const
   {
      return this->mNt;
   }

   int WFTMesher::nPhys3D() const
   {
      return this->mNz;
   }

   int WFTMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int WFTMesher::nSpec2D() const
   {
      const int& M = this->mDims.at(1);
      return M + 1;
   }

   int WFTMesher::nSpec3D() const
   {
      const int& Z = this->mDims.at(2);
      return Z + 1;
   }

   int WFTMesher::nDealias1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int WFTMesher::nDealias2D() const
   {
      return this->mNt/2 + 1;
   }

   int WFTMesher::nDealias3D() const
   {
      return this->mNz;
   }

} // SpatialScheme
} // QuICC
