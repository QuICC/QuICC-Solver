/**
 * @file SphereRadialProfilesWriter.cpp
 * @brief Source of the implementation of the ASCII Radial Profiles in a sphere
 */

// Format of the  .dat file produced here:

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/SphereRadialProfilesWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/RadialProfiles.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereRadialProfilesWriter::SphereRadialProfilesWriter(const std::string& prefix, const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction>> pF)
      : IVariableAsciiWriter(prefix + Tags::RadialProfiles::BASENAME, Tags::RadialProfiles::EXTENSION, prefix + Tags::RadialProfiles::HEADER, type, Tags::RadialProfiles::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mPF(std::move(pF))
   {
      assert(this->mPF.size() >= 1 && "SphereRadialProfilesWriter requires at least 1 element in pF vector");
   }

   SphereRadialProfilesWriter::~SphereRadialProfilesWriter()
   {
   }

   void SphereRadialProfilesWriter::init()
   {
      int tmpSize = 2 * this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Generate grid
      int size = std::max(32, tmpSize); // at least 32 points
      Internal::Array weights;
      QuICC::Polynomial::Quadrature::WorlandChebyshevRule quad;
      quad.computeQuadrature(this->mGrid, weights, size);

      IVariableAsciiWriter::init();
   }

   void SphereRadialProfilesWriter::writeContent()
   {
      this->mProfiles.resize(this->mGrid.size(), this->mPF.size());
      for(int i = 0; i < static_cast<int>(this->mPF.size()); i++)
      {
         this->mProfiles.col(i) = this->mPF[i]->evaluate(this->mGrid, 0, 0);
      }
      // Create file
      this->preWrite();

      // Get the "global" profiles from MPI code, not sure if needed
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mGrid.data(), this->mGrid.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->mProfiles.data(), this->mProfiles.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      Array rGrid     = this->mGrid.cast<MHDFloat>();
      Matrix rProfiles = this->mProfiles.cast<MHDFloat>();

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;

         // header
         this->mFile << "#" << "\t" << "Radius" << "\t";
         for(int i = 0; i < static_cast<int>(this->mPF.size()); i++)
         {
            this->mFile << this->mPF[i]->getName() << "\t";
         }
         this->mFile << "\n";

         // data
         for(int i = 0; i < rGrid.size(); i++)
         {
            this->mFile << ioFW(ioPrec) << rGrid(i);
            for(int j = 0; j < rProfiles.cols(); j++)
            {
               this->mFile << "\t" << ioFW(ioPrec) << rProfiles(i, j);
            }
            this->mFile << "\n";
         }
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if profile is NaN
      if(std::isnan(rProfiles.sum()))
      {
         QuICCEnv().abort("Some spherical Radial Profiles are NaN!");
      }
   }

}
}
}
