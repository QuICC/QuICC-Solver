/**
 * @file ISphericalTorPolEnstrophyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics enstrophy calculation for toroidal/poloidal field in a spherical geometry
 */

// Configuration includes
//

// System includes
//
#include <iomanip>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/ISphericalTorPolEnstrophyWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Enstrophy.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalTorPolEnstrophyWriter::ISphericalTorPolEnstrophyWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyBaseWriter(prefix + Tags::Enstrophy::BASENAME, Tags::Enstrophy::EXTENSION, prefix + Tags::Enstrophy::HEADER, type, Tags::Enstrophy::VERSION, Dimensions::Space::SPECTRAL), mTorEnstrophy(2), mPolEnstrophy(2)
   {
      this->mTorEnstrophy.setConstant(-1);
      this->mPolEnstrophy.setConstant(-1);
   }

   ISphericalTorPolEnstrophyWriter::~ISphericalTorPolEnstrophyWriter()
   {
   }

   void ISphericalTorPolEnstrophyWriter::initializeEnstrophy()
   {
      this->mTorEnstrophy.setZero();
      this->mPolEnstrophy.setZero();
   }

   void ISphericalTorPolEnstrophyWriter::storeQEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      if((l - m)%2 == 0)
      {
         this->mTorEnstrophy(0) += enstrophy;
      } else
      {
         this->mTorEnstrophy(1) += enstrophy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::storeSEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      if((l - m)%2 == 0)
      {
         this->mTorEnstrophy(0) += enstrophy;
      } else
      {
         this->mTorEnstrophy(1) += enstrophy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::storeTEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      if((l - m)%2 == 1)
      {
         this->mPolEnstrophy(0) += enstrophy;
      } else
      {
         this->mPolEnstrophy(1) += enstrophy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::writeContent()
   {
      // Normalize by the volume
      this->mTorEnstrophy /= this->mVolume;
      this->mPolEnstrophy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic enstrophy from MPI code
      #ifdef QUICC_MPI
         Array enstrophy(4);

         enstrophy(0) = this->mTorEnstrophy(0);
         enstrophy(1) = this->mTorEnstrophy(1);
         enstrophy(2) = this->mPolEnstrophy(0);
         enstrophy(3) = this->mPolEnstrophy(1);

         MPI_Allreduce(MPI_IN_PLACE, enstrophy.data(), enstrophy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnstrophy(0) = enstrophy(0);
         this->mTorEnstrophy(1) = enstrophy(1);
         this->mPolEnstrophy(0) = enstrophy(2);
         this->mPolEnstrophy(1) = enstrophy(3);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         // Total
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mTorEnstrophy.sum() + this->mPolEnstrophy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mTorEnstrophy(0) + this->mPolEnstrophy(0) << "\t" << ioFW(ioPrec) << this->mTorEnstrophy(1) + this->mPolEnstrophy(1);
         }

         // Toroidal
         this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mTorEnstrophy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mTorEnstrophy(0) << "\t" << ioFW(ioPrec) << this->mTorEnstrophy(1);
         }

         // Poloidal
         this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mPolEnstrophy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mPolEnstrophy(0) << "\t" << ioFW(ioPrec) << this->mPolEnstrophy(1);
         }

         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic enstrophy is NaN
      if(std::isnan(this->mTorEnstrophy.sum()) || std::isnan(this->mPolEnstrophy.sum()))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Toroidal/Poloidal enstrophy is NaN!");
      }
   }

}
}
}
