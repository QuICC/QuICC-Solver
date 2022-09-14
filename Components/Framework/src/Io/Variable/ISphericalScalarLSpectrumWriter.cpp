/**
 * @file ISphericalScalarLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics L energy spectrum calculation for scalar field in a spherical geometry
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
#include "QuICC/Io/Variable/ISphericalScalarLSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/Spectrum.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalScalarLSpectrumWriter::ISphericalScalarLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarEnergyBaseWriter(prefix + Tags::Spectrum::LBASENAME, Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type, Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mEnergy(0)
   {
   }

   ISphericalScalarLSpectrumWriter::~ISphericalScalarLSpectrumWriter()
   {
   }

   void ISphericalScalarLSpectrumWriter::init()
   {
      this->mEnergy = Array::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));

      ISphericalScalarEnergyBaseWriter::init();
   }

   void ISphericalScalarLSpectrumWriter::resetEnergy()
   {
      this->mEnergy.setZero();
   }

   void ISphericalScalarLSpectrumWriter::storeEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mEnergy(l) += energy;
   }

   void ISphericalScalarLSpectrumWriter::writeContent()
   {
      // Normalize by the volume
      this->mEnergy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mEnergy.data(), this->mEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      using Tools::Formatter::ioIW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << "# time: " << std::setprecision(ioPrec);
         this->mFile << ioFW(ioPrec) << this->mTime << std::endl;

         this->mFile << "# energy: " << std::setprecision(ioPrec);
         MHDFloat t = this->mEnergy.sum();
         this->mFile << ioFW(ioPrec) << t << std::endl;

         this->mFile << "# SD: standard deviation" << std::endl;
         this->mFile << "# characteristic L: " << std::setprecision(ioPrec);
         Array w = Array::LinSpaced(this->mEnergy.size(), 0, this->mEnergy.size()-1);
         MHDFloat e = (w.array()*this->mEnergy.array()).sum()/t;
         MHDFloat v = ((w.array() - e).array().pow(2)*this->mEnergy.array()).sum()/t;
         this->mFile << std::setprecision(2) << std::fixed;
         this->mFile << e << " (SD=" << std::sqrt(v) << ")" << std::left;
         this->mFile << std::endl;
         this->mFile << std::scientific;

         this->mFile << std::left << ioIW() << "# l" << "\t";
         this->mFile << ioFW(ioPrec) << "total" << std::endl;

         // Total
         for(int i = 0; i < this->mEnergy.size(); i++)
         {
            this->mFile << std::left << ioIW() << i << "\t" << std::setprecision(ioPrec);
            this->mFile << ioFW(ioPrec) << this->mEnergy(i);
            this->mFile << std::endl;
         }

         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mEnergy.sum()))
      {
         QuICCEnv().abort("Scalar energy is NaN!");
      }
   }

}
}
}
