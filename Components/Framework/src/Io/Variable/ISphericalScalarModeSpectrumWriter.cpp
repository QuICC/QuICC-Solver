/**
 * @file ISphericalScalarModeSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics mode
 * energy spectrum calculation for scalar field in a spherical geometry
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/ISphericalScalarModeSpectrumWriter.hpp"
#include "QuICC/Io/Variable/Tags/Spectrum.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

ISphericalScalarModeSpectrumWriter::ISphericalScalarModeSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalScalarEnergyBaseWriter(prefix + Tags::Spectrum::MODEBASENAME,
       Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type,
       Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE),
    mEnergy(0, 0)
{}

void ISphericalScalarModeSpectrumWriter::init()
{
   this->mEnergy =
      Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,
                      Dimensions::Space::SPECTRAL),
         this->res().sim().dim(Dimensions::Simulation::SIM2D,
            Dimensions::Space::SPECTRAL));

   ISphericalScalarEnergyBaseWriter::init();
}

void ISphericalScalarModeSpectrumWriter::resetEnergy()
{
   this->mEnergy.setZero();
}

void ISphericalScalarModeSpectrumWriter::storeEnergy(const int l, const int m,
   const MHDFloat energy)
{
   if (m <= l)
   {
      this->mEnergy(l, m) += energy;
   }
}

void ISphericalScalarModeSpectrumWriter::writeContent()
{
   // Normalize by the volume
   this->mEnergy /= this->mVolume;

   // Create file
   this->preWrite();

// Get the "global" Kinetic energy from MPI code
#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, this->mEnergy.data(), this->mEnergy.size(),
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif // QUICC_MPI

   using Tools::Formatter::ioFW;
   using Tools::Formatter::ioIW;
   int ioPrec = 14;

   // Check if the workflow allows IO to be performed
   if (QuICCEnv().allowsIO())
   {
      this->mFile << std::scientific;
      this->mFile << "# time: " << std::setprecision(ioPrec);
      this->mFile << ioFW(ioPrec) << this->mTime << std::endl;

      this->mFile << "# energy: " << std::setprecision(ioPrec);
      MHDFloat t = this->mEnergy.sum();
      this->mFile << ioFW(ioPrec) << t << std::endl;

      this->mFile << std::left << ioIW() << "# l"
                  << "\t"
                  << "m"
                  << "\t";
      this->mFile << ioFW(ioPrec) << "total" << std::endl;

      // Total
      for (int l = 0; l < this->mEnergy.rows(); l++)
      {
         for (int m = 0; m <= l; m++)
         {
            this->mFile << std::left << ioIW() << l << "\t" << m << "\t"
                        << std::setprecision(ioPrec);
            this->mFile << ioFW(ioPrec) << this->mEnergy(l, m);
            this->mFile << std::endl;
         }
      }
   }

   // Close file
   this->postWrite();

   // Abort if kinetic energy is NaN
   if (std::isnan(this->mEnergy.sum()))
   {
      QuICCEnv().abort("Scalar  modeenergy spectrum is NaN!");
   }
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
