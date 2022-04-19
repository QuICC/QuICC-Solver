/**
 * @file Cartesian1DScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) energy calculation for scalar field
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
#include "QuICC/Io/Variable/Cartesian1DScalarEnergyWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/EnergyTags.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian1DScalarEnergyWriter::Cartesian1DScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Energy::BASENAME, Tags::Energy::EXTENSION, prefix + Tags::Energy::HEADER, type, Tags::Energy::VERSION, Dimensions::Space::SPECTRAL), mEnergy(-Array::Ones(2))
   {
   }

   Cartesian1DScalarEnergyWriter::~Cartesian1DScalarEnergyWriter()
   {
   }

   void Cartesian1DScalarEnergyWriter::init()
   {
      // Normalize by Cartesian volume V = (2*pi*Box1D/k1D)*(2*pi*Box2D/k2D)*(zo-zi) but FFT already includes 1/(2*pi)
      MHDFloat zi = this->mPhysical.find(NonDimensional::Lower1D::id())->second->value();
      MHDFloat zo = this->mPhysical.find(NonDimensional::Upper1D::id())->second->value();
      this->mVolume = (zo-zi)/(this->res().sim().boxScale(Dimensions::Simulation::SIM2D)*this->res().sim().boxScale(Dimensions::Simulation::SIM3D));

      IVariableAsciiWriter::init();
   }

   void Cartesian1DScalarEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      // Initialize the energy
      this->mEnergy.setConstant(0.0);

      // Dealias variable data
      coord.communicator().dealiasSpectral(sRange.first->second->rDom(0).rTotal());

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute energy reduction
      Matrix spectrum(rInVar.data().cols(), 1);
      coord.transform1D().reduce(spectrum, rInVar.data(), Transform::Reductor::Energy::id());

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      int start = 0;
      for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         int n = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k);
         if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
         {
            // Include ignored complex conjugate
            if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
            {
               // Energy in mean
               this->mEnergy(0) += spectrum(start, 0);

            } else if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
            {
               this->mEnergy(1) += 2.0*spectrum(start,0);
            }
            start += 1;
            n -= 1;
         }
         this->mEnergy(1) += 2.0*spectrum.col(0).segment(start, n).sum();
         start += n;
      }

      // Normalize by the Cartesian volume
      this->mEnergy.array() /= this->mVolume;
   }

   void Cartesian1DScalarEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mEnergy.data(), this->mEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mEnergy.sum() << "\t" << this->mEnergy(0) << "\t" << this->mEnergy(1) << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mEnergy.sum()))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Scalar energy is NaN!");
      }
   }

}
}
}
