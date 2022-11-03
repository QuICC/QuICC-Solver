/**
 * @file Cartesian1DTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (toroidal/poloidal formulation)
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
#include "QuICC/Io/Variable/Cartesian1DTorPolEnergyWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyD1.hpp"
#include "QuICC/Io/Variable/EnergyTags.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian1DTorPolEnergyWriter::Cartesian1DTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Energy::BASENAME, Tags::Energy::EXTENSION, prefix + Tags::Energy::HEADER, type, Tags::Energy::VERSION, Dimensions::Space::SPECTRAL), mXEnergy(-1.0), mYEnergy(-1.0), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   Cartesian1DTorPolEnergyWriter::~Cartesian1DTorPolEnergyWriter()
   {
   }

   void Cartesian1DTorPolEnergyWriter::init()
   {
      // Normalize by Cartesian volume V = (2*pi*Box1D/k1D)*(2*pi*Box2D/k2D)*(zo-zi) but FFT already includes 1/(2*pi)
      MHDFloat zi = this->mPhysical.find(NonDimensional::Lower1d::id())->second->value();
      MHDFloat zo = this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();
      this->mVolume = (zo-zi)/(this->res().sim().boxScale(Dimensions::Simulation::SIM2D)*this->res().sim().boxScale(Dimensions::Simulation::SIM3D));

      IVariableAsciiWriter::init();
   }

   void Cartesian1DTorPolEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);

      // Initialize the energy
      this->mXEnergy = 0.0;
      this->mYEnergy = 0.0;
      this->mTorEnergy = 0.0;
      this->mPolEnergy = 0.0;

      // Dealias toroidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute energy reduction
      Matrix spectrum(rInVarTor.data().cols(), 1);
      coord.transform1D().reduce(spectrum, rInVarTor.data(), Transform::Reductor::Energy::id());

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      int start = 0;
      int nK = this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         MHDFloat k_ = static_cast<MHDFloat>(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k));
         // Convert to double Fourier mode
         if(k_ >= nK/2 + (nK % 2))
         {
            k_ = k_ - nK;
         }
         // Include boxscale
         k_ *= this->res().sim().boxScale(Dimensions::Simulation::SIM3D);


         for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            MHDFloat j_ = static_cast<MHDFloat>(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k));
            // Include boxscale
            j_ *= this->res().sim().boxScale(Dimensions::Simulation::SIM2D);
            MHDFloat factor = (k_*k_ + j_*j_);

            if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0)
            {
               // Include ignored complex conjugate
               if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
               {
                  this->mXEnergy += spectrum(start, 0);
               } else if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  this->mTorEnergy += 2.0*factor*spectrum(start,0);
               }
            } else
            {
               this->mTorEnergy += 2.0*factor*spectrum(start,0);
            }
            start += 1;
         }
      }

      // Dealias poloidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute energy reduction
      spectrum.resize(rInVarPol.data().cols(), 1);
      spectrum.setZero();
      coord.transform1D().reduce(spectrum, rInVarPol.data(), Transform::Reductor::Energy::id());

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPol);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      start = 0;
      for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         MHDFloat k_ = static_cast<MHDFloat>(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k));
         // Convert to double Fourier mode
         if(k_ >= nK/2 + (nK % 2))
         {
            k_ = k_ - nK;
         }
         // Include boxscale
         k_ *= this->res().sim().boxScale(Dimensions::Simulation::SIM3D);

         for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            MHDFloat j_ = static_cast<MHDFloat>(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k));
            // Include boxscale
            j_ *= this->res().sim().boxScale(Dimensions::Simulation::SIM2D);
            MHDFloat factor = std::pow((k_*k_ + j_*j_),2);

            if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0)
            {
               // Include ignored complex conjugate
               if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
               {
                  this->mYEnergy += spectrum(start, 0);
               } else if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  this->mPolEnergy += 2.0*factor*spectrum(start, 0);
               }
            } else
            {
               this->mPolEnergy += 2.0*factor*spectrum(start, 0);
            }
            start += 1;
         }
      }

      // Dealias poloidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolDz = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute energy reduction
      spectrum.resize(rInVarPolDz.data().cols(), 1);
      spectrum.setZero();
      coord.transform1D().reduce(spectrum, rInVarPolDz.data(), Transform::Reductor::EnergyD1::id());

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolDz);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      start = 0;
      for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         MHDFloat k_ = static_cast<MHDFloat>(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k));
         // Convert to double Fourier mode
         if(k_ >= nK/2 + (nK % 2))
         {
            k_ = k_ - nK;
         }
         // Include boxscale
         k_ *= this->res().sim().boxScale(Dimensions::Simulation::SIM3D);

         for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            MHDFloat j_ = static_cast<MHDFloat>(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k));
            // Include boxscale
            j_ *= this->res().sim().boxScale(Dimensions::Simulation::SIM2D);
            MHDFloat factor = (k_*k_ + j_*j_);

            if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0)
            {
               // Include ignored complex conjugate
               if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
               {
                  // Do nothing
               } else if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  this->mPolEnergy += 2.0*factor*spectrum(start, 0);
               }
            } else
            {
               this->mPolEnergy += 2.0*factor*spectrum(start, 0);
            }
         }
      }

      // Normalize by the Cartesian volume
      this->mXEnergy /= 2.0*this->mVolume;
      this->mYEnergy /= 2.0*this->mVolume;
      this->mTorEnergy /= 2.0*this->mVolume;
      this->mPolEnergy /= 2.0*this->mVolume;
   }

   void Cartesian1DTorPolEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(4);

         energy(0) = this->mXEnergy;
         energy(1) = this->mYEnergy;
         energy(2) = this->mTorEnergy;
         energy(3) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), energy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mXEnergy = energy(0);
         this->mYEnergy = energy(1);
         this->mTorEnergy = energy(2);
         this->mPolEnergy = energy(3);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mXEnergy + this->mYEnergy + this->mTorEnergy + this->mPolEnergy << "\t" << this->mXEnergy << "\t" << this->mYEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mXEnergy) || std::isnan(this->mYEnergy) || std::isnan(this->mTorEnergy) || std::isnan(this->mPolEnergy))
      {
         QuICCEnv().abort("Kinetic energy is NaN!");
      }
   }

}
}
}
