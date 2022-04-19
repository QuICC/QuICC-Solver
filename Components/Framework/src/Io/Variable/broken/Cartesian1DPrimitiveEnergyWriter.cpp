/**
 * @file Cartesian1DPrimitiveEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (primitive formulation)
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
#include "QuICC/Io/Variable/Cartesian1DPrimitiveEnergyWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Io/Variable/EnergyTags.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian1DPrimitiveEnergyWriter::Cartesian1DPrimitiveEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Energy::BASENAME, Tags::Energy::EXTENSION, prefix + Tags::Energy::HEADER, type, Tags::Energy::VERSION, Dimensions::Space::SPECTRAL), mXEnergy(-1.0), mYEnergy(-1.0), mZEnergy(-1.0)
   {
   }

   Cartesian1DPrimitiveEnergyWriter::~Cartesian1DPrimitiveEnergyWriter()
   {
   }

   void Cartesian1DPrimitiveEnergyWriter::init()
   {
      // Normalize by Cartesian volume V = (2*pi*Box1D/k1D)*(2*pi*Box2D/k2D)*(zo-zi) but FFT already includes 1/(2*pi)
      MHDFloat zi = this->mPhysical.find(NonDimensional::Lower1D::id())->second->value();
      MHDFloat zo = this->mPhysical.find(NonDimensional::Upper1D::id())->second->value();
      this->mVolume = (zo-zi)/(this->res().sim().boxScale(Dimensions::Simulation::SIM2D)*this->res().sim().boxScale(Dimensions::Simulation::SIM3D));

      IVariableAsciiWriter::init();
   }

   void Cartesian1DPrimitiveEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);

      // Initialize the energy
      this->mXEnergy = 0.0;
      this->mYEnergy = 0.0;
      this->mZEnergy = 0.0;

      std::vector<FieldComponents::Spectral::Id> comps;
      comps.push_back(FieldComponents::Spectral::X);
      comps.push_back(FieldComponents::Spectral::Y);
      comps.push_back(FieldComponents::Spectral::Z);

      for(auto cIt = comps.begin(); cIt != comps.end(); ++cIt)
      {
         // Dealias toroidal variable data
         coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(*cIt));

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Compute energy reduction
         Matrix spectrum(rInVar.data().cols(), 1);
         coord.transform1D().reduce(spectrum, rInVar.data(), Transform::Reductor::Energy::id());

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

         MHDFloat *pEnergy;
         if(*cIt == FieldComponents::Spectral::X)
         {
            pEnergy = &this->mXEnergy;
         } else if(*cIt == FieldComponents::Spectral::Y)
         {
            pEnergy = &this->mYEnergy;
         } else
         {
            pEnergy = &this->mZEnergy;
         }

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
                  *pEnergy += spectrum(start, 0);
               } else if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  *pEnergy += 2.0*spectrum(start, 0);
               }
               start += 1;
               n -= 1;
            }
            *pEnergy += 2.0*spectrum.col(0).segment(start, n).sum();
            start += n;
         }

         // Normalize by the Cartesian volume
         *pEnergy /= this->mVolume;
      }
   }

   void Cartesian1DPrimitiveEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(3);

         energy(0) = this->mXEnergy;
         energy(1) = this->mYEnergy;
         energy(2) = this->mZEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mXEnergy = energy(0);
         this->mYEnergy = energy(1);
         this->mZEnergy = energy(2);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mXEnergy + this->mYEnergy + this->mZEnergy << "\t" << this->mXEnergy << "\t" << this->mYEnergy << "\t" << this->mZEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mXEnergy) || std::isnan(this->mYEnergy) || std::isnan(this->mZEnergy))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Kinetic energy is NaN!");
      }
   }

}
}
}
