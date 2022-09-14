/**
 * @file ISphericalScalarPowerBaseWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics power calculation for scalar field in a spherical geometry
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
#include "QuICC/Io/Variable/ISphericalScalarPowerBaseWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Transform/Reductor/PowerR2.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/Power.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalScalarPowerBaseWriter::ISphericalScalarPowerBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode)
      : IVariableAsciiWriter(name, ext ,header, type, version, id, mode), mHasMOrdering(false), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false)
   {
   }

   ISphericalScalarPowerBaseWriter::~ISphericalScalarPowerBaseWriter()
   {
   }

   void ISphericalScalarPowerBaseWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ISphericalScalarPowerBaseWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      // Dealias variable data
      std::visit([&](auto&& p){coord.communicator().dealiasSpectral(p->rDom(0).rTotal());}, sRange.first->second);

      // Recover dealiased BWD data
      auto pInVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd(pInVar);

      // Size of spectrum
      auto size = std::visit([](auto&& p)->std::pair<int,int>{return std::make_pair(0,p->data().cols());}, pInVar);
      int l_;
      if(this->mHasMOrdering)
      {
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               l_ = std::min(l_, this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k));
            }
         }
      }
      else
      {
         l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
      }
      size.first = this->res().counter().dimensions(Dimensions::Space::SPECTRAL, l_)(0);

      // Compute power reduction
      Matrix spectrum(size.first, size.second);
      std::visit([&](auto&& p){coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::PowerR2::id());}, pInVar);

      this->resetPower();

      MHDFloat factor = 1.0;
      int idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            // m = 0, no factor of two
            if(m_ == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

               for(int i = 0; i < this->res().counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, l_); i++)
               {
                  this->storePower(i, l_, m_, factor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

            // m = 0, no factor of two
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               // m = 0, no factor of two
               if(m_ == 0)
               {
                  factor = 1.0;
               } else
               {
                  factor = 2.0;
               }

               for(int i = 0; i < this->res().counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, l_); i++)
               {
                  this->storePower(i, l_, m_, factor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(pInVar);
   }

}
}
}
