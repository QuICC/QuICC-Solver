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

   void ISphericalScalarPowerBaseWriter::prepareInput(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      auto&& field = sRange.first->second;

      constexpr auto TId = Dimensions::Transform::TRA1D;
      const int packs = 1;
      coord.communicator().converter<TId>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<TId>().prepareBackwardReceive();


      // Dealias variable data
      std::visit(
            [&](auto&& p)
            {
               coord.communicator().transferForward(Dimensions::Transform::SPECTRAL, p->rDom(0).rTotal(), false);
            },
            field);

      coord.communicator().converter<TId>().initiateForwardSend();
   }

   void ISphericalScalarPowerBaseWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      constexpr auto TId = Dimensions::Transform::TRA1D;

      // Prepare spectral data for transform
      this->prepareInput(coord);

      // Recover dealiased BWD data
      auto pInVar = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVar);

      const auto& tRes = *this->res().cpu()->dim(TId);

      // Size of spectrum
      auto size = std::visit([](auto&& p)->std::pair<int,int>{return std::make_pair(0,p->data().cols());}, pInVar);
      int nN = 0;
      if(this->mHasMOrdering)
      {
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               std::max(nN, tRes.dim<Dimensions::Data::DATB1D>(j,k));
            }
         }
      }
      else
      {
         nN = tRes.dim<Dimensions::Data::DATB1D>();
      }
      size.first = nN;

      // Compute power reduction
      Matrix spectrum(size.first, size.second);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::PowerR2::id());
            },
            pInVar);

      this->resetPower();

      MHDFloat factor = 1.0;
      int idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
            // m = 0, no factor of two
            if(m_ == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);

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
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);

            // m = 0, no factor of two
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
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
      coord.communicator().storage<TId>().freeBwd(pInVar);
   }

}
}
}
