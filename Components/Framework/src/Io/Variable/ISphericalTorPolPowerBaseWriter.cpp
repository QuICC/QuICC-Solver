/**
 * @file ISphericalTorPolPowerBaseWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics power spectrum calculation for toroidal/poloidal field in a spherical geometry
 */

// Configuration includes
//

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/ISphericalTorPolPowerBaseWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Transform/Reductor/Power.hpp"
#include "QuICC/Transform/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Reductor/PowerD1R1.hpp"
#include "QuICC/Io/Variable/Tags/Power.hpp"

namespace QuICC {

namespace Io {

namespace Variable {
   ISphericalTorPolPowerBaseWriter::ISphericalTorPolPowerBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode)
      : IVariableAsciiWriter(name, ext, header, type, version, id, mode), mHasMOrdering(false), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false)
   {
   }

   ISphericalTorPolPowerBaseWriter::~ISphericalTorPolPowerBaseWriter()
   {
   }

   void ISphericalTorPolPowerBaseWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ISphericalTorPolPowerBaseWriter::prepareInput(const FieldComponents::Spectral::Id sId, Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      auto&& field = vRange.first->second;
      assert(std::visit([&](auto&& p)->bool{return (p->dom(0).res().sim().ss().spectral().ONE() == FieldComponents::Spectral::TOR);}, field));
      assert(std::visit([&](auto&& p)->bool{return (p->dom(0).res().sim().ss().spectral().TWO() == FieldComponents::Spectral::POL);}, field));

      constexpr auto TId = Dimensions::Transform::TRA1D;
      const int packs = 1;
      coord.communicator().converter<TId>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<TId>().prepareBackwardReceive();

      // Dealias toroidal variable data
      std::visit(
            [&](auto&& p)
            {
               coord.communicator().transferForward(Dimensions::Transform::SPECTRAL, p->rDom(0).rTotal().rComp(sId), false);
            },
            field);

      coord.communicator().converter<TId>().initiateForwardSend();
   }

   void ISphericalTorPolPowerBaseWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      constexpr auto TId = Dimensions::Transform::TRA1D;
      Matrix spectrum;

      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::TOR, coord);

      // Recover dealiased BWD data
      auto pInVarTor = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarTor);

      const auto& tRes = *this->res().cpu()->dim(TId);

      // Size of spectrum
      auto size = std::visit([](auto&& p)->std::pair<int,int>{return std::make_pair(0,p->data().cols());}, pInVarTor);
      int l_ = 0;
      if(this->mHasMOrdering)
      {
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               l_ = std::min(l_, tRes.idx<Dimensions::Data::DAT2D>(j, k));
            }
         }
      }
      else
      {
         l_ = tRes.idx<Dimensions::Data::DAT3D>(0);
      }
      size.first = this->res().counter().dimensions(Dimensions::Space::SPECTRAL, l_)(0);

      // Compute power reduction
      spectrum.resize(size.first, size.second);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::PowerR2::id());
            },
            pInVarTor);

      this->resetPower();

      MHDFloat lfactor = 0.0;
      MHDFloat factor = 1.0;
      int idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
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
               lfactor = l_*(l_+1.0);

               for(int i = 0; i < this->res().counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, l_); i++)
               {
                  this->storeTPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
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
            lfactor = l_*(l_+1.0);
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
                  this->storeTPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarTor);

      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::POL, coord);

      // Recover dealiased BWD data
      auto pInVarPolQ = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarPolQ);

      // Compute power reduction
      spectrum.setZero();
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::Power::id());
            },
            pInVarPolQ);

      // Compute power in Q component of QST decomposition
      idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
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
               lfactor = std::pow(l_*(l_+1.0),2);

               for(int i = 0; i < this->res().counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, l_); i++)
               {
                  this->storeQPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
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
            lfactor = std::pow(l_*(l_+1.0),2);
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
                  this->storeQPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarPolQ);

      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::POL, coord);

      // Recover dealiased BWD data
      auto pInVarPolS = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarPolS);

      // Compute power reduction
      spectrum.setZero();
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::PowerD1R1::id());
            },
            pInVarPolS);

      // Compute power in S component of QST decomposition
      idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
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
               lfactor = l_*(l_+1.0);

               for(int i = 0; i < this->res().counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, l_); i++)
               {
                  this->storeSPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
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
            lfactor = l_*(l_+1.0);
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
                  this->storeSPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarPolS);
   }

}
}
}
