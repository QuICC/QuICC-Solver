/**
 * @file ISphericalTorPolEnergyBaseWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for toroidal/poloidal field in a spherical geometry
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolEnergyBaseWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"
#include "QuICC/Transform/Path/InsulatingTorPol.hpp"
#include "QuICC/Transform/Path/NoSlipTorPol.hpp"
#include "QuICC/Transform/Path/NoPenetrationTorPol.hpp"
#include "QuICC/Transform/Path/StressFreeTorPol.hpp"
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/ValueEnergy.hpp"
#include "QuICC/Transform/Reductor/ValueEnergyR2.hpp"
#include "QuICC/Transform/Reductor/ValueEnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/InsulatingEnergy.hpp"
#include "QuICC/Transform/Reductor/InsulatingEnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/InsulatingEnergyR2.hpp"
#include "QuICC/Transform/Reductor/NoSlipEnergy.hpp"
#include "QuICC/Transform/Reductor/NoSlipEnergyD1R1.hpp"
#include "QuICC/Io/Variable/Tags/Energy.hpp"

namespace QuICC {

namespace Io {

namespace Variable {
   ISphericalTorPolEnergyBaseWriter::ISphericalTorPolEnergyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode)
      : IVariableAsciiWriter(name, ext, header, type, version, id, mode), mHasMOrdering(false), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false)
   {
      // Set default path
      this->setTransformPath(Transform::Path::TorPol::id());
   }

   void ISphericalTorPolEnergyBaseWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ISphericalTorPolEnergyBaseWriter::prepareInput(const FieldComponents::Spectral::Id sId, Transform::TransformCoordinatorType& coord)
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

      // Dealias variable data
      std::visit(
            [&](auto&& p)
            {
               coord.communicator().transferForward(Dimensions::Transform::SPECTRAL, p->rDom(0).rTotal().rComp(sId), false);
            },
            field);

      coord.communicator().converter<TId>().initiateForwardSend();
   }

   void ISphericalTorPolEnergyBaseWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      DebuggerMacro_msg("ISphericalTorPolEnergyBaseWriter::compute" ,4);

      std::size_t torEnergyR2Id;
      std::size_t polEnergyId;
      std::size_t polEnergyD1R1Id;

      if(this->mPathId == Transform::Path::TorPol::id())
      {
         torEnergyR2Id = Transform::Reductor::EnergyR2::id();
         polEnergyId = Transform::Reductor::Energy::id();
         polEnergyD1R1Id = Transform::Reductor::EnergyD1R1::id();
      }
      else if(this->mPathId == Transform::Path::InsulatingTorPol::id())
      {
         torEnergyR2Id = Transform::Reductor::ValueEnergyR2::id();
         polEnergyId = Transform::Reductor::InsulatingEnergy::id();
         polEnergyD1R1Id = Transform::Reductor::InsulatingEnergyD1R1::id();
      }
      else if(this->mPathId == Transform::Path::NoSlipTorPol::id())
      {
         torEnergyR2Id = Transform::Reductor::ValueEnergyR2::id();
         polEnergyId = Transform::Reductor::NoSlipEnergy::id();
         polEnergyD1R1Id = Transform::Reductor::NoSlipEnergyD1R1::id();
      }
      else if(this->mPathId == Transform::Path::NoPenetrationTorPol::id())
      {
         torEnergyR2Id = Transform::Reductor::InsulatingEnergyR2::id();
         polEnergyId = Transform::Reductor::ValueEnergy::id();
         polEnergyD1R1Id = Transform::Reductor::ValueEnergyD1R1::id();
      }
      else
      {
         throw std::logic_error("Unknown energy transform reductor path (" + std::to_string(this->mPathId) + ") requested for Toroidal/Poloidal");
      }

      constexpr auto TId = Dimensions::Transform::TRA1D;
      Matrix spectrum;

      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::TOR, coord);

      // Recover dealiased BWD data
      auto pInVarTor = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarTor);

      // Compute energy reduction
      spectrum.resize(std::visit([](auto&& p)->int{return p->data().cols();}, pInVarTor), 1);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), torEnergyR2Id);
            },
            pInVarTor);

      this->resetEnergy();

      const auto& tRes = *this->res().cpu()->dim(TId);

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

               this->storeTEnergy(l_, m_, factor*lfactor*spectrum(idx, 0));
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

               this->storeTEnergy(l_, m_, factor*lfactor*spectrum(idx, 0));
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

      // Compute energy reduction
      spectrum.resize(std::visit([](auto&& p)->int{return p->data().cols();}, pInVarPolQ), 1);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), polEnergyId);
            },
            pInVarPolQ);

      // Compute energy in Q component of QST decomposition
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

               this->storeQEnergy(l_, m_, factor*lfactor*spectrum(idx,0));
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

               this->storeQEnergy(l_, m_, factor*lfactor*spectrum(idx,0));
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

      // Compute energy reduction
      spectrum.resize(std::visit([](auto&& p)->int{return p->data().cols();}, pInVarPolS), 1);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), polEnergyD1R1Id);
            },
            pInVarPolS);

      // Compute energy in S component of QST decomposition
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

               this->storeSEnergy(l_, m_, factor*lfactor*spectrum(idx,0));
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

               this->storeSEnergy(l_, m_, factor*lfactor*spectrum(idx,0));
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarPolS);
   }

} // Variable
} // Io
} // QuICC
