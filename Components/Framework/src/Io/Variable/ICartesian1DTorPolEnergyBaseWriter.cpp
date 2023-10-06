/**
 * @file ICartesian1DTorPolEnergyBaseWriter.cpp
 * @brief Source of the implementation of the ASCII Chebyshev for toroidal/poloidal field in a plane layer
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Io/Variable/ICartesian1DTorPolEnergyBaseWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyD1.hpp"
#include "QuICC/Io/Variable/Tags/Energy.hpp"

namespace QuICC {

namespace Io {

namespace Variable {
   ICartesian1DTorPolEnergyBaseWriter::ICartesian1DTorPolEnergyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode)
      : IVariableAsciiWriter(name, ext, header, type, version, id, mode), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false)
   {
   }

   void ICartesian1DTorPolEnergyBaseWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ICartesian1DTorPolEnergyBaseWriter::prepareInput(const FieldComponents::Spectral::Id sId, Transform::TransformCoordinatorType& coord)
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

   void ICartesian1DTorPolEnergyBaseWriter::compute(Transform::TransformCoordinatorType& coord)
   {
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
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::Energy::id());
            },
            pInVarTor);

      this->resetEnergy();

      const auto& tRes = *this->res().cpu()->dim(TId);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      int idx = 0;
      const auto box2D = this->res().sim().boxScale(Dimensions::Simulation::SIM2D);
      const auto box3D = this->res().sim().boxScale(Dimensions::Simulation::SIM3D);
      const int nK = this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         const auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         auto ky = static_cast<MHDFloat>(k_);
         // Convert to double Fourier mode
         if(ky >= nK/2 + (nK % 2))
         {
            ky = ky - nK;
         }
         // Include boxscale
         ky *= box3D;

         for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            const auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
            auto kx = static_cast<MHDFloat>(j_);
            // Include boxscale
            kx *= box2D;
            MHDFloat factor = (ky*ky + kx*kx);

            if(j_ == 0)
            {
               // Include ignored complex conjugate
               if(k_ == 0)
               {
                  factor = 1.0;
               } 
               else if(k_ <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  factor = 2.0;
               }
               else
               {
                  factor = 0.0;
               }
            }
            else
            {
               factor = 2.0;
            }
            this->storeTEnergy(j_, k_, factor*spectrum(idx, 0));
            idx += 1;
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
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::Energy::id());
            },
            pInVarPolQ);

      // Compute energy in Q component of QST decomposition
      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      idx = 0;
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         const auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         auto ky = static_cast<MHDFloat>(k_);
         // Convert to double Fourier mode
         if(ky >= nK/2 + (nK % 2))
         {
            ky = ky - nK;
         }
         // Include boxscale
         ky *= box3D;

         for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            const auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
            auto kx = static_cast<MHDFloat>(j_);
            // Include boxscale
            kx *= box2D;
            MHDFloat factor = std::pow((ky*ky + kx*kx),2);

            if(j_ == 0)
            {
               // Include ignored complex conjugate
               if(k_ == 0)
               {
                  factor = 1.0;
               }
               else if(k_ <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  factor = 2.0;
               }
            } else
            {
               factor = 2.0;
            }
            this->storeQEnergy(j_, k_, factor*spectrum(idx, 0));
            idx += 1;
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
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::EnergyD1::id());
            },
            pInVarPolS);

      // Compute energy in S component of QST decomposition
      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      idx = 0;
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         const auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         auto ky = static_cast<MHDFloat>(k_);
         // Convert to double Fourier mode
         if(ky >= nK/2 + (nK % 2))
         {
            ky = ky - nK;
         }
         // Include boxscale
         ky *= box3D;

         for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            const auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
            auto kx = static_cast<MHDFloat>(j_);
            // Include boxscale
            kx *= box2D;
            MHDFloat factor = (ky*ky + kx*kx);

            if(j_ == 0)
            {
               // Include ignored complex conjugate
               if(k_ == 0)
               {
                  // Do nothing
                  factor = 0.0;
               }
               else if(k_ <= this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  factor = 2.0;
               }
            }
            else
            {
               factor = 2.0;
            }

            this->storeSEnergy(j_, k_, factor*spectrum(idx, 0));
            idx += 1;
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarPolS);
   }

} // Variable
} // Io
} // QuICC
