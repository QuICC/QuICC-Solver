/**
 * @file ICartesian1DScalarEnergyBaseWriter.cpp
 * @brief Source of the implementation of the ASCII Chebyshev energy calculation for scalar field in a plane layer
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/ICartesian1DScalarEnergyBaseWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/Energy.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ICartesian1DScalarEnergyBaseWriter::ICartesian1DScalarEnergyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode)
      : IVariableAsciiWriter(name, ext ,header, type, version, id, mode), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false)
   {
   }

   void ICartesian1DScalarEnergyBaseWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ICartesian1DScalarEnergyBaseWriter::prepareInput(Transform::TransformCoordinatorType& coord)
   {
      // Get field data
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

   void ICartesian1DScalarEnergyBaseWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      DebuggerMacro_msg("ICartesian1DScalarEnergyBaseWriter::compute" ,4);

      constexpr auto TId = Dimensions::Transform::TRA1D;

      // Prepare spectral data for transform
      this->prepareInput(coord);

      // Recover dealiased BWD data
      auto pInVar = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVar);

      // Compute energy reduction
      Matrix spectrum(std::visit([](auto&& p)->int{return p->data().cols();}, pInVar), 1);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::Energy::id());
            },
            pInVar);

      this->resetEnergy();

      const auto& tRes = *this->res().cpu()->dim(TId);

      MHDFloat factor = 1.0;
      int idx = 0;
      // Loop over kx
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         int kx = tRes.idx<Dimensions::Data::DAT3D>(k);
         // Loop over ky
         for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
         {
            int ky = tRes.idx<Dimensions::Data::DAT2D>(j,k);
            if(ky == 0)
            {
               if(kx == 0)
               {
                  factor = 1.0;
               }
               else if(kx < this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  factor = 2.0;
               }
               else
               {
                  factor = 2.0;
               }
            }

            this->storeEnergy(kx, ky, factor*spectrum(idx, 0));
            idx += 1;
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVar);
   }

} // Variable
} // Io
} // QuICC
