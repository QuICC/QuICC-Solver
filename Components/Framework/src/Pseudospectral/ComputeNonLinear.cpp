/**
 * @file Coordinator.cpp
 * @brief Source of the high level pseudospectral coordinator
 */

// System includes
//
#include <algorithm>
#include <stdexcept>
#include <type_traits>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Pseudospectral/Utils.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"
#include "View/View.hpp"
#include "View/ViewUtils.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "Profiler/Interface.hpp"


namespace QuICC {

namespace Pseudospectral {


   void Coordinator::computeNonlinear(const int it)
   {
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::computeNonlinear");

      if (mJitter.get() == nullptr)
      {
         // Use old backward tree + non linear terms + forward tree

         Profiler::RegionStart<2>("Pseudospectral::Coordinator::nlOld");
         // Compute backward transform
         this->updatePhysical(it);

         // compute nonlinear interaction and forward transform
         this->updateSpectral(it);
         Profiler::RegionStop<2>("Pseudospectral::Coordinator::nlOld");

      }
      else
      {
         // Use new graph

         // Copy to view
         // Temperature
         const auto& jwRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA1D);
         auto& temp = mScalarVariables[PhysicalNames::Temperature::id()];
         auto tempVarv = mId2View[PhysicalNames::Temperature::id()];
         details::copyScalar2View(tempVarv, temp, jwRes);

         // Spectral Velocity
         auto& vecVel = mVectorVariables[PhysicalNames::Velocity::id()];
         std::size_t hVelTor = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR);
         std::size_t hVelPol = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL);
         auto& TorVelVarv = mId2View[hVelTor];
         auto& PolVelVarv = mId2View[hVelPol];
         details::copyVector2View(TorVelVarv, PolVelVarv, vecVel, jwRes);

         // Magnetic field
         Graph::varData_t TorMagVarv;
         Graph::varData_t PolMagVarv;
         if (mIsMag)
         {
            auto& vecMag = mVectorVariables[PhysicalNames::Magnetic::id()];
            std::size_t hMagTor = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::TOR);
            std::size_t hMagPol = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::POL);
            TorMagVarv = mId2View[hMagTor];
            PolMagVarv = mId2View[hMagPol];
            details::copyVector2View(TorMagVarv, PolMagVarv, vecMag, jwRes);
         }

         // Physical velocity
         std::size_t hVelR = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Physical::R);
         std::size_t hVelTheta = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Physical::THETA);
         std::size_t hVelPhi = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Physical::PHI);

         auto& UrVarv = mId2View[hVelR];
         auto& UthetaVarv = mId2View[hVelTheta];
         auto& UphiVarv = mId2View[hVelPhi];
         // Physical magnetic field
         std::size_t hMagR = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Physical::R);
         std::size_t hMagTheta = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Physical::THETA);
         std::size_t hMagPhi = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Physical::PHI);
         auto& BrVarv = mId2View[hMagR];
         auto& BthetaVarv = mId2View[hMagTheta];
         auto& BphiVarv = mId2View[hMagPhi];

         // #ifndef NDEBUG
         // Profiler::RegionStart<2>("Pseudospectral::Coordinator::nlOld");
         // // Compute backward transform
         // this->updatePhysical(it);

         // // compute nonlinear interaction and forward transform
         // this->updateSpectral(it);
         // Profiler::RegionStop<2>("Pseudospectral::Coordinator::nlOld");
         // #endif

         Profiler::RegionStart<2>("Pseudospectral::Coordinator::nlNew");
         // Call graph

         // Physical space has always the same layout
         using namespace QuICC::Graph;
         auto& Urv = std::get<R_DCCSC3D_t>(UrVarv);
         auto& Uthetav = std::get<R_DCCSC3D_t>(UthetaVarv);
         auto& Uphiv = std::get<R_DCCSC3D_t>(UphiVarv);

         // Magnetic field
         R_DCCSC3D_t Brv;
         R_DCCSC3D_t Bthetav;
         R_DCCSC3D_t Bphiv;
         if (mIsMag)
         {
            Brv = std::get<R_DCCSC3D_t>(BrVarv);
            Bthetav = std::get<R_DCCSC3D_t>(BthetaVarv);
            Bphiv = std::get<R_DCCSC3D_t>(BphiVarv);
         }

         // manually unpacking, the compilers can't handle so many
         // templated variants with a reasonable amount of time/memory
         /// \todo hide ugly switch
         if (std::holds_alternative<C_DCCSC3D_t>(tempVarv))
         {
            // cpu
            auto& Tv = std::get<C_DCCSC3D_t>(tempVarv);
            auto& TorVelv = std::get<C_DCCSC3D_t>(TorVelVarv);
            auto& PolVelv= std::get<C_DCCSC3D_t>(PolVelVarv);


            if (!mIsMag)
            {
               mJitter->apply(Tv, TorVelv, PolVelv,
                     Urv, Uthetav, Uphiv,
                     Tv, TorVelv, PolVelv);
            }
            else
            {
               auto& TorMagv = std::get<C_DCCSC3D_t>(TorMagVarv);
               auto& PolMagv= std::get<C_DCCSC3D_t>(PolMagVarv);
               mJitter->apply(Tv, TorVelv, PolVelv, TorMagv, PolMagv,
                     Urv, Uthetav, Uphiv, Brv, Bthetav, Bphiv,
                     Tv, TorVelv, PolVelv, TorMagv, PolMagv);
            }
         }
         else if (std::holds_alternative<C_DCCSC3DJIK_t>(tempVarv))
         {
            // gpu layout
            auto& Tv = std::get<C_DCCSC3DJIK_t>(tempVarv);
            auto& TorVelv = std::get<C_DCCSC3DJIK_t>(TorVelVarv);
            auto& PolVelv= std::get<C_DCCSC3DJIK_t>(PolVelVarv);

            if (!mIsMag)
            {
               mJitter->apply(Tv, TorVelv, PolVelv,
                     Urv, Uthetav, Uphiv,
                     Tv, TorVelv, PolVelv);
            }
            else
            {
               auto& TorMagv = std::get<C_DCCSC3DJIK_t>(TorMagVarv);
               auto& PolMagv= std::get<C_DCCSC3DJIK_t>(PolMagVarv);
               mJitter->apply(Tv, TorVelv, PolVelv, TorMagv, PolMagv,
                     Urv, Uthetav, Uphiv, Brv, Bthetav, Bphiv,
                     Tv, TorVelv, PolVelv, TorMagv, PolMagv);
            }
         }
         else
         {
            throw std::logic_error("unsupported type");
         }

         Profiler::RegionStop<2>("Pseudospectral::Coordinator::nlNew");

         // Copy back spectral coeff
         details::copyView2Scalar(temp, tempVarv, jwRes);
         details::copyView2Vector(vecVel, TorVelVarv, PolVelVarv, jwRes);

         // Copy back physical vel for cfl computation
         /// \todo move cfl computation in the graph and remove these copies
         const auto& ftRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA3D);
         details::copyView2Vector(vecVel, UrVarv, UthetaVarv, UphiVarv, ftRes);

         if (mIsMag)
         {
            auto& vecMag = mVectorVariables[PhysicalNames::Magnetic::id()];
            // Copy back spectral coeff
            details::copyView2Vector(vecMag, TorMagVarv, PolMagVarv, jwRes);
            // Copy back physical magnetic field for cfl computation
            details::copyView2Vector(vecMag, BrVarv, BthetaVarv, BphiVarv, ftRes);
         }
      }
   }

} // Pseudospectral
} // QuICC
