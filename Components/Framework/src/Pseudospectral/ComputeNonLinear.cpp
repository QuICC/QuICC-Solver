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
// #include "QuICC/Debug/DebuggerMacro.h"
// #include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
// #include "Environment/QuICCEnv.hpp"
// #include "QuICC/QuICCTimer.hpp"
// #include "QuICC/ModelOperator/ExplicitLinear.hpp"
// #include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
// #include "QuICC/ModelOperator/ExplicitNextstep.hpp"
// #include "QuICC/SolveTiming/After.hpp"
// #include "QuICC/SolveTiming/Before.hpp"
// #include "QuICC/SolveTiming/Prognostic.hpp"
// #include "QuICC/Variables/RequirementTools.hpp"
// #include "QuICC/TransformCoordinators/TransformCoordinatorTools.hpp"
// #include "QuICC/Equations/Tools/EquationTools.hpp"
// #include "QuICC/Simulation/SimulationIoTools.hpp"
// #include "QuICC/PseudospectralTag/Diagnostic.hpp"
// #include "QuICC/PseudospectralTag/Prognostic.hpp"
// #include "QuICC/PseudospectralTag/Trivial.hpp"
// #include "QuICC/PseudospectralTag/Uninitialized.hpp"
// #include "QuICC/PseudospectralTag/Wrapper.hpp"
// #include "QuICC/PhysicalNames/Coordinator.hpp"
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
         auto& UrVarv = mId2View[FieldComponents::Physical::R];
         auto& UthetaVarv = mId2View[FieldComponents::Physical::THETA];
         auto& UphiVarv = mId2View[FieldComponents::Physical::PHI];

         #ifndef NDEBUG
         Profiler::RegionStart<2>("Pseudospectral::Coordinator::nlOld");
         // Compute backward transform
         this->updatePhysical(it);

         // compute nonlinear interaction and forward transform
         this->updateSpectral(it);
         Profiler::RegionStop<2>("Pseudospectral::Coordinator::nlOld");
         #endif

         Profiler::RegionStart<2>("Pseudospectral::Coordinator::nlNew");
         // Call graph

         // Physical space has always the same layout
         using namespace QuICC::Graph;
         auto& Urv = std::get<R_DCCSC3D_t>(UrVarv);
         auto& Uthetav = std::get<R_DCCSC3D_t>(UthetaVarv);
         auto& Uphiv = std::get<R_DCCSC3D_t>(UphiVarv);

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
                     Urv, Uthetav, Uphiv,
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
                     Urv, Uthetav, Uphiv,
                     Tv, TorVelv, PolVelv, TorMagv, PolMagv);
            }
         }
         else
         {
            throw std::logic_error("unsopported type");
         }

         Profiler::RegionStop<2>("Pseudospectral::Coordinator::nlNew");

         // Copy back spectral coeff
         details::copyView2Scalar(temp, tempVarv, jwRes);
         details::copyView2Vector(vecVel, TorVelVarv, PolVelVarv, jwRes);

         // Copy back physical vel for cfl computation
         /// \todo move cfl comp in graph and remove
         const auto& ftRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA3D);
         details::copyView2Vector(vecVel, UrVarv, UthetaVarv, UphiVarv, ftRes);
      }
   }

} // Pseudospectral
} // QuICC
