/**
 * @file CartesianExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cartesian geometries
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/CartesianExactVectorState.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/NonDimensional/Lower2D.hpp"
#include "QuICC/NonDimensional/Upper2D.hpp"
#include "QuICC/NonDimensional/Lower3D.hpp"
#include "QuICC/NonDimensional/Upper3D.hpp"

#include <iostream>
namespace QuICC {

namespace Equations {

   CartesianExactVectorState::CartesianExactVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IVectorEquation(spEqParams,spScheme)
   {
   }

   CartesianExactVectorState::~CartesianExactVectorState()
   {
   }

   void CartesianExactVectorState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CartesianExactVectorState::setStateType(const CartesianExactStateIds::Id id)
   {
      if(this->ss().physical().ONE() != FieldComponents::Physical::NOTUSED)
      {
         this->mTypeId.insert(std::make_pair(this->ss().physical().ONE(), id));
      }

      if(this->ss().physical().TWO() != FieldComponents::Physical::NOTUSED)
      {
         this->mTypeId.insert(std::make_pair(this->ss().physical().TWO(), id));
      }

      if(this->ss().physical().THREE() != FieldComponents::Physical::NOTUSED)
      {
         this->mTypeId.insert(std::make_pair(this->ss().physical().THREE(), id));
      }
   }

   void CartesianExactVectorState::setStateType(const FieldComponents::Physical::Id compId, const CartesianExactStateIds::Id id)
   {
      this->mTypeId.insert(std::make_pair(compId, id));
   }

   void CartesianExactVectorState::setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2)
   {
      Array tmp(2);
      tmp(0) = a1;
      tmp(1) = a2;
      this->mModeA.insert(std::make_pair(compId, tmp));

      tmp(0) = k1;
      tmp(1) = k2;
      this->mModeK.insert(std::make_pair(compId, tmp));
   }

   void CartesianExactVectorState::setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      Array tmp(3);
      tmp(0) = a1;
      tmp(1) = a2;
      tmp(2) = a3;
      this->mModeA.insert(std::make_pair(compId, tmp));

      tmp(0) = k1;
      tmp(1) = k2;
      tmp(2) = k3;
      this->mModeK.insert(std::make_pair(compId, tmp));
   }

   void CartesianExactVectorState::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;
      features.at(CouplingFeature::AllowExplicit) = false;

      if(this->ss().spectral().ONE() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().ONE(), CouplingInformation::TRIVIAL, 0, features);
      }

      if(this->ss().spectral().TWO() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().TWO(), CouplingInformation::TRIVIAL, 0, features);
      }

      if(this->ss().spectral().THREE() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().THREE(), CouplingInformation::TRIVIAL, 0, features);
      }
   }

   void CartesianExactVectorState::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
   {
      CartesianExactStateIds::Id typeId = this->mTypeId.find(compId)->second;

      if(typeId == CartesianExactStateIds::CONSTANT)
      {
         Array modeA = this->mModeA.find(compId)->second;
         rNLComp.rData().setConstant(modeA.prod());

      // Generate constant unit field divergence free state for toroidal/poloidal decomposition test
      } else if(typeId == CartesianExactStateIds::TORPOLCNST)
      {
         #ifdef QUICC_SPATIALDIMENSION_3D
            Array gI, gJ, gK;
            this->buildGrid(gI, gJ, gK);
         #else
            Array gI, gJ;
            this->buildGrid(gI, gJ);
         #endif //QUICC_SPATIALDIMENSION_3D

         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
         for(int iK = 0; iK < nK; ++iK)
         {
            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
            for(int iJ = 0; iJ < nJ; ++iJ)
            {
               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
               for(int iI = 0; iI < nI; ++iI)
               {
                  if(compId == FieldComponents::Physical::X || compId == FieldComponents::Physical::Y)
                  {
                     rNLComp.setPoint(1.0, iI, iJ, iK);
                  }
               }
            }
         }

      // Generate divergence free state for toroidal/poloidal decomposition test
      } else if(typeId == CartesianExactStateIds::TORPOLTFF)
      {
         #ifdef QUICC_SPATIALDIMENSION_3D
            Array gI, gJ, gK;
            this->buildGrid(gI, gJ, gK);
         #else
            Array gI, gJ;
            this->buildGrid(gI, gJ);
         #endif //QUICC_SPATIALDIMENSION_3D

         MHDFloat k_;
         MHDFloat j_;
         MHDFloat i_;

         int nSI = 3;
         int nSJ = 3;
         int nSK = 2;
         MHDFloat bSI = this->res().sim().boxScale(Dimensions::Simulation::SIM2D);
         MHDFloat bSJ = this->res().sim().boxScale(Dimensions::Simulation::SIM3D);

         Array aT(5); aT(0) = -3.0; aT(1) = 2.0; aT(2) = 3.0; aT(3) = -1.0; aT(4) = 5.0;
         Array mT(5); mT(0) = 2.0; mT(1) = 1.0; mT(2) = 2.0; mT(3) = 1.0; mT(4) = 3.0;
         Array aP(5); aP(0) = 6.0; aP(1) = -7.0; aP(2) = 5.0; aP(3) = 2.0; aP(4) = 5.0;
         Array mP(5); mP(0) = -5.0; mP(1) = 4.0; mP(2) = 3.0; mP(3) = -3.0; mP(4) = 1.0;

         // Chebyshev rescaling
         MHDFloat scale = 2.0/(this->eqParams().nd(NonDimensional::Upper1D::id()) - this->eqParams().nd(NonDimensional::Lower1D::id()));

         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
         for(int iK = 0; iK < nK; ++iK)
         {
            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
            for(int iJ = 0; iJ < nJ; ++iJ)
            {
               j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
               for(int iI = 0; iI < nI; ++iI)
               {
                  i_ = gI(iI);

                  MHDFloat val = 0.0;

                  if(compId == FieldComponents::Physical::X)
                  {
                     // Toroidal component
                     for(int sI = 0; sI < nSI; sI++)
                     {
                        MHDFloat sI_ = static_cast<MHDFloat>(sI*bSI);
                        MHDFloat valI = sI_*(-std::sin(sI*i_)+std::cos(sI*i_));
                        for(int sJ = 0; sJ < nSJ; sJ++)
                        {
                           MHDFloat sJ_ = static_cast<MHDFloat>(sJ*bSJ);
                           MHDFloat valJ = (std::cos(sJ*j_)+std::sin(sJ*j_));

                           for(int sK = 0; sK < nSK; sK++)
                           {
                              val += CartesianExactStateIds::chebyshev(aT(sK),sK,k_)*valI*valJ;
                           }
                        }
                     }

                     // Poloidal component
                     for(int sI = 0; sI < nSI; sI++)
                     {
                        MHDFloat sI_ = static_cast<MHDFloat>(sI*bSI);
                        MHDFloat valI = (std::cos(sI*i_)+std::sin(sI*i_));

                        for(int sJ = 0; sJ < nSJ; sJ++)
                        {
                           MHDFloat sJ_ = static_cast<MHDFloat>(sJ*bSJ);
                           MHDFloat valJ = scale*sJ_*(-std::sin(sJ*j_)+std::cos(sJ*j_));

                           if(nSK > 1)
                           {
                              val += (aP(1))*valI*valJ;
                              if(nSK > 2)
                              {
                                 val += (4.0*aP(2)*k_)*valI*valJ;
                                 if(nSK > 3)
                                 {
                                    val += aP(3)*(-3.0 + 12.0*k_*k_)*valI*valJ;
                                    if(nSK > 4)
                                    {
                                       val += aP(4)*(-16.0*k_ + 32.0*k_*k_*k_)*valI*valJ;
                                    }
                                 }
                              }
                           }
                        }
                     }

                     // X mean component
                     for(int sK = 0; sK < nSK; sK++)
                     {
                        val += CartesianExactStateIds::chebyshev(mT(sK),sK,k_);
                     }

                  } else if(compId == FieldComponents::Physical::Y)
                  {
                     // Toroidal component
                     for(int sI = 0; sI < nSI; sI++)
                     {
                        MHDFloat sI_ = static_cast<MHDFloat>(sI*bSI);
                        MHDFloat valI = (std::cos(sI*i_)+std::sin(sI*i_));
                        for(int sJ = 0; sJ < nSJ; sJ++)
                        {
                           MHDFloat sJ_ = static_cast<MHDFloat>(sJ*bSJ);
                           MHDFloat valJ = -sJ_*(-std::sin(sJ*j_)+std::cos(sJ*j_));

                           for(int sK = 0; sK < nSK; sK++)
                           {
                              val += CartesianExactStateIds::chebyshev(aT(sK),sK,k_)*valI*valJ;
                           }
                        }
                     }

                     // Poloidal component
                     for(int sI = 0; sI < nSI; sI++)
                     {
                        MHDFloat sI_ = static_cast<MHDFloat>(sI*bSI);
                        MHDFloat valI = sI_*(-std::sin(sI*i_)+std::cos(sI*i_));

                        for(int sJ = 0; sJ < nSJ; sJ++)
                        {
                           MHDFloat sJ_ = static_cast<MHDFloat>(sJ*bSJ);
                           MHDFloat valJ = scale*(std::cos(sJ*j_)+std::sin(sJ*j_));

                           if(nSK > 1)
                           {
                              val += (aP(1))*valI*valJ;
                              if(nSK > 2)
                              {
                                 val += (4.0*aP(2)*k_)*valI*valJ;
                                 if(nSK > 3)
                                 {
                                    val += aP(3)*(-3.0 + 12.0*k_*k_)*valI*valJ;
                                    if(nSK > 4)
                                    {
                                       val += aP(4)*(-16.0*k_ + 32.0*k_*k_*k_)*valI*valJ;
                                    }
                                 }
                              }
                           }
                        }
                     }

                     // Y mean component
                     for(int sK = 0; sK < nSK; sK++)
                     {
                        val += CartesianExactStateIds::chebyshev(mP(sK),sK,k_);
                     }

                  } else if(compId == FieldComponents::Physical::Z)
                  {
                     // Poloidal component
                     for(int sI = 0; sI < nSI; sI++)
                     {
                        MHDFloat sI_ = static_cast<MHDFloat>(sI*bSI);
                        MHDFloat valI = (std::cos(sI*i_)+std::sin(sI*i_));
                        for(int sJ = 0; sJ < nSJ; sJ++)
                        {
                           MHDFloat sJ_ = static_cast<MHDFloat>(sJ*bSJ);
                           MHDFloat valJ = (sI_*sI_ + sJ_*sJ_)*(std::cos(sJ*j_)+std::sin(sJ*j_));

                           for(int sK = 0; sK < nSK; sK++)
                           {
                              val += CartesianExactStateIds::chebyshev(aP(sK),sK,k_)*valI*valJ;
                           }
                        }
                     }
                  }

                  rNLComp.setPoint(val, iI, iJ, iK);
               }
            }
         }

      } else
      {
         Array modeA = this->mModeA.find(compId)->second;
         Array modeK = this->mModeK.find(compId)->second;
         #ifdef QUICC_SPATIALDIMENSION_3D
            Array gI, gJ, gK;
            this->buildGrid(gI, gJ, gK);
         #else
            Array gI, gJ;
            this->buildGrid(gI, gJ);
         #endif //QUICC_SPATIALDIMENSION_3D

         MHDFloat k_;
         MHDFloat j_;
         MHDFloat i_;
         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
         for(int iK = 0; iK < nK; ++iK)
         {
            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
            for(int iJ = 0; iJ < nJ; ++iJ)
            {
               j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
               for(int iI = 0; iI < nI; ++iI)
               {
                  i_ = gI(iI);

                  MHDFloat val = 0.0;
                  if(static_cast<int>(typeId) < 100)
                  {
                     Array grid(3);
                     grid(0) = k_;
                     grid(1) = j_;
                     grid(2) = i_;
                     val = CartesianExactStateIds::exact3D(typeId, modeA, modeK, grid);
                  } else if(static_cast<int>(typeId) >= 100)
                  {
                     Array grid(2);
                     grid(0) = j_;
                     grid(1) = i_;
                     val = CartesianExactStateIds::exact2D(typeId, modeA, modeK, grid);
                  }

                  rNLComp.setPoint(val, iI, iJ, iK);
               }
            }
         }
      }
   }

   MHDVariant CartesianExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      return MHDComplex(0);
   }

   void CartesianExactVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(false, ss.spectral(), ss.physical()));
      req.enableSpectral();
      req.enablePhysical();
   }

   void CartesianExactVectorState::setNLComponents()
   {
      if(this->ss().spectral().ONE() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().ONE(), 0);
      }

      if(this->ss().spectral().TWO() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().TWO(), 0);
      }

      if(this->ss().spectral().THREE() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().THREE(), 0);
      }
   }

   void CartesianExactVectorState::buildGrid(Array& g1D, Array& g2D, Array& g3D) const
   {
      int nK = this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->res().sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      g3D = this->mspMesh->at(2);
      g2D = this->mspMesh->at(1);
      g1D = this->mspMesh->at(0);
   }

}
}
