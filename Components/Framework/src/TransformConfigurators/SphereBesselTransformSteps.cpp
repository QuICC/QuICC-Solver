/**
 * @file SphereBesselTransformSteps.cpp
 * @brief Source of the implementation of the physical <-> spectral transform steps in a whole shell
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/TransformConfigurators/SphereBesselTransformSteps.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Path/ValueScalar.hpp"
#include "QuICC/Transform/Path/InsulatingScalar.hpp"
#include "QuICC/Transform/Path/ValueTorPol.hpp"
#include "QuICC/Transform/Path/InsulatingTorPol.hpp"
#include "QuICC/Transform/Path/NoSlipTorPol.hpp"
#include "QuICC/Transform/Path/StressFreeTorPol.hpp"
#include "QuICC/Transform/Path/ValueScalarNl.hpp"
#include "QuICC/Transform/Path/InsulatingScalarNl.hpp"
#include "QuICC/Transform/Path/ValueCurlNl.hpp"
#include "QuICC/Transform/Path/InsulatingCurlNl.hpp"
#include "QuICC/Transform/Path/StressFreeCurlNl.hpp"
#include "QuICC/Transform/Path/ValueCurlCurlNl.hpp"
#include "QuICC/Transform/Path/InsulatingCurlCurlNl.hpp"
#include "QuICC/Transform/Path/ValueNegCurlCurlNl.hpp"
#include "QuICC/Transform/Path/ValueBc1NegCurlCurlNl.hpp"
#include "QuICC/Transform/Path/InsulatingNegCurlCurlNl.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/ValueP.hpp"
#include "QuICC/Transform/Forward/InsulatingP.hpp"
#include "QuICC/Transform/Forward/StressFreeP.hpp"
#include "QuICC/Transform/Forward/Overlaplh.hpp"
#include "QuICC/Transform/Forward/OverlaplhD1.hpp"
#include "QuICC/Transform/Forward/OverlaplhOversinDphi.hpp"
#include "QuICC/Transform/Forward/ValuePol.hpp"
#include "QuICC/Transform/Forward/InsulatingPol.hpp"
#include "QuICC/Transform/Forward/ValueQ.hpp"
#include "QuICC/Transform/Forward/InsulatingQ.hpp"
#include "QuICC/Transform/Forward/ValueBc1Q.hpp"
#include "QuICC/Transform/Forward/ValueS.hpp"
#include "QuICC/Transform/Forward/InsulatingS.hpp"
#include "QuICC/Transform/Forward/ValueBc1S.hpp"
#include "QuICC/Transform/Forward/ValueT.hpp"
#include "QuICC/Transform/Forward/InsulatingT.hpp"
#include "QuICC/Transform/Forward/StressFreeT.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/ValueP.hpp"
#include "QuICC/Transform/Backward/InsulatingP.hpp"
#include "QuICC/Transform/Backward/StressFreeP.hpp"
#include "QuICC/Transform/Backward/ValueOverr1.hpp"
#include "QuICC/Transform/Backward/InsulatingOverr1.hpp"
#include "QuICC/Transform/Backward/StressFreeOverr1.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/ValueD1.hpp"
#include "QuICC/Transform/Backward/InsulatingD1.hpp"
#include "QuICC/Transform/Backward/ValueOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/InsulatingOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/StressFreeOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/OversinDphi.hpp"
#include "QuICC/Transform/Backward/ValueSlapl.hpp"
#include "QuICC/Transform/Backward/OversinD1Sin.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"

// Value Torpol = Tor Value, Pol Value
// No-slip Torpol = Tor Value, Pol Value
// Stress-free Torpol = SF Value, Pol Value
// Insulating Torpol = Tor Value, Pol Insulating
namespace QuICC {

namespace Transform {

   SphereBesselTransformSteps::SphereBesselTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
      : ITransformSteps(spScheme)
   {
   }

   bool SphereBesselTransformSteps::applicable(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
   {
      bool ret = false;
      if(spScheme->tag() == "JLFl" || spScheme->tag() == "JLFm")
      {
         ret = true;
      }

      return ret;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::forwardScalar(const std::vector<PathId >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;
      auto flag = components.at(0).second;

      std::size_t rId;
      // Value BC
      if(flag == Path::ValueScalar::id())
      {
         rId = Forward::ValueP::id();
      }
      // Insulating BC
      else if(flag == Path::InsulatingScalar::id())
      {
         rId = Forward::InsulatingP::id();
      }
      else
      {
         throw std::logic_error("Unknown backward scalar transform path");
      }

      transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(rId, scalId, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::forwardNLScalar(const std::vector<PathId >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;
      auto flag = components.at(0).second;

      std::size_t rId;
      // Value BC
      if(flag == Path::ValueScalarNl::id())
      {
         rId = Forward::ValueP::id();
      }
      // Insulating BC
      else if(flag == Path::InsulatingScalarNl::id())
      {
         rId = Forward::InsulatingP::id();
      }
      else
      {
         throw std::logic_error("Requested an unknown nonlinear scalar forward transform (ID = " + std::to_string(flag) + ")");
      }

      DebuggerMacro_msg("Using SphereBesselTransformSteps for forwardNLScalar", 1);

      transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(rId, scalId, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::forwardVector(const std::vector<PathId >& components) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         assert(components.size() == 2);
         FieldComponents::Spectral::Id curlId = components.at(0).first;
         auto curlFlag = components.at(0).second;
         FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
         auto curlcurlFlag = components.at(1).second;

         // Value BC for toroidal and Value BC for poloidal
         if((curlFlag == Path::ValueTorPol::id() && curlcurlFlag == curlFlag) || 
            (curlFlag == Path::NoSlipTorPol::id() && curlcurlFlag == curlFlag) || 
               (curlFlag == Path::InsulatingTorPol::id() && curlcurlFlag == curlFlag))
         {
            std::size_t tId;
            std::size_t pId;
            // Value BC
            if(curlFlag == Path::ValueTorPol::id() || curlFlag == Path::NoSlipTorPol::id())
            {
               tId = Forward::ValueP::id();
               pId = Forward::ValuePol::id();
            }
            // Stress-free BC
            else if(curlFlag == Path::StressFreeTorPol::id())
            {
               tId = Forward::StressFreeP::id();
               pId = Forward::ValuePol::id();
            }
            // Insulating BC
            else if(curlFlag == Path::InsulatingTorPol::id())
            {
               tId = Forward::ValueP::id();
               pId = Forward::InsulatingPol::id();
            }
            // Compute Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(tId, curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(tId, curlId, Arithmetics::Sub::id());

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Overlaplh::id());
            transform.back().addEdge(pId, curlcurlId, Arithmetics::Add::id());
         }
         else
         {
            throw std::logic_error("Requested an unknown vector forward transform (ID = " + std::to_string(curlFlag) + ")");
         }
      } else
      {
         throw std::logic_error("SphereBesseTransformSteps for primitive forwardVector not implemented");
      }

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::forwardNLVector(const std::vector<PathId >& components) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         DebuggerMacro_msg("Using SphereBesselTransformSteps for Toroidal/Poloidal forwardNLVector", 1);

         assert(components.size() == 2);
         FieldComponents::Spectral::Id curlId = components.at(0).first;
         auto curlFlag = components.at(0).second;
         FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
         auto curlcurlFlag = components.at(1).second;

         // CurlNl BC
         if(curlFlag == Path::ValueCurlNl::id() || 
               curlFlag == Path::StressFreeCurlNl::id() || 
               curlFlag == Path::InsulatingCurlNl::id())
         {
            std::size_t tId;
            // Value BC
            if(curlFlag == Path::ValueCurlNl::id())
            {
               tId = Forward::ValueT::id();
            }
            // Insulating BC
            else if(curlFlag == Path::InsulatingCurlNl::id())
            {
               tId = Forward::InsulatingT::id();
            }
            // Stress-Free BC
            else if(curlFlag == Path::StressFreeCurlNl::id())
            {
               tId = Forward::StressFreeT::id();
            }
            else
            {
               throw std::logic_error("Unknown path requested");
            }
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(tId, curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(tId, curlId, Arithmetics::Sub::id());
         }
         else
         {
            throw std::logic_error("Requested an unknown curl nonlinear vector forward transform (ID = " + std::to_string(curlFlag) + ")");
         }

         // CurlCurlNl BC 
         if(curlcurlFlag == Path::ValueCurlCurlNl::id() || curlcurlFlag == Path::InsulatingCurlCurlNl::id())
         {
            std::size_t qId;
            std::size_t sId;
            // Value BC 
            if(curlcurlFlag == Path::ValueCurlCurlNl::id())
            {
               qId = Forward::ValueQ::id();
               sId = Forward::ValueS::id();
            }
            // Insulating BC
            else if(curlcurlFlag == Path::InsulatingCurlCurlNl::id())
            {
               qId = Forward::InsulatingQ::id();
               sId = Forward::InsulatingS::id();
            }
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(qId, curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(sId, curlcurlId, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(sId, curlcurlId, Arithmetics::Sub::id());
         }
         // Negative CurlCurlNl
         else if(curlcurlFlag == Path::ValueNegCurlCurlNl::id() || curlcurlFlag == Path::ValueBc1NegCurlCurlNl::id())
         {
            std::size_t qId;
            std::size_t sId;
            // Value BC
            if(curlcurlFlag == Path::ValueNegCurlCurlNl::id())
            {
               qId = Forward::ValueQ::id();
               sId = Forward::ValueS::id();
            }
            // Value BC with 1 additional BC
            else if(curlcurlFlag == Path::ValueBc1NegCurlCurlNl::id())
            {
               qId = Forward::ValueBc1Q::id();
               sId = Forward::ValueBc1S::id();
            }
            // Insulating BC
            else if(curlcurlFlag == Path::InsulatingNegCurlCurlNl::id())
            {
               qId = Forward::InsulatingQ::id();
               sId = Forward::InsulatingS::id();
            }
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(qId, curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(sId, curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(sId, curlcurlId, Arithmetics::Add::id());
         }
         else
         {
            throw std::logic_error("Requested an unknown curlcurl vector forward transform (ID = " + std::to_string(curlcurlFlag) + ")");
         }

      }
      // The following assumes the physical values are obtained from a primitive formulation
      else
      {
         throw std::logic_error("SphereBesseTransformSteps for primitive forwardNLVector not implemented");
      }
      DebuggerMacro_msg("... done", 1);

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardScalar(const PhysPathId& req) const
   {
      std::vector<TransformPath> transform;
      std::size_t rId;

      // Value BC
      if(req.find(FieldComponents::Physical::SCALAR)->second == Path::ValueScalar::id())
      {
         rId = Backward::ValueP::id();
      }
      // Insulating BC
      else if(req.find(FieldComponents::Physical::SCALAR)->second == Path::InsulatingScalar::id())
      {
         rId = Backward::InsulatingP::id();
      }
      else
      {
         throw std::logic_error("Unknown backward scalar transform path");
      }

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(rId);
      transform.back().addEdge(Backward::P::id());
      transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardGradient(const PhysPathId& req) const
   {
      std::vector<TransformPath> transform;
      std::size_t rId;

      // Value BC
      if(req.find(FieldComponents::Physical::R)->second == Path::ValueScalar::id())
      {
         rId = Backward::ValueD1::id();
      }
      // Insulating BC
      else if(req.find(FieldComponents::Physical::R)->second == Path::InsulatingScalar::id())
      {
         rId = Backward::InsulatingD1::id();
      }
      else
      {
         throw std::logic_error("Unknown path transform path");
      }
      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
      transform.back().addEdge(rId);
      transform.back().addEdge(Backward::P::id());
      transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());

      // Value BC
      if(req.find(FieldComponents::Physical::THETA)->second == Path::ValueScalar::id())
      {
         rId = Backward::ValueOverr1::id();
      }
      // Insulating BC
      else if(req.find(FieldComponents::Physical::THETA)->second == Path::InsulatingScalar::id())
      {
         rId = Backward::InsulatingOverr1::id();
      }
      else
      {
         throw std::logic_error("Unknown path transform path");
      }
      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
      transform.back().addEdge(rId);
      transform.back().addEdge(Backward::D1::id());
      transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

      // Value BC
      if(req.find(FieldComponents::Physical::PHI)->second == Path::ValueScalar::id())
      {
         rId = Backward::ValueOverr1::id();
      }
      // Insulating BC
      else if(req.find(FieldComponents::Physical::PHI)->second == Path::InsulatingScalar::id())
      {
         rId = Backward::InsulatingOverr1::id();
      }
      else
      {
         throw std::logic_error("Unknown path transform path");
      }
      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
      transform.back().addEdge(rId);
      transform.back().addEdge(Backward::OversinDphi::id());
      transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardGradient2(const Grad2PathId& req) const
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      throw std::logic_error("Second derivative is not implementated yet!");

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardVector(const PhysPathId& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         std::size_t rTorId;
         std::size_t rPolId;

         auto cR = req.find(FieldComponents::Physical::R)->second;
         if(cR == Path::ValueTorPol::id() || cR == Path::NoSlipTorPol::id() || cR == Path::StressFreeTorPol::id())
         {
            rPolId = Backward::ValueOverr1::id();
         }
         else if(cR == Path::InsulatingTorPol::id())
         {
            rPolId = Backward::InsulatingOverr1::id();
         }
         else
         {
            throw std::logic_error("Unknown path transform path");
         }
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(rPolId);
         transform.back().addEdge(Backward::Laplh::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());

         auto cT = req.find(FieldComponents::Physical::THETA)->second;
         if(cT == Path::ValueTorPol::id() || cT == Path::NoSlipTorPol::id())
         {
            rTorId = Backward::ValueP::id();
            rPolId = Backward::ValueOverr1D1R1::id();
         }
         else if(cT == Path::StressFreeTorPol::id())
         {
            rTorId = Backward::StressFreeP::id();
            rPolId = Backward::ValueOverr1D1R1::id();
         }
         else if(cT == Path::InsulatingTorPol::id())
         {
            rTorId = Backward::InsulatingP::id();
            rPolId = Backward::ValueOverr1D1R1::id();
         }
         else
         {
            throw std::logic_error("Unknown path transform path");
         }
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(rTorId);
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(rPolId);
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

         auto cP = req.find(FieldComponents::Physical::PHI)->second;
         if(cP == Path::ValueTorPol::id() || cP == Path::NoSlipTorPol::id())
         {
            rTorId = Backward::ValueP::id();
            rPolId = Backward::ValueOverr1D1R1::id();
         }
         else if(cP == Path::StressFreeTorPol::id())
         {
            rTorId = Backward::StressFreeP::id();
            rPolId = Backward::ValueOverr1D1R1::id();
         }
         else if(cP == Path::InsulatingTorPol::id())
         {
            rTorId = Backward::InsulatingP::id();
            rPolId = Backward::ValueOverr1D1R1::id();
         }
         else
         {
            throw std::logic_error("Unknown path transform path");
         }
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(rTorId);
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());
         
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(rPolId);
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
      }
      else
      {
         throw std::logic_error("Primitive backwardVector for SphereBessel is not implemented");
      }

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardVGradient(FieldComponents::Spectral::Id id, const PhysPathId& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         throw std::logic_error("vector gradient not implemented for SphereBessel");
      } else
      {
         throw std::logic_error("vector gradient of primitive variables not implemented for SphereBessel");
      }

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardCurl(const PhysPathId& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         std::size_t rTorId;
         std::size_t rPolId;

         auto cR = req.find(FieldComponents::Physical::R)->second;
         if(cR == Path::ValueTorPol::id() || cR == Path::NoSlipTorPol::id())
         {
            rTorId = Backward::ValueOverr1::id();
         }
         else if(cR == Path::StressFreeTorPol::id())
         {
            rTorId = Backward::StressFreeOverr1::id();
         }
         else if(cR == Path::InsulatingTorPol::id())
         {
            rTorId = Backward::InsulatingOverr1::id();
         }
         else
         {
            throw std::logic_error("Unknown path transform path");
         }
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(rTorId);
         transform.back().addEdge(Backward::Laplh::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());

         auto cT = req.find(FieldComponents::Physical::THETA)->second;
         if(cT == Path::ValueTorPol::id() || cT == Path::NoSlipTorPol::id())
         {
            rTorId = Backward::ValueOverr1D1R1::id();
            rPolId = Backward::ValueSlapl::id();
         }
         else if(cT == Path::StressFreeTorPol::id())
         {
            rTorId = Backward::StressFreeOverr1D1R1::id();
            rPolId = Backward::ValueSlapl::id();
         }
         else if(cT == Path::InsulatingTorPol::id())
         {
            rTorId = Backward::InsulatingOverr1D1R1::id();
            rPolId = Backward::ValueSlapl::id();
         }
         else
         {
            throw std::logic_error("Unknown path transform path");
         }
         // Toroidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(rTorId);
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

         // Poloidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(rPolId);
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());

         auto cP = req.find(FieldComponents::Physical::PHI)->second;
         if(cP == Path::ValueTorPol::id() || cP == Path::NoSlipTorPol::id())
         {
            rTorId = Backward::ValueOverr1D1R1::id();
            rPolId = Backward::ValueSlapl::id();
         }
         else if(cP == Path::StressFreeTorPol::id())
         {
            rTorId = Backward::StressFreeOverr1D1R1::id();
            rPolId = Backward::ValueSlapl::id();
         }
         else if(cP == Path::InsulatingTorPol::id())
         {
            rTorId = Backward::InsulatingOverr1D1R1::id();
            rPolId = Backward::ValueSlapl::id();
         }
         else
         {
            throw std::logic_error("Unknown path transform path");
         }
         // Toroidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(rTorId);
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());

         // Poloidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(rPolId);
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
      }
      else
      {
         throw std::logic_error("Curl of primitive variables not implemented for SphereBessel");
      }

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardDivergence() const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         // The divergence is zero be construction in this case!
         throw std::logic_error("Divergence should not be used in Toroidal/Poloidal expansion");
      } else
      {
         throw std::logic_error("Divergence of primitve variables not implemented for SphereBessel");
      }

      return transform;
   }

} // Transform
} // QuICC
