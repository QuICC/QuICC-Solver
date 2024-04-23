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
#include "QuICC/Transform/Path/Value/Scalar.hpp"
#include "QuICC/Transform/Path/Insulating/Scalar.hpp"
#include "QuICC/Transform/Path/Value/TorPol.hpp"
#include "QuICC/Transform/Path/Insulating/TorPol.hpp"
#include "QuICC/Transform/Path/Value/Tor.hpp"
#include "QuICC/Transform/Path/Insulating/Tor.hpp"
#include "QuICC/Transform/Path/Value/Pol.hpp"
#include "QuICC/Transform/Path/Insulating/Pol.hpp"
#include "QuICC/Transform/Path/Value/ScalarNl.hpp"
#include "QuICC/Transform/Path/Insulating/ScalarNl.hpp"
#include "QuICC/Transform/Path/Value/CurlNl.hpp"
#include "QuICC/Transform/Path/Insulating/CurlNl.hpp"
#include "QuICC/Transform/Path/Value/CurlCurlNl.hpp"
#include "QuICC/Transform/Path/Insulating/CurlCurlNl.hpp"
#include "QuICC/Transform/Path/Value/NegCurlCurlNl.hpp"
#include "QuICC/Transform/Path/Insulating/NegCurlCurlNl.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/Value/P.hpp"
#include "QuICC/Transform/Forward/Insulating/P.hpp"
#include "QuICC/Transform/Forward/Overlaplh.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/OverlaplhD1.hpp"
#include "QuICC/Transform/Forward/OverlaplhOversinDphi.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/Value/Pol.hpp"
#include "QuICC/Transform/Forward/Insulating/Pol.hpp"
#include "QuICC/Transform/Forward/Value/Q.hpp"
#include "QuICC/Transform/Forward/Insulating/Q.hpp"
#include "QuICC/Transform/Forward/Value/S.hpp"
#include "QuICC/Transform/Forward/Insulating/S.hpp"
#include "QuICC/Transform/Forward/Value/T.hpp"
#include "QuICC/Transform/Forward/Insulating/T.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Value/P.hpp"
#include "QuICC/Transform/Backward/Insulating/P.hpp"
#include "QuICC/Transform/Backward/Value/Overr1.hpp"
#include "QuICC/Transform/Backward/Insulating/Overr1.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/Value/D1.hpp"
#include "QuICC/Transform/Backward/Insulating/D1.hpp"
#include "QuICC/Transform/Backward/Value/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Insulating/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/OversinDphi.hpp"
#include "QuICC/Transform/Backward/Value/Slapl.hpp"
#include "QuICC/Transform/Backward/Insulating/Slapl.hpp"
#include "QuICC/Transform/Backward/OversinD1Sin.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"

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

      // Value BC
      if(flag == Path::Value::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::Value::P::id(), scalId, Arithmetics::Add::id());
      }
      // Insulating BC
      else if(flag == Path::Insulating::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::Insulating::P::id(), scalId, Arithmetics::Add::id());
      }
      else
      {
         throw std::logic_error("Requested an unknown scalar forward transform (ID = " + std::to_string(flag) + ")");
      }

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::forwardNLScalar(const std::vector<PathId >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;
      auto flag = components.at(0).second;

      DebuggerMacro_msg("Using SphereBesselTransformSteps for forwardNLScalar", 1);

      // Value BC
      if(flag == Path::Value::ScalarNl::id())
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::Value::P::id(), scalId, Arithmetics::Add::id());
      }
      // Insulating BC
      else if(flag == Path::Insulating::ScalarNl::id())
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::Insulating::P::id(), scalId, Arithmetics::Add::id());
      }
      else
      {
         throw std::logic_error("Requested an unknown nonlinear scalar forward transform (ID = " + std::to_string(flag) + ")");
      }

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
         if(curlFlag == Path::Value::TorPol::id() && curlcurlFlag == Path::Value::TorPol::id())
         {
            // Compute Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Value::P::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Value::P::id(), curlId, Arithmetics::Sub::id());

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Overlaplh::id());
            transform.back().addEdge(Forward::Value::Pol::id(), curlcurlId, Arithmetics::Add::id());
         }
         // Value BC for toroidal and Insulating BC for poloidal
         else if(curlFlag == Path::Value::TorPol::id() && curlcurlFlag == Path::Insulating::TorPol::id())
         {
            // Compute Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Value::P::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Value::P::id(), curlId, Arithmetics::Sub::id());

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Overlaplh::id());
            transform.back().addEdge(Forward::Insulating::Pol::id(), curlcurlId, Arithmetics::Add::id());
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

         // Value BC
         if(curlFlag == Path::Value::CurlNl::id())
         {
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Value::T::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Value::T::id(), curlId, Arithmetics::Sub::id());
         }
         // Insulating BC
         else if(curlFlag == Path::Insulating::CurlNl::id())
         {
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Insulating::T::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Insulating::T::id(), curlId, Arithmetics::Sub::id());
         }
         else
         {
            throw std::logic_error("Requested an unknown curl nonlinear vector forward transform (ID = " + std::to_string(curlFlag) + ")");
         }

         // Value BC 
         if(curlcurlFlag == Path::Value::CurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Value::Q::id(), curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Value::S::id(), curlcurlId, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Value::S::id(), curlcurlId, Arithmetics::Sub::id());
         }
         // Insulating BC
         else if(curlcurlFlag == Path::Insulating::CurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Insulating::Q::id(), curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Insulating::S::id(), curlcurlId, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Insulating::S::id(), curlcurlId, Arithmetics::Sub::id());
         }
         // Negative, Value BC
         else if(curlcurlFlag == Path::Value::NegCurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Value::Q::id(), curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Value::S::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Value::S::id(), curlcurlId, Arithmetics::Add::id());
         }
         // Negative, Insulating BC
         else if(curlcurlFlag == Path::Insulating::NegCurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Insulating::Q::id(), curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::Insulating::S::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::Insulating::S::id(), curlcurlId, Arithmetics::Add::id());
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

      if(req.find(FieldComponents::Physical::SCALAR)->second == Path::Value::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Backward::Value::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
      }
      else if(req.find(FieldComponents::Physical::SCALAR)->second == Path::Insulating::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Backward::Insulating::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  SphereBesselTransformSteps::backwardGradient(const PhysPathId& req) const
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second == Path::Value::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Value::D1::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
      }
      else if(req.find(FieldComponents::Physical::R)->second == Path::Insulating::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Insulating::D1::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
      }

      if(req.find(FieldComponents::Physical::THETA)->second == Path::Value::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Value::Overr1::id());
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
      }
      else if(req.find(FieldComponents::Physical::THETA)->second == Path::Insulating::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Insulating::Overr1::id());
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
      }

      if(req.find(FieldComponents::Physical::PHI)->second == Path::Value::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Value::Overr1::id());
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
      }
      else if(req.find(FieldComponents::Physical::PHI)->second == Path::Insulating::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Insulating::Overr1::id());
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
      }

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
         if(req.find(FieldComponents::Physical::R)->second == Path::Value::Pol::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Value::Overr1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::R)->second == Path::Insulating::Pol::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Insulating::Overr1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second == Path::Value::Tor::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::Value::P::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::THETA)->second == Path::Insulating::Tor::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::Insulating::P::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::THETA)->second == Path::Value::Pol::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Value::Overr1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::THETA)->second == Path::Insulating::Pol::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Insulating::Overr1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second == Path::Value::Tor::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::Value::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());
         }
         else if(req.find(FieldComponents::Physical::PHI)->second == Path::Insulating::Tor::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::Insulating::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());
         }
         else if(req.find(FieldComponents::Physical::PHI)->second == Path::Value::Pol::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Value::Overr1D1R1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::PHI)->second == Path::Insulating::Pol::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Insulating::Overr1D1R1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
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
         if(req.find(FieldComponents::Physical::R)->second == Path::Value::Tor::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Value::Overr1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::R)->second == Path::Insulating::Tor::id())
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Insulating::Overr1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second == Path::Value::Tor::id())
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Value::Overr1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::THETA)->second == Path::Insulating::Tor::id())
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Insulating::Overr1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::THETA)->second == Path::Value::Pol::id())
         {
            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Value::Slapl::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());
         }
         else if(req.find(FieldComponents::Physical::THETA)->second == Path::Insulating::Pol::id())
         {
            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Insulating::Slapl::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second == Path::Value::Tor::id())
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Value::Overr1D1R1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::PHI)->second == Path::Insulating::Tor::id())
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Insulating::Overr1D1R1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::PHI)->second == Path::Value::Pol::id())
         {
            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Value::Slapl::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
         else if(req.find(FieldComponents::Physical::PHI)->second == Path::Insulating::Pol::id())
         {
            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Insulating::Slapl::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      } else
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
