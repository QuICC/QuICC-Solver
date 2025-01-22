/**
 * @file ShellTransformSteps.cpp
 * @brief Source of the implementation of the physical <-> spectral transform steps in a spherical shell
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/ShellTransformSteps.hpp"

// Project includes
//
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Path/Scalar.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"
#include "QuICC/Transform/Path/ScalarNl.hpp"
#include "QuICC/Transform/Path/CurlNl.hpp"
#include "QuICC/Transform/Path/CurlCurlNl.hpp"
#include "QuICC/Transform/Path/I2ScalarNl.hpp"
#include "QuICC/Transform/Path/I2CurlNl.hpp"
#include "QuICC/Transform/Path/I2CurlCurlNl.hpp"
#include "QuICC/Transform/Path/NegI2CurlCurlNl.hpp"
#include "QuICC/Transform/Path/NegI2rCurlCurlNl.hpp"
#include "QuICC/Transform/Path/NegI4CurlCurlNl.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/I2P.hpp"
#include "QuICC/Transform/Forward/Overlaplh.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/OverlaplhD1.hpp"
#include "QuICC/Transform/Forward/OverlaplhOversinDphi.hpp"
#include "QuICC/Transform/Forward/I4Q.hpp"
#include "QuICC/Transform/Forward/I4S.hpp"
#include "QuICC/Transform/Forward/I2Q.hpp"
#include "QuICC/Transform/Forward/I2S.hpp"
#include "QuICC/Transform/Forward/I2rQ.hpp"
#include "QuICC/Transform/Forward/I2rS.hpp"
#include "QuICC/Transform/Forward/I2T.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/Overr2.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1Laplh.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/OversinDphi.hpp"
#include "QuICC/Transform/Backward/OversinLaplhDphi.hpp"
#include "QuICC/Transform/Backward/Slaplr.hpp"
#include "QuICC/Transform/Backward/OversinD1Sin.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"

namespace QuICC {

namespace Transform {

   ShellTransformSteps::ShellTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
      : ITransformSteps(spScheme)
   {
   }

   bool ShellTransformSteps::applicable(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
   {
      return spScheme->has(SpatialScheme::Feature::ShellGeometry);
   }

   std::vector<TransformPath>  ShellTransformSteps::forwardScalar(const std::vector<PathId >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;
      auto flag = components.at(0).second;

      if(flag == Path::Scalar::id())
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), scalId, Arithmetics::Add::id());
      }
      else
      {
         throw std::logic_error("Requested an unknown scalar forward transform");
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::forwardNLScalar(const std::vector<PathId >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;
      auto flag = components.at(0).second;

      // Without quasi-inverse
      if(flag == Path::ScalarNl::id())
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), scalId, Arithmetics::Add::id());
      }
      else
      {
         throw std::logic_error("Requested an unknown nonlinear scalar forward transform");
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::forwardVector(const std::vector<PathId >& components) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         assert(components.size() == 2);
         FieldComponents::Spectral::Id curlId = components.at(0).first;
         auto curlFlag = components.at(0).second;
         FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
         auto curlcurlFlag = components.at(1).second;

         if(curlFlag == Path::TorPol::id() && curlcurlFlag == Path::TorPol::id())
         {
            // Compute Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::P::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::P::id(), curlId, Arithmetics::Sub::id());

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Overlaplh::id());
            transform.back().addEdge(Forward::Pol::id(), curlcurlId, Arithmetics::Add::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }
      }
      else
      {
         assert(components.size() == 3);

         transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::R, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::THETA, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::PHI, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::forwardNLVector(const std::vector<PathId >& components) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         assert(components.size() == 2);
         FieldComponents::Spectral::Id curlId = components.at(0).first;
         auto curlFlag = components.at(0).second;
         FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
         auto curlcurlFlag = components.at(1).second;

         // Integrate for standard second order equation
         if(curlFlag == Path::I2CurlNl::id())
         {
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::I2T::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::I2T::id(), curlId, Arithmetics::Sub::id());
         } else
         {
            throw std::logic_error("Requested an unknown curl nonlinear vector forward transform");
         }

         // Integrate for standard fourth order spherical equation with negative sign
         if(curlcurlFlag == Path::NegI4CurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I4Q::id(), curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::I4S::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::I4S::id(), curlcurlId, Arithmetics::Add::id());
         }
         // Integrate for standard second order spherical equation with negative sign
         else if(curlcurlFlag == Path::NegI2CurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I2Q::id(), curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::I2S::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::I2S::id(), curlcurlId, Arithmetics::Add::id());
         }
         // Integrate for standard second order spherical equation with negative sign and an additional factor of r
         else if(curlcurlFlag == Path::NegI2rCurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I2rQ::id(), curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::I2rS::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::I2rS::id(), curlcurlId, Arithmetics::Add::id());
         }
         // Integrate for second order spherical equation
         else if(curlcurlFlag == Path::I2CurlCurlNl::id())
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I2Q::id(), curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhD1::id());
            transform.back().addEdge(Forward::I2S::id(), curlcurlId, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::OverlaplhOversinDphi::id());
            transform.back().addEdge(Forward::I2S::id(), curlcurlId, Arithmetics::Sub::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }
      } else
      {
         assert(components.size() == 3);

         transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::R, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::THETA, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::PHI, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Backward::P::id());
      transform.back().addEdge(Backward::P::id());
      transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Overr1::id());
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::Overr1::id());
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req) const
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      throw std::logic_error("Second derivative is not implementated yet!");

      pairId = std::make_pair(FieldComponents::Physical::R,FieldComponents::Physical::R);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::R,FieldComponents::Physical::THETA);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::R,FieldComponents::Physical::PHI);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::THETA,FieldComponents::Physical::THETA);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::THETA,FieldComponents::Physical::PHI);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::PHI,FieldComponents::Physical::PHI);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Overr1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::Overr1D1R1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            // Poloidal part 1
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Slaplr::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());

            // Poloidal part 2
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Overr2::id());
            transform.back().addEdge(Backward::OversinLaplhDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1D1R1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Slaplr::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::Overr2::id());
            transform.back().addEdge(Backward::D1Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::OversinD1Sin::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::OversinDphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
            transform.back().addEdge(Backward::Overr1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  ShellTransformSteps::backwardDivergence() const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         // The divergence is zero be construction in this case!
         throw std::logic_error("Divergence should not be used in Toroidal/Poloidal expansion");
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::Overr1::id());
         transform.back().addEdge(Backward::OversinD1Sin::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::Overr1::id());
         transform.back().addEdge(Backward::OversinDphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
      }

      return transform;
   }

}
}
