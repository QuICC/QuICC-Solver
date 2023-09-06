/**
 * @file IVariableAsciiWriter.cpp
 * @brief Source of the implementation of the generic variable to ASCII file writer
 */

// System includes
//

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Hasher.hpp"
#include "QuICC/NonDimensional/Coordinator.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   IVariableAsciiWriter::IVariableAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode)
      : IAsciiWriter(name, ext, header, type, version, mode), mTime(-1.0), mTimestep(-1.0), mSpaceId(id)
   {
   }

   Dimensions::Space::Id IVariableAsciiWriter::space() const
   {
      return this->mSpaceId;
   }

   void IVariableAsciiWriter::setPhysical(const std::map<std::string,MHDFloat>& parameters, const std::map<std::string,std::size_t>& boundary)
   {
      // Convert parameters to NonDimensional numbers
      for(auto it = parameters.cbegin(); it != parameters.cend(); ++it)
      {
         size_t nd = Hasher::makeId(it->first);
         this->mPhysical.insert(std::make_pair(nd,NonDimensional::Coordinator::map().find(nd)->second->create(it->second)));
      }

      // Save boundary flags
      this->mBoundary = boundary;
   }

   const Resolution& IVariableAsciiWriter::res() const
   {
      return *this->mspRes;
   }

   void IVariableAsciiWriter::setMesh(const std::vector<Array>& mesh)
   {
      this->mMesh = mesh;
   }

   void IVariableAsciiWriter::setSimTime(const MHDFloat time, const MHDFloat timestep)
   {
      this->mTime = time;

      this->mTimestep = timestep;
   }

   void IVariableAsciiWriter::setResolution(SharedResolution spRes)
   {
      // Store resolution object
      this->mspRes = spRes;
   }

   void IVariableAsciiWriter::expect(const std::size_t id)
   {
      this->mExpected.insert(id);
   }

   bool IVariableAsciiWriter::isFull() const
   {
      bool status = true;

      // Check that all expected scalars and vectors are present
      status = status && (this->mScalars.size() + this->mVectors.size() == this->mExpected.size());

      // Check that the resolution has been set
      status = status && this->mspRes;

      return status;
   }

   void IVariableAsciiWriter::addScalar(const std::pair<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalar)
   {
      // Only add the variable if it is listed as expected
      if(this->mExpected.count(scalar.first))
      {
         this->mScalars.insert(scalar);

         // Resolution is not yet set extract from scalar
         if(!this->mspRes)
         {
            std::visit([&](auto&& p){this->setResolution(p->dom(0).spRes());}, scalar.second);
         }
      }
   }

   void IVariableAsciiWriter::addVector(const std::pair<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vector)
   {
      // Only add the variable if it is listed as expected
      if(this->mExpected.count(vector.first))
      {
         this->mVectors.insert(vector);

         // Resolution is not yet set extract from vector
         if(!this->mspRes)
         {
            std::visit([&](auto&& p){this->setResolution(p->dom(0).spRes());}, vector.second);
         }
      }
   }

   IVariableAsciiWriter::scalar_iterator_range  IVariableAsciiWriter::scalarRange()
   {
      return std::make_pair(this->mScalars.begin(), this->mScalars.end());
   }

   IVariableAsciiWriter::vector_iterator_range  IVariableAsciiWriter::vectorRange()
   {
      return std::make_pair(this->mVectors.begin(), this->mVectors.end());
   }

   void IVariableAsciiWriter::compute(Transform::TransformCoordinatorType&)
   {
   }

}
}
}
