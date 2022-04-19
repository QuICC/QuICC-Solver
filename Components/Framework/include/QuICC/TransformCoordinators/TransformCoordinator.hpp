/**
 * @file TransformCoordinator.hpp
 * @brief Implementation of a transform coordinator
 */

#ifndef QUICC_TRANSFORMCOORDINATOR_HPP
#define QUICC_TRANSFORMCOORDINATOR_HPP

// Debug includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//
#include<set>
#include<map>
#include<memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/Transform/ITransform.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

   /**
    *
    * @brief Implementation of a transform coordinator
    *
    * The transform coordinator overlooks the whole transform process. Providing storage, communication
    * and the actual transforms
    */
   template <typename TComm> class TransformCoordinator
   {
      public:
         /// Typedef for the communicator type
         typedef TComm CommunicatorType;

         /**
          * @brief Very basic constructor
          */
         TransformCoordinator();

         /**
          * @brief Destructor
          */
         virtual ~TransformCoordinator();

         /**
          * @brief Add a transform
          */
         void addTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::ITransform>  spTransform);

         /**
          * @brief Initialise transforms
          *
          * @param spRes         Resolution information object
          * @param forwardTree   Transform forward tree
          * @param backwardTree  Transform backward tree
          */
         void defineTransforms(const std::vector<Transform::TransformTree>& forwardTree, const std::vector<Transform::TransformTree>& backwardTree, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Initialise the data communicator
          *
          * @param spRes   Resolution information object
          */
         void initCommunicator(SharedResolution spRes);

         /**
          * @brief Get list of required options for transform
          */
         void requiredOptions(std::set<std::size_t>& list) const;

         /**
          * @brief Set options of transform
          */
         void setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options);

         /**
          * @brief Get the communicator
          */
         CommunicatorType&  communicator();

         /**
          * @brief Get the transform integrator tree
          */
         const std::vector<Transform::TransformTree>& forwardTree() const;

         /**
          * @brief Get the transform projector tree
          */
         const std::vector<Transform::TransformTree>& backwardTree() const;

         /**
          * @brief Get grid array(s) of the mesh
          */
         std::vector<Array>  mesh() const;

         /**
          * @brief Get sahred grid array(s) of the mesh
          */
         std::shared_ptr<std::vector<Array> >  spMesh();

         /**
          * @brief Get the communicator
          */
         Transform::ITransform&  transform(const int id);

         /**
          * @brief Get the communicator
          */
         Transform::ITransform&  transform(const Dimensions::Transform::Id id);

         /**
          * @brief Get the communicator
          */
         Transform::ITransform&  transform1D();

         /**
          * @brief Get the communicator
          */
         Transform::ITransform&  transform2D();

         /**
          * @brief Get the communicator
          */
         Transform::ITransform&  transform3D();

         /**
          * @brief Get the communicator
          */
         Transform::ITransform&  transformND();

         /**
          * @brief Get spatial scheme
          */
         const SpatialScheme::ISpatialScheme& ss() const;

         /**
          * @brief Profile storage used by transform coordinator
          */
         void profileStorage() const;

      protected:

      private:
         /**
          * @brief Shared spatial scheme
          */
         SpatialScheme::SharedCISpatialScheme mspSS;

         /**
          * @brief The communicator
          */
         CommunicatorType  mCommunicator;

         /**
          * @brief The transforms
          */
         std::map<Dimensions::Transform::Id,std::shared_ptr<Transform::ITransform> >  mTransforms;

         /**
          * @brief Transform backward tree
          */
         std::vector<Transform::TransformTree> mBackwardTree;

         /**
          * @brief Transform forward tree
          */
         std::vector<Transform::TransformTree> mForwardTree;
   };

   template <typename TComm>
      inline const SpatialScheme::ISpatialScheme& TransformCoordinator<TComm>::ss() const
   {
      return *this->mspSS;
   }

   template <typename TComm>
      inline Transform::ITransform& TransformCoordinator<TComm>::transform(const Dimensions::Transform::Id id)
   {
      return *this->mTransforms.find(id)->second;
   }

   template <typename TComm>
      inline Transform::ITransform&  TransformCoordinator<TComm>::transform1D()
   {
      return this->transform(Dimensions::Transform::TRA1D);
   }

   template <typename TComm>
      inline Transform::ITransform&  TransformCoordinator<TComm>::transform2D()
   {
      return this->transform(Dimensions::Transform::TRA2D);
   }

   template <typename TComm>
      inline Transform::ITransform& TransformCoordinator<TComm>::transform3D()
   {
      return this->transform(Dimensions::Transform::TRA3D);
   }

   template <typename TComm>
      inline Transform::ITransform&  TransformCoordinator<TComm>::transformND()
   {
      return *this->mTransforms.rbegin()->second;
   }

   template <typename TComm>
      inline typename TransformCoordinator<TComm>::CommunicatorType&  TransformCoordinator<TComm>::communicator()
   {
      return this->mCommunicator;
   }

   template <typename TComm>
      inline const std::vector<Transform::TransformTree>& TransformCoordinator<TComm>::forwardTree() const
   {
      return this->mForwardTree;
   }

   template <typename TComm>
      inline const std::vector<Transform::TransformTree>& TransformCoordinator<TComm>::backwardTree() const
   {
      return this->mBackwardTree;
   }

   template <typename TComm>
      TransformCoordinator<TComm>::TransformCoordinator()
   {
   }

   template <typename TComm>
      TransformCoordinator<TComm>::~TransformCoordinator()
   {
   }

   template <typename TComm>
      void TransformCoordinator<TComm>::addTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::ITransform> spTransform)
   {
      this->mTransforms.insert(std::make_pair(id, spTransform));
   }

   template <typename TComm>
      void TransformCoordinator<TComm>::defineTransforms(const std::vector<Transform::TransformTree>& forwardTree, const std::vector<Transform::TransformTree>& backwardTree, SpatialScheme::SharedCISpatialScheme spScheme)
   {
      // Store the forward transform tree
      this->mForwardTree = forwardTree;

      // Store the backward transform tree
      this->mBackwardTree = backwardTree;

      // Store spatial sscheme
      this->mspSS = spScheme;
   }

   template <typename TComm>
      void TransformCoordinator<TComm>::initCommunicator(SharedResolution spRes)
   {
      // initialise the communicator
      for(auto it = this->mTransforms.cbegin(); it != this->mTransforms.cend(); ++it)
      {
         this->mCommunicator.init(it->first, spRes->spFwdSetup(it->first), spRes->spBwdSetup(it->first));
      }
   }

   template <typename TComm>
      void TransformCoordinator<TComm>::requiredOptions(std::set<std::size_t>& list) const
   {
      assert(this->mTransforms.size() > 0);

      for(auto it = this->mTransforms.cbegin(); it != this->mTransforms.cend(); ++it)
      {
         it->second->requiredOptions(list, it->first);
      }
   }

   template <typename TComm>
      void TransformCoordinator<TComm>::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options)
   {
      assert(this->mTransforms.size() > 0);

      for(auto it = this->mTransforms.begin(); it != this->mTransforms.end(); ++it)
      {
         it->second->setOptions(options, it->first);
      }
   }

   template <typename TComm>
      std::vector<Array> TransformCoordinator<TComm>::mesh() const
   {
      std::vector<Array>   mesh;

      // Add grid of first dimension
      for(auto it = this->mTransforms.cbegin(); it != this->mTransforms.cend(); ++it)
      {
         mesh.push_back(it->second->meshGrid());
      }

      return mesh;
   }

   template <typename TComm>
      void TransformCoordinator<TComm>::profileStorage() const
   {
      #ifdef QUICC_STORAGEPROFILE

      for(auto it = this->mTransforms.cbegin(); it != this->mTransforms.cend(); ++it)
      {
         it->second->profileStorage();
      }

      #endif // QUICC_STORAGEPROFILE
   }

}

#endif // QUICC_TRANSFORMCOORDINATOR_HPP
