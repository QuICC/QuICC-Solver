/** 
 * @file IBuilder.hpp
 * @brief Implementation of the basic components of the spatial scheme builder
 */

#ifndef QUICC_SPATIALSCHEME_IBUILDER_HPP
#define QUICC_SPATIALSCHEME_IBUILDER_HPP

// Configuration includes
//

// System includes
//
#include <assert.h>
#include <vector>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/SpatialScheme/ICosts.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the basic components of the spatial scheme builder
    */
   class IBuilder : public ICosts
   {
      public:
         /**
          * @brief Interpret the configuration dimensions
          */
         static void interpretConfigDimensions(ArrayI& dim);

         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         explicit IBuilder(const int dims, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~IBuilder();

         /**
          * @brief Tune the shared resolution used by simulation
          */
         virtual void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr);

         /**
          * @brief Initialise the scheme
          */
         void init();

         /**
          * @brief Create indexes for a possibly restricted set (to simplify implementation is is defined for 3D cases)
          *
          * @param transId Transform ID
          * @param fwd1D   Storage for forward indexes of first dimension
          * @param bwd1D   Storage for backward indexes of first dimension
          * @param idx2D   Storage for the indexes of second dimension
          * @param idx3D   Storage for forward indexes of third dimension
          * @param id      ID of the bin
          * @param bins    Total number of bins (useful to build efficient pairs)
          * @param n0      Starting index of restricted set
          * @param nN      Length of restricted set
          * @param flag    Flag to specify location of splitting
          */
         virtual int fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id = ArrayI(), const ArrayI& bins = ArrayI(), const ArrayI& n0 = ArrayI(), const ArrayI& nN = ArrayI(), const Splitting::Locations::Id flag = Splitting::Locations::NONE) = 0;

         /**
          * @brief Get total of splittable indexes 
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag) = 0;

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const = 0;

         /**
          * @brief Add scheme specific transform setups to resolution
          */
         virtual void addTransformSetups(SharedResolution spRes) const = 0;

         /**
          * @brief Get spatial scheme dimensions
          */
         virtual ArrayI resolution() const = 0;

         /**
          * @brief Get the simulation wide spectral array dimensions (can be different from spectral resolution)
          */
         const ArrayI& getTransformSpace() const;

         /**
          * @brief Add index counter to shared resolution
          */
         virtual void addIndexCounter(SharedResolution spRes);
         
      protected:
         /**
          * @brief Main purpose of the grid values
          */
         GridPurpose::Id purpose() const;

         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions() = 0;

         /**
          * @brief Get forward dimension of given transform
          *
          * @param i Transform index
          */
         int dim(const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId) const;

         /**
          * @brief Get forward dimension of given transform
          *
          * @param i Transform index
          */
         void setDimension(int d, const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId);

         /**
          * @brief Set transform dimension
          */
         void setTransformSpace(const ArrayI& dim);

         /**
          * @brief Get dimension of the domain
          */
         int dims() const;

         /**
          * @brief Tune resolution with MPI related conditions
          */
         void tuneMpiResolution(const Parallel::SplittingDescription& descr);

      private:
         /**
          * @brief Main purpose of the grid values
          */
         const GridPurpose::Id mPurpose;

         /**
          * @brief Dimension of the domain
          */
         int   mDims;

         /**
          * @brief Array size of the spectral arrays (can be different from spectral resolution)
          */
         ArrayI mTransformSpace;

         /**
          * @brief Full dimensions of the domain
          */
         std::vector<ArrayI>   mDimensions;
   };

   /// Typedef for a shared pointer to a IBuilder object
   typedef std::shared_ptr<IBuilder>   SharedIBuilder;

   /**
    * @brief Template builder factory
    */
   template <typename TBuilder> std::shared_ptr<TBuilder> makeBuilder(ArrayI& dim, const GridPurpose::Id purpose, const bool needInterpretation);

   //
   //
   //
   template <typename TBuilder> std::shared_ptr<TBuilder> makeBuilder(ArrayI& dim, const GridPurpose::Id purpose, const bool needInterpretation)
   {
      if(needInterpretation)
      {
         TBuilder::interpretConfigDimensions(dim);
      }

      std::shared_ptr<TBuilder> spBuilder = std::make_shared<TBuilder>(dim, purpose);
      spBuilder->init();
      return spBuilder;
   }
}
}

#endif // QUICC_SPATIALSCHEME_IBUILDER_HPP
