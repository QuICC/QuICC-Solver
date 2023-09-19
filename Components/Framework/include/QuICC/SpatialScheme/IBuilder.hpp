/**
 * @file IBuilder.hpp
 * @brief Implementation of the basic components of the spatial scheme builder
 */

#ifndef QUICC_SPATIALSCHEME_IBUILDER_HPP
#define QUICC_SPATIALSCHEME_IBUILDER_HPP

// System includes
//
#include <assert.h>
#include <vector>
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/SpatialScheme/ICosts.hpp"
#include "QuICC/SpatialScheme/IMesher.hpp"
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
          * @param purpose Grid purpose
          * @param options Scheme options
          */
         explicit IBuilder(const int dims, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         virtual ~IBuilder() = default;

         /**
          * @brief Tune the shared resolution used by simulation
          */
         virtual void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr);

         /**
          * @brief Initialise the scheme
          */
         void init();

         /**
          * @brief Create indexes for a possibly restricted set (to simplify implementation it is defined for 3D cases)
          *
          * @param transId Transform ID
          * @param fwd1D   Storage for forward indexes of first dimension (the one acted on by the transform and contigous in memory)
          * @param bwd1D   Storage for backward indexes of first dimension (the one acted on by the transform and contigous in memory)
          * @param idx2D   Storage for the indexes of second dimension
          * @param idx3D   Storage for forward indexes of third dimension
          * @param id      ID of the bin (process/rank)
          * @param bins    Total number of bins (useful to build efficient pairs)
          */
         virtual int fillIndexes(const Dimensions::Transform::Id transId, std::vector<std::vector<std::vector<int> > >& fwd1D, std::vector<std::vector<std::vector<int> > >& bwd1D, std::vector<std::vector<int> >& idx2D, std::vector<int>& idx3D, const std::vector<int>& id, const std::vector<int>& bins) = 0;

         /**
          * @brief Get total of splittable indexes
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag) = 0;

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

         /**
          * @brief Mesher
          */
         void  setMesher(std::shared_ptr<IMesher> m, const bool isCustom);

         /**
          * @brief Spectral space and transform have different ordering?
          */
         virtual bool sameSpectralOrdering() const;

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

         /**
          * @brief Mesher
          */
         const IMesher& mesher() const;

         /**
          * @brief Mesher
          */
         IMesher& mesher();

         /**
          * @brief Scheme options (Poly vs FFT, Uniform vs Triangular, etc)
          */
         std::map<std::size_t,std::vector<std::size_t>> mOptions;

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

         /**
          * @brief Mesher
          */
         std::shared_ptr<IMesher> mspMesher;
   };

   /// Typedef for a shared pointer to a IBuilder object
   typedef std::shared_ptr<IBuilder>   SharedIBuilder;

   /**
    * @brief Template builder factory
    */
   template <typename TBuilder> std::shared_ptr<TBuilder> makeBuilder(ArrayI& dim, const GridPurpose::Id purpose, const bool needInterpretation, const std::map<std::size_t,std::vector<std::size_t>>& options, std::shared_ptr<IMesher> m);

   //
   //
   //
   template <typename TBuilder> std::shared_ptr<TBuilder> makeBuilder(ArrayI& dim, const GridPurpose::Id purpose, const bool needInterpretation, const std::map<std::size_t,std::vector<std::size_t>>& options, std::shared_ptr<IMesher> m)
   {
      if(needInterpretation)
      {
         TBuilder::interpretConfigDimensions(dim);
      }

      std::shared_ptr<TBuilder> spBuilder = std::make_shared<TBuilder>(dim, purpose, options);
      if(m)
      {
         spBuilder->setMesher(m, true);
      }
      spBuilder->init();
      return spBuilder;
   }

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_IBUILDER_HPP
