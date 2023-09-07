/** 
 * @file I3DBuilder.hpp
 * @brief Implementation of a generic regular 3D scheme
 */

#ifndef QUICC_SPATIALSCHEME_I3DBUILDER_HPP
#define QUICC_SPATIALSCHEME_I3DBUILDER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of a generic regular 3D scheme
    */
   class I3DBuilder: public IBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations 
          * @param purpose Setup purpose: simulation, visualization
          * @param options Options for builder
          */
         explicit I3DBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         virtual ~I3DBuilder() = default;

         /**
          * @brief Create indexes for a possibly restricted set
          *
          * @param transId Transform ID
          * @param fwd1D   Storage for forward indexes of first dimension
          * @param bwd1D   Storage for backward indexes of first dimension
          * @param idx2D   Storage for the indexes of second dimension
          * @param idx3D   Storage for forward indexes of third dimension
          * @param id      ID of the bin
          * @param bins    Total number of bins (useful to build efficient pairs)
          */
         virtual int fillIndexes(const Dimensions::Transform::Id transId, std::vector<std::vector<std::vector<int> > >& fwd1D, std::vector<std::vector<std::vector<int> > >& bwd1D, std::vector<std::vector<int> >& idx2D, std::vector<int>& idx3D, const std::vector<int>& id, const std::vector<int>& bins) override;

         /**
          * @brief Get total of splittable indexes 
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag) override;
         
      protected:
         /**
          * @brief Initialise the domain dimensions with default setup
          */
         virtual void setDimensions() override;

         /**
          * @brief Validate mode splitting
          *
          * @param modes   Modes
          * @param transId Transform ID
          */
         virtual int checkModes(const std::multimap<int,int>& modes, const Dimensions::Transform::Id transId) const;

         /**
          * @brief Get truncation tool
          */
         virtual std::shared_ptr<Tools::IBase> truncationTools(const Dimensions::Transform::Id transId) const = 0;

         /**
          * @brief Compute mode distribution
          *
          * @param modes   Map of all modes
          * @param id      List of local ID in splitting groups
          * @param bins    Number of bins in each spitting group
          * @param transId Transform stage to compute splittin for
          */
         void split(std::multimap<int,int>& modes, const std::vector<int>& id, const std::vector<int>& bins, const Dimensions::Transform::Id transId);

         /**
          * @brief Set transform costs
          */
         virtual void setCosts() override;

         /**
          * @brief Set transform scalings
          */
         virtual void setScalings() override;

         /**
          * @brief Set transform memory footprint
          */
         virtual void setMemoryScore() override;

      private:
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_I3DBUILDER_HPP
