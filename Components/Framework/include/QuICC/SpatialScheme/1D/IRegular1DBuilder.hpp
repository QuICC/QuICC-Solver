/** 
 * @file IRegular1DBuilder.hpp
 * @brief Implementation of a generic regular 1D scheme
 */

#ifndef QUICC_SPATIALSCHEME_IREGULAR1DBUILDER_HPP
#define QUICC_SPATIALSCHEME_IREGULAR1DBUILDER_HPP

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
    * @brief Implementation of a generic regular 1D scheme
    */
   class IRegular1DBuilder: public IBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations 
          * @param purpose Setup purpose: simulation, visualization
          * @param options Options for builder
          */
         explicit IRegular1DBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         virtual ~IRegular1DBuilder() = default;

         /**
          * @brief Get spatial scheme dimensions
          */
         virtual ArrayI resolution() const override;

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
         virtual int fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::vector<int>& id, const std::vector<int>& bins) override;

         /**
          * @brief Get total of splittable indexes 
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag) override;
         
      protected:
         /**
          * @brief Get truncation tool
          */
         std::shared_ptr<Tools::IBase> truncationTools(const Dimensions::Transform::Id transId) const;

         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions() override;

         /**
          * @brief Set transform costs
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setCosts() override;

         /**
          * @brief Set transform scalings
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setScalings() override;

         /**
          * @brief Set transform memory footprint
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setMemoryScore() override;

         /**
          * @brief First truncation
          */
         int   mI;

      private:
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_IREGULAR1DBUILDER_HPP
