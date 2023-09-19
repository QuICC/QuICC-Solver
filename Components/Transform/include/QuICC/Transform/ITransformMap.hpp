/**
 * @file ITransformMap.hpp
 * @brief Implementation of the Worland transform in a sphere
 */

#ifndef QUICC_TRANSFORM_ITRANSFORMMAP_HPP
#define QUICC_TRANSFORM_ITRANSFORMMAP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Worland transform in a sphere
    */
   template <typename TMapOp> class ITransformMap
   {
      public:
         typedef std::map<std::size_t,std::shared_ptr<TMapOp> > MapType;

         /**
          * @brief Constructor
          */
         ITransformMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~ITransformMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         virtual void operator()(MapType& m) const = 0;

         /**
          * @brief Map transform operator to ID
          */
         template <typename TOp> void addOperator(MapType& m, const std::size_t id = 0) const;

      private:
   };

   template <typename TMapOp> template <typename TOp> void ITransformMap<TMapOp>::addOperator(MapType& m, const std::size_t id) const
   {
      std::size_t lid;

      if(id == 0)
      {
         lid = m.size();
      }
      else
      {
         lid = id;
      }

      std::shared_ptr<TOp> spOp = std::make_shared<TOp>();
      auto ret = m.try_emplace(lid, spOp);

      if(!ret.second)
      {
         throw std::logic_error("Operator could not be added to transform operator map");
      }
   }

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_ITRANSFORMMAP_HPP
