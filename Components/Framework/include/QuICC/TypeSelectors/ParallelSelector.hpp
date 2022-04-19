/**
 * @file ParallelSelector.hpp
 * @brief Typedefs for the parallelisation algorithms
 */

#ifndef QUICC_TRANSFORM_PARALLELSELECTOR_HPP
#define QUICC_TRANSFORM_PARALLELSELECTOR_HPP

// Configuration includes
//

// System includes
//
#include <memory>
#include <stdexcept>

// External includes
//

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"
#include "QuICC/TransformConfigurators/ForwardSerialConfigurator.hpp"
#include "QuICC/TransformConfigurators/BackwardSerialConfigurator.hpp"
#include "QuICC/TransformConfigurators/ForwardSingle1DConfigurator.hpp"
#include "QuICC/TransformConfigurators/BackwardSingle1DConfigurator.hpp"
#include "QuICC/TransformConfigurators/ForwardSingle2DConfigurator.hpp"
#include "QuICC/TransformConfigurators/BackwardSingle2DConfigurator.hpp"
#include "QuICC/TransformConfigurators/ForwardTubularConfigurator.hpp"
#include "QuICC/TransformConfigurators/BackwardTubularConfigurator.hpp"
#include "QuICC/TransformGroupers/ForwardEquationGrouper.hpp"
#include "QuICC/TransformGroupers/BackwardEquationGrouper.hpp"
#include "QuICC/TransformGroupers/ForwardSingle1DGrouper.hpp"
#include "QuICC/TransformGroupers/BackwardSingle1DGrouper.hpp"
#include "QuICC/TransformGroupers/ForwardSingle2DGrouper.hpp"
#include "QuICC/TransformGroupers/BackwardSingle2DGrouper.hpp"
#include "QuICC/TransformGroupers/ForwardTransformGrouper.hpp"
#include "QuICC/TransformGroupers/BackwardTransformGrouper.hpp"
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"

namespace QuICC {

   namespace Transform {

      /// Transform configurator selector template
      template <Splitting::Algorithms::Id TAlgo> struct ConfigSelector;

      /// Transform configurator selector specialised for SERIAL case
      template <> struct ConfigSelector<Splitting::Algorithms::SERIAL>
      {
         /// Typedef for forward configurator
         typedef ForwardSerialConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSerialConfigurator  BwdConfigType;
      };

      /// Transform configurator selector specialised for SINGLE1D case
      template <> struct ConfigSelector<Splitting::Algorithms::SINGLE1D>
      {
         /// Typedef for forward configurator
         typedef ForwardSingle1DConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSingle1DConfigurator  BwdConfigType;
      };

      /// Transform configurator selector specialised for SINGLE2D case
      template <> struct ConfigSelector<Splitting::Algorithms::SINGLE2D>
      {
         /// Typedef for forward configurator
         typedef ForwardSingle2DConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSingle2DConfigurator  BwdConfigType;
      };

      /// Transform configurator selector specialised for TUBULAR case
      template <> struct ConfigSelector<Splitting::Algorithms::TUBULAR>
      {
         /// Typedef for forward configurator
         typedef ForwardTubularConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardTubularConfigurator  BwdConfigType;
      };

      /// Transform configurator selector specialised for COUPLED2D case
      template <> struct ConfigSelector<Splitting::Algorithms::COUPLED2D>
      {
         /// Typedef for forward configurator
         typedef ForwardSingle1DConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSingle1DConfigurator  BwdConfigType;
      };

      /// Transform grouper selector template
      template <Splitting::Groupers::Id TGrouper,Splitting::Algorithms::Id TAlgo> struct GrouperSelector;

      /// Transform grouper selector for EQUATION grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::EQUATION,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardEquationGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardEquationGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };

      /// Transform grouper selector for SINGLE1D grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::SINGLE1D,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardSingle1DGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardSingle1DGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };

      /// Transform grouper selector for SINGLE2D grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::SINGLE2D,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardSingle2DGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardSingle2DGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };

      /// Transform grouper selector for TRANSFORM grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::TRANSFORM,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardTransformGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardTransformGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };
   }

   namespace Parallel
   {
      /**
       * @brief Set the transform grouper depending on setup
       */
      void setGrouper(const SplittingDescription& descr, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper);

      /**
       * @brief Templated grouper selector
       */
      template <Splitting::Groupers::Id TGroup, Splitting::Algorithms::Id TAlgo> void setGrouper(Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper);

      /**
       * @brief Templated grouper selector
       */
      template <Splitting::Groupers::Id TGroup> void setGrouper(const Splitting::Algorithms::Id algo, const int dims, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper);

      template <Splitting::Groupers::Id TGroup, Splitting::Algorithms::Id TAlgo> void setGrouper(Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper)
      {
         std::shared_ptr<typename Transform::GrouperSelector<TGroup,TAlgo>::FwdGrouperType> spFwd = std::make_shared<typename Transform::GrouperSelector<TGroup,TAlgo>::FwdGrouperType>();
         std::shared_ptr<typename Transform::GrouperSelector<TGroup,TAlgo>::BwdGrouperType> spBwd = std::make_shared<typename Transform::GrouperSelector<TGroup,TAlgo>::BwdGrouperType>();

         spFwdGrouper = spFwd;
         spBwdGrouper = spBwd;
      }

      template <Splitting::Groupers::Id TGroup> void setGrouper(const Splitting::Algorithms::Id algo, const int dims, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper)
      {
         if(algo == Splitting::Algorithms::SERIAL)
         {
            setGrouper<TGroup,Splitting::Algorithms::SERIAL>(spFwdGrouper, spBwdGrouper);

      #ifdef QUICC_MPI
         #ifdef QUICC_MPIALGO_SINGLE1D
         } else if(algo == Splitting::Algorithms::SINGLE1D)
         {
            setGrouper<TGroup,Splitting::Algorithms::SINGLE1D>(spFwdGrouper, spBwdGrouper);
         #endif //QUICC_MPIALGO_SINGLE1D
         #if defined QUICC_MPIALGO_SINGLE2D
         } else if(dims == 3 && algo == Splitting::Algorithms::SINGLE2D)
         {
            setGrouper<TGroup,Splitting::Algorithms::SINGLE2D>(spFwdGrouper, spBwdGrouper);
         #endif //defined QUICC_MPIALGO_SINGLE2D
         #if defined QUICC_MPIALGO_TUBULAR
         } else if(dims == 3 && algo == Splitting::Algorithms::TUBULAR)
         {
            setGrouper<TGroup,Splitting::Algorithms::TUBULAR>(spFwdGrouper, spBwdGrouper);
         #endif //defined QUICC_MPIALGO_TUBULAR
         #if defined QUICC_MPIALGO_COUPLED2D
         } else if(dims == 3 && algo == Splitting::Algorithms::COUPLED2D)
         {
            setGrouper<TGroup,Splitting::Algorithms::COUPLED2D>(spFwdGrouper, spBwdGrouper);
         #endif //defined QUICC_MPIALGO_COUPLED2D
      #endif //QUICC_MPI
         } else
         {
            throw std::logic_error("Unknown algorithm for transform grouper setup");
         }
      }
   }

}

#endif // QUICC_TRANSFORM_PARALLELSELECTOR_HPP
