/**
 * @file CylinderTransformSteps.hpp
 * @brief Interface for generating transform steps
 */

#ifndef QUICC_TRANSFORM_CYLINDERTRANSFORMSTEPS_HPP
#define QUICC_TRANSFORM_CYLINDERTRANSFORMSTEPS_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/TransformConfigurators/ITransformSteps.hpp"
#include "QuICC/TransformConfigurators/TransformPath.hpp"

namespace QuICC {

   namespace Transform {

      /**
       * @brief Interface for generating transform steps involved in the physical <-> spectral transforms
       */
      class CylinderTransformSteps: public ITransformSteps
      {
         public:
            /**
             * @brief Constructor
             */
            CylinderTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme);

            /**
             * @brief Constructor
             */
            virtual ~CylinderTransformSteps();

            /**
             * @brief Is implementation applicable to scheme?
             */
            static bool applicable(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme);

            /**
             * @brief Generate the list of branches in scalar integration transform
             */
            virtual std::vector<TransformPath>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const;

            /**
             * @brief Generate the list of branches in scalar integration transform
             */
            virtual std::vector<TransformPath>  forwardNLScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const;

            /**
             * @brief Generate the list of branches in vector integration transform
             *
             * @param components Spectral components where to store results: 0: \f$\vec r \nabla\wedge N\f$ 1: \f$\vec r \nabla\wedge\nabla\wedge N\f$
             */
            virtual std::vector<TransformPath>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const;

            /**
             * @brief Generate the list of branches in vector integration transform
             *
             * @param components Spectral components where to store results: 0: \f$\vec r \nabla\wedge N\f$ 1: \f$\vec r \nabla\wedge\nabla\wedge N\f$
             */
            virtual std::vector<TransformPath>  forwardNLVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const;

            /**
             * @brief Generate the list of branches in scalar projection transform
             */
            virtual std::vector<TransformPath>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req) const;

            /**
             * @brief Generate the list of branches in scalar gradient transform
             */
            virtual std::vector<TransformPath>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req) const;

            /**
             * @brief Generate the list of branches in scalar 2nd order gradient transform
             */
            virtual std::vector<TransformPath>  backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req) const;

            /**
             * @brief Generate the list of branches in vector projection transform
             */
            virtual std::vector<TransformPath>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req) const;

            /**
             * @brief Generate the list of branches in vector gradient transform
             */
            virtual std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req) const;

            /**
             * @brief Generate the list of branches in vector curl transform
             */
            virtual std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req) const;

            /**
             * @brief Generate the list of branches in vector divergence transform
             */
            virtual std::vector<TransformPath>  backwardDivergence() const;
      };
   }
}

#endif // QUICC_TRANSFORM_CYLINDERTRANSFORMSTEPS_HPP
