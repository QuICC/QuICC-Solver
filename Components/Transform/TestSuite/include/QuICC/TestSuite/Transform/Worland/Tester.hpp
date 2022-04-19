/**
 * @file Tester.hpp
 * @brief Tester for Worland transforms
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_WORLAND_TESTER_HPP
#define QUICC_TESTSUITE_TRANSFORM_WORLAND_TESTER_HPP


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/TestSuite/Transform/TesterBase.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"
#include "QuICC/Transform/Fft/Worland/Tools.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"

namespace transf = ::QuICC::Transform;

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace Worland {

   template <typename TOp, typename TOp2 = void> class Tester: public Transform::TesterBase<TOp>
   {
      public:
         /// Typedef for parameter type
         typedef typename Transform::TesterBase<TOp>::ParameterType ParameterType;

         /**
          * @brief Constructor
          */
         Tester(const std::string& fname, const bool keepData);

         /*
          * @brief Destructor
          */
         virtual ~Tester() = default;

      protected:
         /// Typedef ContentType from base
         typedef typename Transform::TesterBase<TOp>::ContentType ContentType;

         /**
          * @brief Build filename extension with resolution information
          */
         virtual std::string resname(const ParameterType& param) const override;

         /**
          * @brief Test operator
          */
         virtual Matrix applyOperator(const ParameterType& param, const TestType type) const override;

         /**
          * @brief Format the parameters
          */
         virtual std::string formatParameter(const ParameterType& param) const override;

      private:
         /**
          * @brief Append specific path
          */
         void appendPath();

         /**
          * @brief Test projector
          */
         virtual Matrix applyProjector(const ParameterType& param) const;

         /**
          * @brief Test integrator
          */
         virtual Matrix applyIntegrator(const ParameterType& param) const;

         /**
          * @brief Test reductor
          */
         virtual Matrix applyReductor(const ParameterType& param) const;

         /**
          * @brief Test backward-forward loop
          */
         virtual Matrix applyBFLoop(const ParameterType& param) const;

         /**
          * @brief Initialize Poly operator
          */
         template <typename T> void initOperator(T& op, internal::Array& igrid, const transf::Poly::SharedSetup spSetup) const;

         /**
          * @brief Initialize FFT operator
          */
         template <typename T> void initOperator(T& op, internal::Array& igrid, const transf::Fft::Worland::SharedSetup spSetup) const;

         /**
          * @brief Build transform operator setup
          */
         std::shared_ptr<typename TOp::SetupType> buildSetup(const ParameterType& param, const TestType type) const;

         /**
          * @brief Reverse data order
          */
         template <typename T> void reverseData(T& data) const;

   };

   template <typename TOp, typename TOp2> Tester<TOp,TOp2>::Tester(const std::string& fname, const bool keepData)
      : Transform::TesterBase<TOp>(fname, keepData)
   {
      this->appendPath();
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::appendPath()
   {
      this->mPath += "Worland/";
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyOperator(const ParameterType& param, const TestType type) const
   {
      Matrix outData;
      switch(type)
      {
         case TestType::PROJECTOR:
            outData = this->applyProjector(param);
            break;
         case TestType::INTEGRATOR:
            outData = this->applyIntegrator(param);
            break;
         case TestType::REDUCTOR:
            outData = this->applyReductor(param);
            break;
         case TestType::BFLOOP:
            outData = this->applyBFLoop(param);
            break;
         default:
            throw std::logic_error("Test type not implemented");
            break;
      }

      return outData;
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyProjector(const ParameterType& param) const
   {
      const TestType type = TestType::PROJECTOR;

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Input data
      MatrixZ inData(spSetup->specSize(),spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp op;
      internal::Array igrid;
      this->initOperator(op, igrid, spSetup);

      MatrixZ outData(op.outRows(), op.outCols());

      op.transform(outData, inData);
      bool isReversed = (igrid(igrid.size()-1) < igrid(0));
      if(isReversed)
      {
         this->reverseData(outData);
      }

      Matrix out(2*outData.rows(),outData.cols());
      out.topRows(outData.rows()) = outData.real();
      out.bottomRows(outData.rows()) = outData.imag();

      return out;
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyIntegrator(const ParameterType& param) const
   {
      const TestType type = TestType::INTEGRATOR;

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Input data
      MatrixZ inData(spSetup->fwdSize(),spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp op;
      internal::Array igrid;
      this->initOperator(op, igrid, spSetup);

      MatrixZ outData(op.outRows(), op.outCols());

      bool isReversed = (igrid(igrid.size()-1) < igrid(0));
      if(isReversed)
      {
         this->reverseData(inData);
      }

      op.transform(outData, inData);

      Matrix out(2*outData.rows(),outData.cols());
      out.topRows(outData.rows()) = outData.real();
      out.bottomRows(outData.rows()) = outData.imag();

      return out;
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyReductor(const ParameterType& param) const
   {
      const TestType type = TestType::REDUCTOR;

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Input data
      Matrix inData(spSetup->specSize(), spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp op;
      internal::Array igrid;
      this->initOperator(op, igrid, spSetup);

      Matrix outData(op.outRows(), op.outCols());

      op.transform(outData, inData);

      return outData;
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyBFLoop(const ParameterType& param) const
   {
      if constexpr(std::is_same_v<TOp2,void>)
      {
         throw std::logic_error("Bacward-forward loop can only be computed if second operator type is given");
      }
      else
      {
         const TestType type = TestType::BFLOOP;

         // Create setup
         auto spSetup = this->buildSetup(param, type);

         // Input data
         MatrixZ inData(spSetup->specSize(),spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp opB;
         internal::Array igrid;
         this->initOperator(opB, igrid, spSetup);

         MatrixZ tmpData(opB.outRows(), opB.outCols());

         opB.transform(tmpData, inData);

         TOp2 opF;
         this->initOperator(opF, igrid, spSetup);

         MatrixZ outData(opF.outRows(), opF.outCols());

         opF.transform(outData, tmpData);

         Matrix out(2*outData.rows(),outData.cols());
         out.topRows(outData.rows()) = outData.real();
         out.bottomRows(outData.rows()) = outData.imag();

         return out;
      }
   }

   template <typename TOp, typename TOp2> std::string Tester<TOp,TOp2>::resname(const ParameterType& param) const
   {
      int id = static_cast<int>(param.at(0));

      std::string s = "_id" + std::to_string(id);
      return s;
   }

   template <typename TOp, typename TOp2> std::string Tester<TOp,TOp2>::formatParameter(const ParameterType& param) const
   {
      int id = static_cast<int>(param.at(0));

      std::stringstream ss;
      ss << "id: " << id;

      return ss.str();
   }

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::initOperator(T& op, internal::Array& igrid, const transf::Poly::SharedSetup spSetup) const
   {
      // Create quadrature
      internal::Array iweights;
      ::QuICC::Polynomial::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, spSetup->fwdSize());

      op.init(spSetup, igrid, iweights);
   }

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::initOperator(T& op, internal::Array& igrid, const transf::Fft::Worland::SharedSetup spSetup) const
   {
      ::QuICC::Transform::Fft::Worland::Tools::computeGrid(igrid, spSetup->fwdSize());

      op.init(spSetup);
   }

   template <typename TOp, typename TOp2> std::shared_ptr<typename TOp::SetupType> Tester<TOp,TOp2>::buildSetup(const ParameterType& param, const TestType type) const
   {
      // Read metadata
      Array meta(0);
      std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
      readList(meta, fullname);

      // Create setup
      int nMeta = 2;
      int specN = meta(0);
      int physN = meta(1);
      auto spSetup = std::make_shared<typename TOp::SetupType>(physN, specN, GridPurpose::SIMULATION);

      // Gather indices
      std::map<int,int> indices;
      for(int i = nMeta; i < meta.size(); i++)
      {
         int l = static_cast<int>(meta(i));
         indices[l]++;
      }

      // Add indices with multiplier
      for(const auto& [l, mult]: indices)
      {
         spSetup->addIndex(l, mult);
      }
      spSetup->lock();

      return spSetup;
   }

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::reverseData(T& data) const
   {
      data.colwise().reverseInPlace();
   }

}
}
}
}

#endif //QUICC_TESTSUITE_TRANSFORM_WORLAND_TESTER_HPP
