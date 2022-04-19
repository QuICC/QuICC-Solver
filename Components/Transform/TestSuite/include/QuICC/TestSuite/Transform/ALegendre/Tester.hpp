/**
 * @file Tester.hpp
 * @brief Tester for ALegendre transforms
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_ALEGENDRE_TESTER_HPP
#define QUICC_TESTSUITE_TRANSFORM_ALEGENDRE_TESTER_HPP


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>
#include <sstream>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/TestSuite/Transform/TesterBase.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
//#include "QuICC/Transform/Fft/ALegendre/Setup.hpp"
//#include "QuICC/Transform/Fft/ALegendre/Tools.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"

namespace transf = ::QuICC::Transform;

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace ALegendre {

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
          * @brief Read real data from file
          */
         virtual void readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

         /**
          * @brief Read complex data from file
          */
         virtual void readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

         /**
          * @brief Read data from database file
          */
         template <typename TData> void dbReadFile(TData& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

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
         template <typename T> void initOperator(T& op, internal::Array& igrid, const transf::Poly::ALegendre::SharedSetup spSetup) const;

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
      this->mPath += "ALegendre/";
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      if(param.size() == 1)
      {
         TesterBase<TOp>::basicReadFile(data, param, type, ctype);
      }
      else
      {
         this->dbReadFile(data, param, type, ctype);
      }
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      if(param.size() == 1)
      {
         TesterBase<TOp>::basicReadFile(data, param, type, ctype);
      }
      else
      {
         this->dbReadFile(data, param, type, ctype);
      }
   }

   template <typename TOp, typename TOp2> template <typename TData> void Tester<TOp,TOp2>::dbReadFile(TData& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      // Read database file
      ParameterType dbParam = {param.at(0)};
      auto spDbSetup = this->buildSetup(dbParam, type);
      TData dbData = TData::Zero(data.rows(), spDbSetup->slowSize());
      std::string fullname = this->makeFilename(dbParam, this->refRoot(), type, ctype);
      readData(dbData, fullname);

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Loop over indexes
      int col = 0;
      for(int j = 0; j < spSetup->slowSize(); j++)
      {
         int j_ = spSetup->slow(j);
         // Loop over multiplier
         for(int i = 0; i < spSetup->mult(j); i++)
         {
            data.col(col) = dbData.col(j_);
            col++;
         }
      }
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
      MatrixZ inData = MatrixZ::Zero(spSetup->fastSize(0),spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp op;
      internal::Array igrid;
      this->initOperator(op, igrid, spSetup);

      MatrixZ outData = MatrixZ::Zero(op.outRows(), op.outCols());

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
      MatrixZ inData = MatrixZ::Zero(spSetup->fwdSize(),spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp op;
      internal::Array igrid;
      this->initOperator(op, igrid, spSetup);

      MatrixZ outData = MatrixZ::Zero(op.outRows(), op.outCols());

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
         MatrixZ inData = MatrixZ::Zero(spSetup->fastSize(0),spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp opB;
         internal::Array igrid;
         this->initOperator(opB, igrid, spSetup);

         MatrixZ tmpData = MatrixZ::Zero(opB.outRows(), opB.outCols());

         opB.transform(tmpData, inData);

         TOp2 opF;
         this->initOperator(opF, igrid, spSetup);

         MatrixZ outData = MatrixZ::Zero(opF.outRows(), opF.outCols());

         opF.transform(outData, tmpData);

         Matrix out(2*outData.rows(),outData.cols());
         out.topRows(outData.rows()) = outData.real();
         out.bottomRows(outData.rows()) = outData.imag();

         return out;
      }
   }

   template <typename TOp, typename TOp2> std::string Tester<TOp,TOp2>::resname(const ParameterType& param) const
   {
      auto id = param.at(0);

      std::stringstream ss;
      ss.precision(10);
      ss << "_id" << id;

      if(param.size() == 3)
      {
         int np = param.at(1);
         ss << "_np" << np;
         int r = param.at(2);
         ss << "_r" << r;
      }

      return ss.str();
   }

   template <typename TOp, typename TOp2> std::string Tester<TOp,TOp2>::formatParameter(const ParameterType& param) const
   {
      auto id = param.at(0);

      std::stringstream ss;
      ss << "id: " << id;

      if(param.size() == 3)
      {
         int np = param.at(1);
         ss << ", np: " << np;
         int r = param.at(2);
         ss << ", r: " << r;
      }

      return ss.str();
   }

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::initOperator(T& op, internal::Array& igrid, const transf::Poly::ALegendre::SharedSetup spSetup) const
   {
      // Create quadrature
      internal::Array iweights;
      ::QuICC::Polynomial::Quadrature::LegendreRule quad;
      quad.computeQuadrature(igrid, iweights, spSetup->fwdSize());

      op.init(spSetup, igrid, iweights);
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
      if(param.size() == 1)
      {
         for(int i = nMeta; i < meta.size(); i++)
         {
            int k_ = static_cast<int>(meta(i));
            indices[k_]++;
         }
      }
      else
      {
         assert((meta.size() - nMeta) % 2 == 0);
         for(int i = nMeta; i < meta.size(); i += 2)
         {
            int k_ = static_cast<int>(meta(i));
            int mult = static_cast<int>(meta(i+1));
            indices.insert(std::pair(k_,mult));
         }
      }

      // Add indices with multiplier
      for(const auto& [k_, mult]: indices)
      {
         spSetup->addIndex(k_, mult);
      }
      spSetup->lock();

      return spSetup;
   }

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::reverseData(T& data) const
   {
      data.colwise().reverseInPlace();
   }

} // ALegendre
} // Transform
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_TRANSFORM_ALEGENDRE_TESTER_HPP
