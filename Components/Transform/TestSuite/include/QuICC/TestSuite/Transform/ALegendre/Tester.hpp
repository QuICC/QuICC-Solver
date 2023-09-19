/**
 * @file Tester.hpp
 * @brief Tester for ALegendre transforms
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_ALEGENDRE_TESTER_HPP
#define QUICC_TESTSUITE_TRANSFORM_ALEGENDRE_TESTER_HPP

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>
#include <sstream>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/TestSuite/Transform/TesterBase.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
//#include "QuICC/Transform/Fft/ALegendre/Setup.hpp"
//#include "QuICC/Transform/Fft/ALegendre/Tools.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"

template<class>
struct sfinae_true : std::true_type{};

namespace internal{
     template<typename C, typename... Args>
          static auto test_transform(int)
                -> sfinae_true<decltype(std::declval<C>().transform(std::declval<Args>()...))>;
       template<typename, typename... Args>
            static auto test_transform(long) -> std::false_type;
}

template <typename T, typename... Args>
class has_transform: public decltype(internal::test_transform<T, Args...>(0)){};

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
         virtual void readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

         /**
          * @brief Read complex data from file
          */
         virtual void readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

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
      int dbRows;
      if(type == TestType::PROJECTOR && ctype == ContentType::INPUT)
      {
         dbRows = spDbSetup->fastSize(0);
      }
      else
      {
         dbRows = data.rows();
      }

      TData dbData = TData::Zero(dbRows, spDbSetup->slowSize());
      std::string fullname = this->makeFilename(dbParam, this->refRoot(), type, ctype);
      readData(dbData, fullname);

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Loop over indexes
      int col = 0;
      const int dataRows = data.rows();
      for(int j = 0; j < spSetup->slowSize(); j++)
      {
         int j_ = spSetup->slow(j);
         // Loop over multiplier
         for(int i = 0; i < spSetup->mult(j); i++)
         {
            data.block(0, col, dataRows, 1) = dbData.block(0, j_, dataRows, 1);
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

      MatrixZ outData;

      for (unsigned int i = 0; i < this->mIter; ++i)
      {
         outData = MatrixZ::Zero(op.outRows(), op.outCols());
         op.transform(outData, inData);
      }

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

      MatrixZ outData;

      bool isReversed = (igrid(igrid.size()-1) < igrid(0));
      if(isReversed)
      {
         this->reverseData(inData);
      }

      for (unsigned int i = 0; i < this->mIter; ++i)
      {
         outData = MatrixZ::Zero(op.outRows(), op.outCols());
         op.transform(outData, inData);
      }

      Matrix out;
      if(param.size() == 1)
      {
         out.resize(2*outData.rows(),outData.cols());
         out.topRows(outData.rows()) = outData.real();
         out.bottomRows(outData.rows()) = outData.imag();
      }
      else
      {
         out = Matrix::Zero(2*spSetup->specSize(),outData.cols());
         out.topRows(outData.rows()) = outData.real();
         out.block(spSetup->specSize(), 0, outData.rows(), outData.cols()) = outData.imag();
      }

      return out;
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyReductor(const ParameterType& param) const
   {
      if constexpr(has_transform<TOp, Matrix&, const MatrixZ&>::value)
      {
         const TestType type = TestType::REDUCTOR;

         // Create setup
         auto spSetup = this->buildSetup(param, type);

         // Input data
         MatrixZ inData(spSetup->specSize(), spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp op;
         internal::Array igrid;
         this->initOperator(op, igrid, spSetup);

         Matrix outData(op.outRows(), op.outCols());

         op.transform(outData, inData);

         return outData;
      }
      else
      {
         throw std::logic_error("This operator is not a reductor");

         Matrix out;
         return out;
      }
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyBFLoop(const ParameterType& param) const
   {
      if constexpr(std::is_same_v<TOp2,void>)
      {
         throw std::logic_error("Backward-forward loop can only be computed if second operator type is given");
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

      // Distributed data meta file
      if(param.size() == 3)
      {
         int np = param.at(1);
         ss << "_np" << np;
         int r = param.at(2);
         ss << "_r" << r;
         ss << "_stage1";
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
      // Read DB metadata
      ParameterType dbParam(param.begin(), param.begin()+1);
      Array dbMeta(0);
      std::string fullname = this->makeFilename(dbParam, this->refRoot(), type, ContentType::META);
      readList(dbMeta, fullname);

      // Read (distributed) metadata
      Array meta(0);
      fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
      readList(meta, fullname);

      // Create setup, three special lines: specN, physN, nModes
      int nMeta = 3;
      if(dbMeta(0) != meta(0) || dbMeta(1) != meta(1))
      {
         throw std::logic_error("Distributed data doesn't match database");
      }
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
