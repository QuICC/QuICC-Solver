/**
 * @file Tester.hpp
 * @brief Tester for Chebyshev linear map transforms
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_CHEBYSHEV_LINEARMAP_TESTER_HPP
#define QUICC_TESTSUITE_TRANSFORM_CHEBYSHEV_LINEARMAP_TESTER_HPP

// System includes
//
#include <catch2/catch.hpp>
#include <set>
#include <sstream>
#include <string>

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/TestSuite/Transform/TesterBase.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"
#include "Types/Typedefs.hpp"

template <typename T, typename... Args> class has_transform
{
   template <typename C,
      typename = decltype(std::declval<C>().transform(std::declval<Args>()...))>
   static std::true_type test(int);
   template <typename C> static std::false_type test(...);

public:
   static constexpr bool value = decltype(test<T>(0))::value;
};

namespace transf = ::QuICC::Transform::Fft::Chebyshev;

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace Chebyshev {

namespace LinearMap {

template <typename TOp, typename TOp2 = void>
class Tester : public Transform::TesterBase<TOp>
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

   /// Typedef for Forward input data
   typedef MatrixZ FwdType;

   /// Typedef for Backward input data
   typedef MatrixZ BwdType;

   /**
    * @brief Build filename extension with resolution information
    */
   virtual std::string resname(const ParameterType& param) const override;

   /**
    * @brief Read data from file
    */
   virtual void readFile(Matrix& data, const ParameterType& param,
      const TestType type, const ContentType ctype) const override;

   /**
    * @brief Read complex data from file
    */
   virtual void readFile(MatrixZ& data, const ParameterType& param,
      const TestType type, const ContentType ctype) const override;

   /**
    * @brief Read data from database file
    */
   template <typename TData>
   void dbReadFile(TData& data, const ParameterType& param, const TestType type,
      const ContentType ctype) const;

   /**
    * @brief Test operator
    */
   virtual Matrix applyOperator(const ParameterType& param,
      const TestType type) const override;

   /**
    * @brief Format the parameters
    */
   virtual std::string formatParameter(
      const ParameterType& param) const override;

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
   template <typename T>
   void initOperator(T& op, const transf::SharedSetup spSetup) const;

   /**
    * @brief Build transform operator setup
    */
   std::shared_ptr<typename TOp::SetupType> buildSetup(
      const ParameterType& param, const TestType type) const;
};

template <typename TOp, typename TOp2>
Tester<TOp, TOp2>::Tester(const std::string& fname, const bool keepData) :
    Transform::TesterBase<TOp>(fname, keepData)
{
   this->appendPath();
}

template <typename TOp, typename TOp2> void Tester<TOp, TOp2>::appendPath()
{
   this->mPath += "Chebyshev/LinearMap/";
}

template <typename TOp, typename TOp2>
void Tester<TOp, TOp2>::readFile(Matrix& data, const ParameterType& param,
   const TestType type, const ContentType ctype) const
{
   if (param.size() == 1)
   {
      TesterBase<TOp>::basicReadFile(data, param, type, ctype);
   }
   else
   {
      this->dbReadFile(data, param, type, ctype);
   }
}

template <typename TOp, typename TOp2>
void Tester<TOp, TOp2>::readFile(MatrixZ& data, const ParameterType& param,
   const TestType type, const ContentType ctype) const
{
   if (param.size() == 1)
   {
      TesterBase<TOp>::basicReadFile(data, param, type, ctype);
   }
   else
   {
      this->dbReadFile(data, param, type, ctype);
   }
}

template <typename TOp, typename TOp2>
template <typename TData>
void Tester<TOp, TOp2>::dbReadFile(TData& data, const ParameterType& param,
   const TestType type, const ContentType ctype) const
{
   // Read database file
   ParameterType dbParam = {param.at(0)};
   auto spDbSetup = this->buildSetup(dbParam, type);
   int dbRows = data.rows();
   int dbCols = spDbSetup->slowSize();
   if (type == TestType::PROJECTOR && ctype == ContentType::INPUT)
   {
      dbRows = spDbSetup->fastSize(0);
   }

   // Read database file
   TData dbData = TData::Zero(dbRows, dbCols);
   std::string fullname =
      this->makeFilename(dbParam, this->refRoot(), type, ctype);
   readData(dbData, fullname);

   // Create setup
   auto spSetup = this->buildSetup(param, type);

   // Count modes
   int nModes = 0;
   for (int j = 0; j < spSetup->slowSize(); j++)
   {
      nModes += spSetup->mult(j);
   }

   std::function<void(TData&, const TData&, const int, const int)> fillData =
      [](TData& data, const TData& db, const int idx, const int j_)
   {
      const int dataRows = data.rows();
      data.block(0, idx, dataRows, 1) = db.block(0, j_, dataRows, 1);
   };

   // Special case for energy reduction
   if (type == TestType::REDUCTOR && ctype == ContentType::REFERENCE &&
       data.rows() == nModes && data.cols() == 1)
   {
      fillData = [](TData& data, const TData& db, const int idx, const int j_)
      { data(idx, 0) = db(0, j_); };
   }

   // Loop over indexes
   int idx = 0;
   for (int j = 0; j < spSetup->slowSize(); j++)
   {
      int j_ = spSetup->slow(j);
      // Loop over multiplier
      for (int i = 0; i < spSetup->mult(j); i++)
      {
         fillData(data, dbData, idx, j_);
         idx++;
      }
   }
}

template <typename TOp, typename TOp2>
Matrix Tester<TOp, TOp2>::applyOperator(const ParameterType& param,
   const TestType type) const
{
   Matrix outData;
   switch (type)
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

template <typename TOp, typename TOp2>
Matrix Tester<TOp, TOp2>::applyProjector(const ParameterType& param) const
{
   if constexpr (has_transform<TOp, FwdType&, const BwdType&>::value)
   {
      const TestType type = TestType::PROJECTOR;

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Input data
      BwdType inData(spSetup->specSize(), spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp op;
      this->initOperator(op, spSetup);

      FwdType outData(op.outRows(), op.outCols());

      op.transform(outData, inData);

      if constexpr (std::is_same_v<FwdType, MatrixZ>)
      {
         Matrix out(2 * outData.rows(), outData.cols());
         out.topRows(outData.rows()) = outData.real();
         out.bottomRows(outData.rows()) = outData.imag();

         return out;
      }
      else
      {
         return outData;
      }
   }
   else
   {
      throw std::logic_error("This operator is not an projector");

      Matrix out;
      return out;
   }
}

template <typename TOp, typename TOp2>
Matrix Tester<TOp, TOp2>::applyIntegrator(const ParameterType& param) const
{
   if constexpr (has_transform<TOp, BwdType&, const FwdType&>::value)
   {
      const TestType type = TestType::INTEGRATOR;

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Input data
      FwdType inData(spSetup->fwdSize(), spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp op;
      this->initOperator(op, spSetup);

      BwdType outData(op.outRows(), op.outCols());

      op.transform(outData, inData);

      Matrix out(2 * outData.rows(), outData.cols());
      out.topRows(outData.rows()) = outData.real();
      out.bottomRows(outData.rows()) = outData.imag();

      return out;
   }
   else
   {
      throw std::logic_error("This operator is not an integrator");

      Matrix out;
      return out;
   }
}

template <typename TOp, typename TOp2>
Matrix Tester<TOp, TOp2>::applyReductor(const ParameterType& param) const
{
//   if constexpr (has_transform<TOp, Matrix&, const BwdType&>::value)
//   {
//      const TestType type = TestType::REDUCTOR;
//
//      // Create setup
//      auto spSetup = this->buildSetup(param, type);
//
//      // Input data
//      BwdType inData(spSetup->specSize(), spSetup->blockSize());
//      this->readFile(inData, param, type, ContentType::INPUT);
//
//      TOp op;
//      this->initOperator(op, spSetup);
//
//      Matrix outData(op.outRows(), op.outCols());
//
//      op.transform(outData, inData);
//
//      return outData;
//   }
//   else
//   {
      throw std::logic_error("This operator is not an reductor");

      Matrix out;
      return out;
//   }
}

template <typename TOp, typename TOp2>
Matrix Tester<TOp, TOp2>::applyBFLoop(const ParameterType& param) const
{
   if constexpr (std::is_same_v<TOp2, void>)
   {
      throw std::logic_error("Bacward-forward loop can only be computed if "
                             "second operator type is given");
   }
   else
   {
      const TestType type = TestType::BFLOOP;

      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Input data
      BwdType inData(spSetup->specSize(), spSetup->blockSize());
      this->readFile(inData, param, type, ContentType::INPUT);

      TOp opB;
      this->initOperator(opB, spSetup);

      FwdType tmpData(opB.outRows(), opB.outCols());

      opB.transform(tmpData, inData);

      TOp2 opF;
      this->initOperator(opF, spSetup);

      BwdType outData(opF.outRows(), opF.outCols());

      opF.transform(outData, tmpData);

      Matrix out(2 * outData.rows(), outData.cols());
      out.topRows(outData.rows()) = outData.real();
      out.bottomRows(outData.rows()) = outData.imag();

      return out;
   }
}

template <typename TOp, typename TOp2>
std::string Tester<TOp, TOp2>::resname(const ParameterType& param) const
{
   auto id = param.at(0);

   std::stringstream ss;
   ss.precision(10);
   ss << "_id" << id;

   if (param.size() == 3)
   {
      int np = param.at(1);
      ss << "_np" << np;
      int r = param.at(2);
      ss << "_r" << r;
   }

   return ss.str();
}

template <typename TOp, typename TOp2>
std::string Tester<TOp, TOp2>::formatParameter(const ParameterType& param) const
{
   auto id = param.at(0);

   std::stringstream ss;
   ss << "id: " << id;

   if (param.size() == 3)
   {
      int np = param.at(1);
      ss << ", np: " << np;
      int r = param.at(2);
      ss << ", r: " << r;
   }

   return ss.str();
}

template <typename TOp, typename TOp2>
template <typename T>
void Tester<TOp, TOp2>::initOperator(T& op,
   const transf::SharedSetup spSetup) const
{
   op.init(spSetup);
}

template <typename TOp, typename TOp2>
std::shared_ptr<typename TOp::SetupType> Tester<TOp, TOp2>::buildSetup(
   const ParameterType& param, const TestType type) const
{
   // Read DB metadata
   ParameterType dbParam(param.begin(), param.begin() + 1);
   Array dbMeta(0);
   std::string fullname =
      this->makeFilename(dbParam, this->refRoot(), type, ContentType::META);
   readList(dbMeta, fullname);

   // Read (distributed) metadata
   Array meta(0);
   fullname =
      this->makeFilename(param, this->refRoot(), type, ContentType::META);
   readList(meta, fullname);

   // Create setup
   int nMeta = 5;
   if (dbMeta(0) != meta(0) || dbMeta(1) != meta(1) || dbMeta(3) != meta(3) ||
       dbMeta(4) != meta(4))
   {
      throw std::logic_error("Distributed data doesn't match database");
   }
   int specN = meta(0);
   int physN = meta(1);
   double lb = meta(2);
   double ub = meta(3);
   int nModes = meta(4);

   // Gather indices
   std::map<int, std::pair<int, int>> indices;
   assert((meta.size() - nMeta - 2 * nModes) % 2 == 0);
   int nModes2D = 0;
   int h = nMeta;

   // Create mode list
   for (int i = 0; i < nModes; i++)
   {
      int k_ = static_cast<int>(meta(h));
      int mult = static_cast<int>(meta(h + 1));
      indices.insert(std::pair(k_, std::make_pair(mult, 0)));
      h += 2;
      nModes2D += mult;
   }

   auto spSetup = std::make_shared<typename TOp::SetupType>(physN, nModes2D,
      specN, GridPurpose::SIMULATION);
   spSetup->setBoxScale(1.0);
   spSetup->setBounds(lb, ub);

   // Check meta data size
   if (meta.size() - nMeta - 2 * nModes - 2 * nModes2D != 0)
   {
      throw std::logic_error(
         "Meta data format is not supported (file: " + fullname + ")");
   }

   // Set truncation
   for (auto& [k_, p]: indices)
   {
      // Get 1D truncation of first 2D mode
      p.second = meta(h + 1);

      // Check all 2D modes have same truncation
      for (int j = 0; j < p.first; j++)
      {
         if (p.second != meta(h + 1))
         {
            throw std::logic_error(
               "Meta data format is not supported (file: " + fullname + ")");
         }
         h += 2;
      }
   }

   // Add indices with multiplier
   for (const auto& [k_, p]: indices)
   {
      spSetup->addIndex(k_, p.first);
   }
   spSetup->lock();

   return spSetup;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_TRANSFORM_CHEBYSHEV_LINEARMAP_TESTER_HPP
