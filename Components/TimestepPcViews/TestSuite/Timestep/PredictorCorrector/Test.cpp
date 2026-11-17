/**
 * @file Test.cpp
 * @brief Source of test
 */

// System includes
//

// Project includes
//
#include "TestSuite/Timestep/PredictorCorrector/Test.hpp"
#include "TestSuite/Io.hpp"
#include "TestSuite/Timestep/PredictorCorrector/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace PredictorCorrector {

Test::Test(const int id) :
    nN(-1),
    basisId(BasisId::WORLAND),
    l(-1),
    dt(-1),
    timeCoeffs(0),
    tEnd(0),
    startRow(1),
    start(0),
    cols(1),
    t0(0.0),
    id(id)
{
   this->configure();
}

void Test::configure()
{
   const std::string refdir = "_refdata/Timestep/PredictorCorrector/";

   switch (this->id)
   {
   case 0:
   case 1:
   case 2:
   case 3:
   case 4:
   case 5:
   case 6:
   case 7:
      this->basisId = BasisId::WORLAND;
      this->equationId = EquationId::DIFFUSION;
      this->fbase = refdir + "Worland/diffusion_id" + std::to_string(this->id);
      break;
   case 10:
   case 11:
   case 12:
   case 13:
   case 14:
   case 15:
   case 16:
   case 17:
      this->basisId = BasisId::WORLAND;
      this->equationId = EquationId::BIDIFFUSION;
      this->fbase =
         refdir + "Worland/bidiffusion_id" + std::to_string(this->id - 10);
      break;
   case 100:
      this->basisId = BasisId::CHEBYSHEV;
      this->equationId = EquationId::DIFFUSION;
      this->fbase =
         refdir + "Chebyshev/diffusion_id" + std::to_string(this->id - 100);
      break;
   default:
      throw std::logic_error("Unknown test ID");
   }

   // Read meta data
   Array meta;
   std::string metaFile = this->fbase + "_meta.dat";
   readList(meta, metaFile);

   this->nN = static_cast<int>(meta(0));
   int s = 0;
   if (this->basisId == BasisId::WORLAND)
   {
      this->l = static_cast<int>(meta(1));
      s = 1;
   }
   this->timeCoeffs.resize(2);
   this->timeCoeffs(0) = meta(1 + s);
   this->timeCoeffs(1) = meta(2 + s);
   this->tEnd = Array::Zero(meta.size() - (3 + s));
   for (int i = 0; i < this->tEnd.size(); i++)
   {
      this->tEnd(i) = std::pow(10.0, meta(3 + s + i));
   }
   this->dt = this->tEnd(0);

   // Other setup
   this->startRow = ArrayI::Zero(1);
}

Matrix Test::getReference()
{
   const std::string refdir = "_refdata/Timestep/PredictorCorrector/";

   std::string refFile = this->fbase + "_ref.dat";

   Matrix ref(this->nN, this->tEnd.size());
   readData(ref, refFile);

   return ref;
}

Matrix Test::getInitial()
{
   const std::string refdir = "_refdata/Timestep/PredictorCorrector/";

   std::string refFile = this->fbase + "_in.dat";

   Matrix ref(this->nN, 1);
   readData(ref, refFile);

   return ref;
}

Matrix Test::getForcing()
{
   const std::string refdir = "_refdata/Timestep/PredictorCorrector/";

   std::string refFile = this->fbase + "_forcing.dat";

   Matrix ref(this->nN, 4);
   readData(ref, refFile);
   ref = -ref;

   return ref;
}

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC
