/**
 * @file ImExSBDF2.cpp
 * @brief Implementation of an implicit/explicit SBDF scheme of order 2
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/ImExSBDF2.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   ImExSBDF2::ImExSBDF2()
   {
   }

   ImExSBDF2::~ImExSBDF2()
   {
   }

   // Scheme requires single step
   int ImExSBDF2::steps() const
   {
      return 1;
   }

   // Scheme order
   int ImExSBDF2::order() const
   {
      return 2;
   }

   // Scheme requires field at t_(n-1)
   int ImExSBDF2::fieldMemory() const
   {
      return 1;
   }

   // Scheme requires nonlinear term at t(n-1)
   int ImExSBDF2::nonlinearMemory() const
   {
      return 1;
   }

   // Name of the scheme
   std::string ImExSBDF2::name() const
   {
      return "ImExSBDF2";
   }

   MHDFloat ImExSBDF2::lhsT(const int step) const
   {
      assert(step < this->steps());

      return 3.0/2.0;
   }

   MHDFloat ImExSBDF2::lhsL(const int step) const
   {
      assert(step < this->steps());

      return 1.0;
   }

   MHDFloat ImExSBDF2::rhsT(const int i, const int step) const
   {
      assert(step < this->steps());
      assert(i > -1);
      assert(i < this->fieldMemory()+1);

      MHDFloat coeff;
      if(i == 0)
      {
         coeff = 2.0;
      } else
      {
         coeff = -1.0/2.0;
      }

      return coeff;
   }

   MHDFloat ImExSBDF2::rhsL(const int i, const int step) const
   {
      assert(step < this->steps());
      assert(i > -1);
      assert(i < this->fieldMemory()+1);

      return 0.0;
   }

   MHDFloat ImExSBDF2::rhsN(const int i, const int step) const
   {
      assert(step < this->steps());
      assert(i > -1);
      assert(i < this->nonlinearMemory()+1);

      MHDFloat coeff;
      if(i == 0)
      {
         coeff = 2.0;
      } else
      {
         coeff = -1.0;
      }

      return coeff;
   }

   MHDFloat ImExSBDF2::cEx(const int i) const
   {
      return 1.0;
   }

   int ImExSBDF2::fieldMemory(const int step) const
   {
      assert(step < this->steps());

      return this->fieldMemory();
   }

   int ImExSBDF2::nonlinearMemory(const int step) const
   {
      assert(step < this->steps());

      return this->nonlinearMemory();
   }

   void ImExSBDF2::init()
   {
   }

   bool ImExSBDF2::hasEmbedded() const
   {
      return false;
   }

}
}
