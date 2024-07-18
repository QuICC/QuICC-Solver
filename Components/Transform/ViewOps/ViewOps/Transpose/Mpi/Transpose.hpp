/**
 * @file Transpose.hpp
 * @brief Transpose operations on Views
 */
#pragma once

// External includes
//

// Project includes
//
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "View/View.hpp"
#include "ViewOps/Transpose/Mpi/Comm.hpp"
#include "ViewOps/Transpose/Mpi/Coordinates.hpp"
#include "ViewOps/Transpose/Tags.hpp"


namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Mpi backends
namespace Mpi {

using namespace QuICC::Operator;

/// @brief Transpose operator
/// @tparam Tout
/// @tparam Tin
template <class Tout, class Tin, class Perm>
class Op : public UnaryBaseOp<Op<Tout, Tin, Perm>, Tout, Tin>
{
public:
   /// @brief default constructor
   Op(std::shared_ptr<Comm<typename Tin::ScalarType>> comm) : _comm(comm){};
   Op() = delete;
   /// @brief dtor
   ~Op() = default;

private:
   /// @brief action implementation
   /// @param out output View
   /// @param in input View
   void applyImpl(Tout& out, const Tin& in);
   /// @brief give access to base class
   friend UnaryBaseOp<Op<Tout, Tin, Perm>, Tout, Tin>;
   /// @brief communicator object
   std::shared_ptr<Comm<typename Tin::ScalarType>> _comm;
};


template <class Tout, class Tin, class Perm>
void Op<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::applyImpl");

   auto comm = _comm.get();
   assert(comm != nullptr);

   if (!comm->isSetup())
   {
      // Get Coo
      std::vector<point_t> cooOld = View::getCoo<Tin, p012_t>(in);
      std::vector<point_t> cooNew = View::getCoo<Tout, Perm>(out);
      assert(cooOld.size() == in.size());
      assert(cooNew.size() == out.size());
      // Setup
      comm->setComm(cooNew, cooOld);
   }

   comm->exchange(out.data(), in.data());
}

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
