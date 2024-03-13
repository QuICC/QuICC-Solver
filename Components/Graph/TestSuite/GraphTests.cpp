
#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

#include <Quiccir-c/Dialects.h>
#include <mlir-c/IR.h>
#include <mlir-c/RegisterEverything.h>

// #include "Memory/Cpu\NewDelete.hpp"
// #include "Memory/Cuda/Malloc.hpp"
// #include "Memory/Memory.hpp"
// #include "Qir/Qir.hpp"
// #include "QirBuilder/Builder.hpp"
// #include "ViewOps/Pointwise/Cpu/Pointwise.hpp"
// #include "Visitor/Transform/Visitor.hpp"
// #include "ViewOps/ViewMemoryUtils.hpp"
#include "Profiler/Interface.hpp"

int main()
{
   QuICC::Profiler::Initialize();

   Catch::Session session; // There must be exactly one instance

   auto returnCode = session.run();

   QuICC::Profiler::Finalize();

   return returnCode;
}

/// @brief \todo only register what's needed
/// @param ctx
static void registerAllUpstreamDialects(MlirContext ctx) {
  MlirDialectRegistry registry = mlirDialectRegistryCreate();
  mlirRegisterAllDialects(registry);
  mlirContextAppendDialectRegistry(ctx, registry);
  mlirDialectRegistryDestroy(registry);
}

TEST_CASE("Simple Tree", "[SimpleTree]")
{
  MlirContext ctx = mlirContextCreate();
  registerAllUpstreamDialects(ctx);

  mlirDialectHandleRegisterDialect(mlirGetDialectHandle__quiccir__(), ctx);

  MlirModule module = mlirModuleCreateParse(
      ctx, mlirStringRefCreateFromCString(
        "!type_umod = !quiccir.view<1x2x3xf64, \"layoutUmod\">\n"
        "!type_uval = !quiccir.view<1x3x3xf64, \"layoutUval\">\n"
        "!type_tumod = tensor<1x2x3xf64, \"layoutUmod\">\n"
        "!type_tuval = tensor<1x3x3xf64, \"layoutUval\">\n"
        "func.func @entry(%thisArr: !llvm.ptr<array<2 x ptr>> {llvm.noalias}, %uout: !type_umod, %umod: !type_umod) {\n"
        "  %tumod = builtin.unrealized_conversion_cast %umod : !type_umod to !type_tumod\n"
        "  %tuval = quiccir.jw.prj %tumod : !type_tumod -> !type_tuval attributes{implptr = 0 :i64}\n"
        "  %ret = quiccir.jw.int %tuval : !type_tuval -> !type_tumod attributes{implptr = 1 :i64}\n"
        "  quiccir.materialize %ret in %uout : (!type_tumod, !type_umod)\n"
        "  return\n"
        "}\n"));
  if (mlirModuleIsNull(module)) {
    printf("ERROR: Could not parse.\n");
    mlirContextDestroy(ctx);
    CHECK(false);
  }
  else {
    MlirOperation moduleOp = mlirModuleGetOperation(module);
    mlirOperationDump(moduleOp);

    mlirModuleDestroy(module);
    mlirContextDestroy(ctx);
  }

}