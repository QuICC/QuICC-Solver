
#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

#include <Quiccir/IR/QuiccirDialect.h>
#include <Quiccir/Transforms/QuiccirPasses.h>
#include <Quiccir/Pipelines/Passes.h>

#include <mlir/InitAllDialects.h>
#include <mlir/Dialect/Func/Extensions/AllExtensions.h>
#include <mlir/IR/MLIRContext.h>
#include <mlir/Parser/Parser.h>

#include <mlir/Pass/PassManager.h>

#include <llvm/Support/CommandLine.h>

// #include "Memory/Cpu\NewDelete.hpp"
// #include "Memory/Cuda/Malloc.hpp"
// #include "Memory/Memory.hpp"
// #include "ViewOps/Pointwise/Cpu/Pointwise.hpp"
// #include "Visitor/Transform/Visitor.hpp"
// #include "ViewOps/ViewMemoryUtils.hpp"
#include "Profiler/Interface.hpp"

int main(int argc, char **argv)
{
  // Register any command line options.
  // mlir::registerAsmPrinterCLOptions();
  // mlir::registerMLIRContextCLOptions();
  // mlir::registerPassManagerCLOptions();

  // llvm::cl::ParseCommandLineOptions(argc, argv, "quiccir jitter\n");

  QuICC::Profiler::Initialize();

  Catch::Session session; // There must be exactly one instance

  auto returnCode = session.run();

  QuICC::Profiler::Finalize();

  return returnCode;
}

TEST_CASE("Simple Tree", "[SimpleTree]")
{
  mlir::DialectRegistry registry;
  mlir::func::registerAllExtensions(registry);
  // Add the following to include *all* MLIR Core dialects, or selectively
  // include what you need like above. You only need to register dialects that
  // will be *parsed* by the tool, not the one generated
  registerAllDialects(registry);

  mlir::MLIRContext ctx(registry);
  // Load our Dialect in this MLIR Context.
  ctx.loadDialect<mlir::quiccir::QuiccirDialect>();

  std::string modStr =
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
    "}\n";

  // Load
  mlir::OwningOpRef<mlir::ModuleOp> module;
  module = mlir::parseSourceString<mlir::ModuleOp>(modStr, &ctx);

  // Dump
  module->dump();

  // Top level (module) pass manager
  mlir::PassManager pm(&ctx);
  mlir::quiccir::quiccLibCallPipelineBuilder(pm);

  // Lower
  if (mlir::failed(pm.run(*module))) {
    CHECK(false);
  }

  // Dump
  module->dump();



}