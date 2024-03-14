
#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

// Passes
#include <Quiccir/IR/QuiccirDialect.h>
#include <Quiccir/Transforms/QuiccirPasses.h>
#include <Quiccir/Pipelines/Passes.h>
#include <Quiccir-c/Utils.h>

#include <mlir/InitAllDialects.h>
#include <mlir/Dialect/Func/Extensions/AllExtensions.h>
#include <mlir/IR/MLIRContext.h>
#include <mlir/Parser/Parser.h>

#include <mlir/Pass/PassManager.h>
// #include <llvm/Support/CommandLine.h>

// JIT
#include <mlir/ExecutionEngine/ExecutionEngine.h>
#include <mlir/ExecutionEngine/OptUtils.h>
#include <llvm/Support/SourceMgr.h>
#include <llvm/Support/TargetSelect.h>
#include <mlir/Target/LLVMIR/Dialect/Builtin/BuiltinToLLVMIRTranslation.h>
#include <mlir/Target/LLVMIR/Dialect/LLVMIR/LLVMToLLVMIRTranslation.h>

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
  // module->dump();

  // Top level (module) pass manager
  mlir::PassManager pm(&ctx);
  mlir::quiccir::quiccLibCallPipelineBuilder(pm);

  // Lower
  if (mlir::failed(pm.run(*module))) {
    CHECK(false);
  }

  // Dump
  // module->dump();

  // JIT

  // Initialize LLVM targets.
  llvm::InitializeNativeTarget();
  llvm::InitializeNativeTargetAsmPrinter();

  // Register the translation from MLIR to LLVM IR, which must happen before we
  // can JIT-compile.
  mlir::registerBuiltinDialectTranslation(*module->getContext());
  mlir::registerLLVMDialectTranslation(*module->getContext());

  // An optimization pipeline to use within the execution engine.
  auto optPipeline = mlir::makeOptimizingTransformer(
      /*optLevel=*/0, /*sizeLevel=*/0,
      /*targetMachine=*/nullptr);

  // Create an MLIR execution engine. The execution engine eagerly JIT-compiles
  // the module.
  mlir::ExecutionEngineOptions engineOptions;
  engineOptions.transformer = optPipeline;
  // engineOptions.sharedLibPaths = executionEngineLibs;
  auto maybeEngine = mlir::ExecutionEngine::create(*module, engineOptions);
  assert(maybeEngine && "failed to construct an execution engine");
  auto &engine = maybeEngine.get();

  // Invoke the JIT-compiled function.
  std::string symbol = "entry";
  auto funSym = engine->lookup(symbol);
  if (auto E = funSym.takeError()) {
    assert(false && "JIT invocation failed");
  }

  // Fake call
  view3_t viewRef_in;
  view3_t viewRef_out;
   std::array<void*, 2> thisArr;
  auto fun = (void (*)(void*, view3_t*, view3_t*))funSym.get();
  fun(thisArr.data(), &viewRef_out, &viewRef_in);

}