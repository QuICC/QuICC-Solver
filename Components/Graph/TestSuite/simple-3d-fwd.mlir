// mod -> phys
// phys -> mod
// Nr x Nphi x Ntheta
!type_uval = !quiccir.view<3x10x6xf64, "R_DCCSC3D_t">
// L x N x M  ( ..., radial, ...)
!type_umod = !quiccir.view<6x2x6xf64, "C_DCCSC3D_t">

!type_tuval = tensor<3x10x6xf64, "R_DCCSC3D_t">
!type_tumod = tensor<6x2x6xf64, "C_DCCSC3D_t">

func.func private @fwd(%Pol: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
  // Pol
  %Pol1 = quiccir.jw.prj %Pol : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
  %Pol1T = quiccir.transpose %Pol1 permutation = [1, 2, 0] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
  %Pol2 = quiccir.al.prj %Pol1T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 2 :i64}
  %Pol2T = quiccir.transpose %Pol2 permutation = [1, 2, 0] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 3 :i64}
  %Pol3 = quiccir.fr.prj %Pol2T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 4 :i64}
  return %Pol3 : tensor<?x?x?xf64>
}

func.func @entry(%thisArr: !llvm.ptr<array<15 x ptr>> {llvm.noalias}, %Rv: !type_uval, %Polv: !type_umod) {
  %Pol = builtin.unrealized_conversion_cast %Polv : !type_umod to !type_tumod
  %Polur = tensor.cast %Pol : !type_tumod to tensor<?x?x?xf64>
  %Rur = call @fwd(%Polur) : (tensor<?x?x?xf64>) -> tensor<?x?x?xf64>
  %R = tensor.cast %Rur : tensor<?x?x?xf64> to !type_tuval
  /// if this is the only consumer write to existing buffer
  quiccir.materialize %R in %Rv : (!type_tuval, !type_uval)
  return
}

// /home/gcastigl/codes/quiccir/build-17/bin/quiccir-opt ../Components/Graph/TestSuite/simple-3d-loop.mlir --inline  --set-quiccir-dims='phys=3,10,6 mods=2,6,6' --set-quiccir-view-lay='lay-ppp2mpp=R_DCCSC3D_t,C_DCCSC3D_t lay-pmp2mmp=C_DCCSC3D_t,C_S1CLCSC3D_t lay-pmm2mmm=C_DCCSC3D_t,C_DCCSC3D_t' --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize
