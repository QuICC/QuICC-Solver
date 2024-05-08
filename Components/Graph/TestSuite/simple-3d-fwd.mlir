// mod -> phys
// phys -> mod
// Nr x Nphi x Ntheta
!type_uval = !quiccir.view<3x10x6xf64, "R_DCCSC3D_t">
// L x N x M  ( ..., radial, ...)
!type_umod = !quiccir.view<6x2x6xf64, "C_DCCSC3D_t">

!type_tuval = tensor<3x10x6xf64, "R_DCCSC3D_t">
!type_tumod = tensor<6x2x6xf64, "C_DCCSC3D_t">

func.func private @fwd(%R: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
  // R
  %R1 = quiccir.fr.int %R : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 10 :i64}
  %R1T = quiccir.transpose %R1 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 11 :i64}
  %R2 = quiccir.al.int %R1T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 12 :i64}
  %R2T = quiccir.transpose %R2 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 13 :i64}
  %R3 = quiccir.jw.int %R2T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 14 :i64}
  return %R3 : tensor<?x?x?xf64>
}

func.func @entry(%thisArr: !llvm.ptr<array<15 x ptr>> {llvm.noalias}, %Polv: !type_umod, %Rv: !type_uval) {
  %R = builtin.unrealized_conversion_cast %Rv : !type_uval to !type_tuval
  %Rur = tensor.cast %R : !type_tuval to tensor<?x?x?xf64>
  %Polur = call @fwd(%Rur) : (tensor<?x?x?xf64>) -> tensor<?x?x?xf64>
  %Pol = tensor.cast %Polur : tensor<?x?x?xf64> to !type_tumod
  /// if this is the only consumer write to existing buffer
  quiccir.materialize %Pol in %Polv : (!type_tumod, !type_umod)
  return
}

// /home/gcastigl/codes/quiccir/build-17/bin/quiccir-opt ../Components/Graph/TestSuite/simple-3d-loop.mlir --inline  --set-quiccir-dims='phys=3,10,6 mods=2,6,6' --set-quiccir-view-lay='lay-ppp2mpp=R_DCCSC3D_t,C_DCCSC3D_t lay-pmp2mmp=C_DCCSC3D_t,C_S1CLCSC3D_t lay-pmm2mmm=C_DCCSC3D_t,C_DCCSC3D_t' --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize
