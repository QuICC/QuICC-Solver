// phys -> mod
// Nr x Nphi x Ntheta
// !type_uval = !quiccir.view<3x10x6xf64, "R_DCCSC3D_t">
// L x N x M  ( ..., radial, ...)
// !type_umod = !quiccir.view<6x2x6xf64, "C_DCCSC3D_t">

func.func @entry(%R: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
  // R
  %R1 = quiccir.fr.int %R : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 10 :i64}
  %R1T = quiccir.transpose %R1 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 11 :i64}
  %R2 = quiccir.al.int %R1T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 12 :i64}
  %R2T = quiccir.transpose %R2 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 13 :i64}
  %R3 = quiccir.jw.int %R2T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 14 :i64}
  return %R3 : tensor<?x?x?xf64>
}

// ./bin/quiccir-opt ../simple-3d-fwd.mlir --inline --quiccir-view-wrapper='dim-rets=7,3,6 dim-args=7,3,6 lay-args=C_DCCSC3D_t lay-rets=C_DCCSC3D_t' --inline --set-quiccir-dims='phys=6,10,10 mods=3,6,7' --set-quiccir-view-lay='lay-ppp2mpp=R_DCCSC3D_t,C_DCCSC3D_t lay-pmp2mmp=C_DCCSC3D_t,C_S1CLCSC3D_t lay-pmm2mmm=C_DCCSC3D_t,C_DCCSC3D_t' --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize
