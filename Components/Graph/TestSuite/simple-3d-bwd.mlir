// mod -> phys
// Nr x Nphi x Ntheta
// !type_uval = !quiccir.view<3x10x6xf64, "R_DCCSC3D_t">
// L x N x M  ( ..., radial, ...)
// !type_umod = !quiccir.view<6x2x6xf64, "C_DCCSC3D_t">


func.func @entry(%Pol: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
  // Pol
  %Pol1 = quiccir.jw.prj %Pol : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
  %Pol1T = quiccir.transpose %Pol1 permutation = [1, 2, 0] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
  %Pol2 = quiccir.al.prj %Pol1T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 2 :i64}
  %Pol2T = quiccir.transpose %Pol2 permutation = [1, 2, 0] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 3 :i64}
  %Pol3 = quiccir.fr.prj %Pol2T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 4 :i64}
  return %Pol3 : tensor<?x?x?xf64>
}

// ./bin/quiccir-opt ../simple-3d-bwd.mlir --inline --quiccir-view-wrapper='dim-rets=7,3,6 dim-args=7,3,6 lay-args=C_DCCSC3D_t lay-rets=C_DCCSC3D_t' --inline --set-quiccir-dims='phys=6,10,10 mods=3,6,7' --set-quiccir-view-lay='lay-ppp2mpp=R_DCCSC3D_t,C_DCCSC3D_t lay-pmp2mmp=C_DCCSC3D_t,C_S1CLCSC3D_t lay-pmm2mmm=C_DCCSC3D_t,C_DCCSC3D_t' --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize