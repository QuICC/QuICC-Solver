import quicc.recurrence.sphere_jacobi as j

layout = j.w_layout_cpp
# Make alpha and beta the Chebyshev versions
j.chooseParameters(j.w_alpha_chebyshev, j.w_beta_minhalf, layout)


# Higher-order plain derivatives (plus)
j.i3D2p()
j.i3D3p()
j.i4D2p()
j.i4D3p()
j.i4D4p()

# Higher-order Lapl-derivative combos (plus)
j.i3laplDp()
j.i4laplDp()
j.i4laplD2p()

# Higher-order plain derivatives (minus)
j.i3D2m()
j.i3D3m()
j.i4D2m()
j.i4D3m()
j.i4D4m()

# Higher-order Lapl-derivative combos (minus)
j.i3laplDm()
j.i4laplDm()
j.i4laplD2m()
