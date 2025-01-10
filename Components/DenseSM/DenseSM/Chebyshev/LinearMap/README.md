# Dense operators

These are linear operators that can include a generic function of radius in their definition. These operators perform a projection from one spectral field to another spectral field.


Consider a field $v$ whose spectral coefficients are $\hat{\bf v} = (\hat{v}_0, \hat{v}_1, \dots, \hat{v}_{N-1})^T$. 

As an example, consider, the ```R2DivR2FD1R1``` operator. In the name ```Div``` here stands for 'division by' and ```F``` indicates a function of radius (a scalar field). This operator is meant to perform the following operation:

$$
\frac{r^2}{r^2} f(r) \partial_r(r v) = f(r) \partial_r(r v).
$$

The ```R2DivR2``` is obviously redundant, but it is there to keep better track of where the operator acts in the original equations.

Considering a Chebyshev expansion for $v$, the above is:
$$
f \partial_r(r v) = \sum_n t_n \hat{v}_n f \partial_r\left(r T_n\right)
$$
where $t_n$ is a factor that takes care of the different normalisation of $T_0$ with respect of $T_n$ with $n>0$.



Spectrally, this operator converts the $\hat{\bf v}$ into a set of $\hat{\bf w}$ so that:

$$
\sum_n t_n \hat{v}_n f(r) \partial_r\left(r T_n\right) = \sum_n t_n \hat{w}_n T_n.
$$

As a shorthand notation, we could write:
$$
w = \texttt{R2DivR2} \cdot v.
$$

Clearly the presence of a generic function $f$ prevents a sparse approach. The above problem is approached pseudo-spectrally:

- First the product $f(r) \partial_r\left(r T_n\right)$ is performed in physical space.

- Then the result is transformed back to spectral space.


# Adding a new operator
The procedure to follow to define and test a new operator is:

1. Implement new operator (for example copying an existing operator and modifying it)
2. Add the operator to `Components/DenseSM/DenseSM/Chebyshev/LinearMap/CMakeLists.txt`

3. Test the operator

See below for some details.

# Example: R1FD2R1 operator

This operator is meant to perform:

$$
r f(r) \partial^2_r(r v)
$$

The ```R1FD2R1``` operator was created modifying the (slightly simpler) ```R2DivR2FD1R1```operator. See comments in `R1FD2R1.cpp` and `R1FD2R1.hpp`.


## Test the new R1FD2R1 operator

The tests are performed as follows:

1. Create a `build` folder in the `QuICC` root directory. It is best to call it exactly just `build`. Not `build_test` or anything like that. Unless one were to modify the Mathematica script to place the references in a specific build folder (see below)
2. Add the new operator to `Components/DenseSM/TestSuite/DenseSM/Tests/Chebyshev/LinearMap/CMakeLists.txt`.

    Remeber to add the ULP and MPULP for the new test. This is essentially the maximum allowed numerical difference between the reference and the test values.

    **If the test fails because of ULP issues, this value can be set to a higher value to make the test pass**

3. Generate reference with the mathematica notebook `Components/DenseSM/TestSuite/Mathematica/ChebyshevTripleHarmonicReference.nb`
4. Build the test
5. Run the test

See below for details on each step:
