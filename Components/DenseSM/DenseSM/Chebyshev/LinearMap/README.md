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


### Mathematica reference

The easiest thing to do is to modify an existing reference. Some caveats/things to consider:

- Be sure that the reference is placed in the `build` directory. Some Mathematica ref files will not have a general path generation but it's hardcoded into the notebook and needs to be changed. 

    This is not currently (January 2025) a problem for the notebook being discussed now

- It is expensive to generate tests for all operators implemented in the notebook. Unless that is what you want, you might want to comment out operators that are not to be tested (comment out these lines in calls to `refTripleHarmonic`).



### Build and run the test

In the `build` directory (or wherever you put the tests):

-       cmake ..  -DQUICC_TESTSUITE_DENSESM=ON

-       make -j 6 DenseSMChebyshevTests

-       ctest -R _R1FD2R1 --output-on-failure


# Mathematica reference: calculation and output

For each of these operators $O_d$ the tests:

- take as input:
    - a radial function $f$
    - maximum degrees $nNr$ and $nNc$ (columns and rows).

- gives in output spectral coefficients, organised in matrix form, resulting from the projection of $O_d(T_j); \ j=0,\dots,nNc -1$ onto a chebyshev basis $T_i; \ i=0,\dots,nNr -1$.

    In other words, take an initial set of $T_j; \ j=0,\dots,nNc -1$, act on it with the dense operator $O_d$, project ont a basis $T_i; \ i=0,\dots,nNr -1$. 


## Calculation
The analytical forms of the operators are given in the functions, e.g., `r1Fd2r1Intg` and similar.

The integrals defined in the Table are of the form
$$
I_1 = \frac{2}{\pi}\int_0^1 \frac{1}{\sqrt{1-x^2}} T_i(x) \Phi dr
$$
where $\Phi(r)$ is the mathematical form of our operator. For `R1FD2R1`:
$$
\Phi = r f \partial_r^2(r T_j(x) \tau_j)
$$
and $\tau_j$ accomodates the Chebyshev polynomials normalisation (1 for $T_0$, 2 otherwise). The above integral needs to be converted in a way that only $x$ appears as the spatial variable, not $r$. To do so, use the following:

$$
r = ax+b.
$$
From which:
$$
\partial_x = a \ \partial_r;\quad dr = a \ dx
$$