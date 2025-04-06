**This code provides all numerical calculations and generates Figures 2–6 in the paper.**  
The code is organized into four parts as follows:

---

**PART I – Function Definitions and Core Computations**  
- `R11(gamma)`: Computes the approximation ratio of `APPROX.1(λ, 0.5, p)`  
- `R12(gamma)`: Computes the approximation ratio of `APPROX.1(λ, 0.6677, p)`  
- `R13(gamma)`: Computes the approximation ratio of `APPROX.1(λ, θ, p)`  
- `R21(gamma)`: Computes the approximation ratio of `APPROX.2`  
- `R22(gamma)`: Computes the maximum value obtained from the linear programs in `APPROX.2(σ = 1)`  
- `lp1(gamma, sigma)`: Computes the value of the first linear program  
- `lp2(gamma, sigma)`: Computes the value of the second linear program  

---

**PART II – Generating Figures 2, 3, and 6**  
This section uses the above functions to produce the data and plots for Figures 2, 3, and 6.

---

**PART III – Generating Figure 4**  
This section performs additional computations and plotting required specifically for Figure 4.

---

**PART IV – Generating Figure 5**  
This section is dedicated to producing Figure 5 based on numerical evaluations.
