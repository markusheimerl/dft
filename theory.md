# **Neural Network Variational Monte Carlo: Rigorous Mathematical Formulation**  

---

## **1. Quantum Many-Body Problem**

### 1.1 Configuration Space  
• Electronic coordinates:                           $\mathbf r=(\mathbf r_1,\dots ,\mathbf r_N)\in\mathbb R^{3N}$  
• Clamped nuclear coordinates (Born–Oppenheimer): $\mathbf R=(\mathbf R_1,\dots ,\mathbf R_M)\in\mathbb R^{3M}$ (constant)  
• Nuclear charges:                                  $\mathbf Z=(Z_1,\dots ,Z_M)\in\mathbb Z_+^M$

### 1.2 Electronic Hamiltonian (in a.u.)  
$\displaystyle  
\hat H=\hat T_e+\hat V_{ee}+\hat V_{eN}+V_{NN}, \qquad  
\hat T_e=-\tfrac12\sum_{i=1}^N\nabla_{\mathbf r_i}^2
$  

Potential terms  
$\displaystyle  
\hat V_{ee}= \sum_{1\le i<j\le N}\frac1{|\mathbf r_i-\mathbf r_j|},\qquad  
\hat V_{eN}= -\sum_{i=1}^N\sum_{A=1}^M\frac{Z_A}{|\mathbf r_i-\mathbf R_A|},\qquad  
V_{NN}= \sum_{1\le A<B\le M}\frac{Z_A Z_B}{|\mathbf R_A-\mathbf R_B|}\;(\text{constant}).  
$

### 1.3 Time-Independent Schrödinger Equation  
$\hat H\Psi(\mathbf r,\boldsymbol\sigma)=E\Psi(\mathbf r,\boldsymbol\sigma)$  
with the fermionic antisymmetry condition under exchange of *both* space and spin coordinates  
$\Psi(\dots,\mathbf r_i\sigma_i,\dots,\mathbf r_j\sigma_j,\dots)= -\Psi(\dots,\mathbf r_j\sigma_j,\dots,\mathbf r_i\sigma_i,\dots)$.

---

## **2. Neural-Network Wave-Function Ansatz**

### 2.1 Parameterised Form  
$\displaystyle  
\Psi_\vartheta(\mathbf r,\boldsymbol\sigma)=\mathcal A\!\bigl[\Phi_\theta(\mathbf r,\boldsymbol\sigma)\bigr]\;
\mathrm e^{\,J_\phi(\mathbf r)},\qquad
\vartheta:=(\theta,\phi)\in\Theta\subset\mathbb R^p.
$  
• $\mathcal A$ builds an antisymmetric Slater-type part.  
• $J_\phi$ is a symmetric Jastrow factor implemented by a neural network.

### 2.2 Slater Determinant Construction  
Spin partition: $N_\uparrow=\lceil N/2\rceil,\;N_\downarrow=\lfloor N/2\rfloor$.  

A neural network $f_\theta:\mathbb R^{3}\!\times\mathbb R^{3N}\!\times\mathbb R^{3M}\!\times\mathbb Z_+^M\to\mathbb R^{K}$ outputs *K* (usually $K\!\ge\!N_\uparrow$) orbital values per electron:

$\displaystyle  
\mathbf o_i=f_\theta(\mathbf r_i;\mathbf r,\mathbf R,\mathbf Z)\in\mathbb R^{K}.
$  

Select the first $N_\uparrow$ (respectively $N_\downarrow$) orbitals to build the square sub-matrices  

$\displaystyle  
\mathbf M_\uparrow=\bigl[o_1^{(1:\!N_\uparrow)};\dots ;o_{N_\uparrow}^{(1:\!N_\uparrow)}\bigr],\qquad  
\mathbf M_\downarrow=\bigl[o_{N_\uparrow+1}^{(1:\!N_\downarrow)};\dots ;o_{N}^{(1:\!N_\downarrow)}\bigr].
$  

Slater part (antisymmetric in same-spin coordinates)  
$\displaystyle  
\mathcal A\!\bigl[\Phi_\theta(\mathbf r,\boldsymbol\sigma)\bigr]
       =\bigl[\det\mathbf M_\uparrow\bigr]\bigl[\det\mathbf M_\downarrow\bigr]\;
         \chi(\boldsymbol\sigma),
$  
where the spin function $\chi$ is the usual product of $N_\uparrow$ spin-up and $N_\downarrow$ spin-down states, guaranteeing full fermionic antisymmetry.

### 2.3 Jastrow Factor  
$\displaystyle  
J_\phi(\mathbf r)=g_\phi\!\left(\{|\mathbf r_i-\mathbf r_j|\}_{i<j},
                               \{|\mathbf r_i-\mathbf R_A|\}_{i,A}\right),
\qquad  
\mathcal J_\phi(\mathbf r)=\mathrm e^{J_\phi(\mathbf r)}.
$  
Short-range two-body terms needed to satisfy electron–electron and electron–nucleus cusp conditions can be hard-wired or learned.

---

## **3. Variational Principle**

For every *normalised* trial state ($\langle\Psi_\vartheta|\Psi_\vartheta\rangle=1$)  
$\displaystyle  
E[\vartheta]=\langle\Psi_\vartheta|\hat H|\Psi_\vartheta\rangle\;\ge\;E_0.
$  

Expectation values may be written with respect to the probability density  
$\displaystyle  
\rho_\vartheta(\mathbf r)=|\Psi_\vartheta(\mathbf r,\boldsymbol\sigma)|^2
$  
(the spin variables are summed over).

Local energy  
$\displaystyle  
E_{\mathrm{loc}}(\mathbf r;\vartheta)=\frac{\hat H\Psi_\vartheta(\mathbf r,\boldsymbol\sigma)}
                                           {\Psi_\vartheta(\mathbf r,\boldsymbol\sigma)}
                                          \in\mathbb R.  
$

---

## **4. Monte-Carlo Sampling**

Generate a Markov chain $\{\mathbf r^{(k)}\}_{k=1}^K$ with stationary density $\rho_\vartheta$.  
Metropolis–Hastings or force-biased Langevin (drift-diffusion) proposals are common:

1. Drifted proposal  
   $\mathbf r'\!=\!\mathbf r^{(k)}+\eta\,\mathbf F(\mathbf r^{(k)})+\sqrt{2\eta}\,\boldsymbol\xi$  
   with force $\mathbf F=\nabla_{\mathbf r}\log|\Psi_\vartheta|$ and $\boldsymbol\xi\sim\mathcal N(0,I)$.  
2. Accept with the usual Metropolis ratio to ensure detailed balance.

Energy estimator  
$\displaystyle  
\hat E_K[\vartheta]=\frac1K\sum_{k=1}^K E_{\mathrm{loc}}(\mathbf r^{(k)};\vartheta).
$  

Ergodic CLT for correlated samples  
$\displaystyle  
\sqrt{K}\bigl(\hat E_K[\vartheta]-E[\vartheta]\bigr)\;\xrightarrow{d}\;
\mathcal N\!\bigl(0,\,\sigma_{\mathrm{eff}}^2[\vartheta]\bigr),\qquad
\sigma_{\mathrm{eff}}^2=2\tau_{\mathrm{int}}\operatorname{Var}_{\rho_\vartheta}\!(E_{\mathrm{loc}}),
$  
where $\tau_{\mathrm{int}}$ is the integrated autocorrelation time.

---

## **5. Energy-Gradient Formula**

Because the Hamiltonian does not depend on $\vartheta$, the minimal-variance “score-function” identity holds:

$\displaystyle  
\nabla_\vartheta E
      =2\; \mathbb E_{\rho_\vartheta}\!\Bigl[(E_{\mathrm{loc}}-E)\;
                                                \nabla_\vartheta\log|\Psi_\vartheta|\Bigr].
$  

Monte-Carlo estimator (using the same configurations)  

$\displaystyle  
\widehat{\nabla_\vartheta E}
      =\frac{2}{K}\sum_{k=1}^K
         \Bigl(E_{\mathrm{loc}}^{(k)}-\hat E_K\Bigr)\;
         \nabla_\vartheta\log|\Psi_\vartheta(\mathbf r^{(k)})|.
$  

No explicit $\nabla_\vartheta E_{\mathrm{loc}}$ term is required.

---

## **6. Stochastic Optimisation**

Stochastic gradient descent (SGD) or adaptive variants (Adam, RMSProp, …):

$\displaystyle  
\vartheta^{(t+1)}=\vartheta^{(t)}-\alpha_t\;
                  \widehat{\nabla_\vartheta E}[\vartheta^{(t)}],
\qquad
\sum_t\alpha_t=\infty,\;
\sum_t\alpha_t^2<\infty.
$

Under these step-size conditions and standard regularity/mixing assumptions, every accumulation point of the sequence $\{\vartheta^{(t)}\}$ is almost surely a stationary point of $E$ (global optimality cannot be guaranteed for the non-convex landscape).

---

## **7. Complete Algorithm**

Input : nuclear data $(\mathbf R,\mathbf Z)$; network definitions $f_\theta$, $g_\phi$  
Output: parameters $\vartheta^\star$, ground-state energy estimate $E^\star$

```
1. Initialisation:
     draw ϑ⁽⁰⁾ = (θ⁽⁰⁾, φ⁽⁰⁾) ~ N(0,σ²I)
     sample electronic configuration r⁽⁰⁾ from some broad distribution
2. for t = 0 … T−1
     a) Generate K correlated samples {r⁽ᵏ⁾} via MCMC targeting ρ_{ϑ⁽ᵗ⁾}
     b) For each sample compute local energy E_loc⁽ᵏ⁾
     c) Energy estimate        Ē = (1/K) Σ_k E_loc⁽ᵏ⁾
     d) Gradient estimate      ∇̂ = (2/K) Σ_k (E_loc⁽ᵏ⁾ − Ē) ∇_ϑ log|Ψ_{ϑ⁽ᵗ⁾}(r⁽ᵏ⁾)|
     e) Parameter update       ϑ⁽ᵗ+¹⁾ = ϑ⁽ᵗ⁾ − α_t ∇̂
3. return ϑ* = ϑ⁽ᵀ⁾,  E* = Ē
```

Optional improvements  
• Enforce cusp conditions in $J_\phi$.  
• Use correlated-sampling reweighting when changing parameters.  
• Employ variance reduction (blocking, batching, control-variates).

---

## **8. Summary and Scope**

The Neural-Network VMC formalism above is strictly *ab initio* within  

1. the clamped-nuclei (Born–Oppenheimer) approximation, and  
2. the expressive power of the chosen parametrisation $\Psi_\vartheta$.  

Apart from these controllable approximations, no empirical data or model bias is introduced; all averages are evaluated stochastically but *exactly* with respect to the underlying quantum Hamiltonian.