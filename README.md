<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [['$', '$'] ],
 displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
 }
 });
</script>


# edr-sd

This repository contains Python codes for constructing an effective discrete representation of a system-bath model. In other words, this codes provides an approximation of the bath correlation function $C(t) for a given spectral density $J(\omega)$
$$\begin{aligned}
C(t)&=\frac{1}{2\pi} \int_{-\infty}^{\infty} \mathrm{d}\omega J(\omega)\left[\mathrm{coth}\left(\frac{\beta \omega}{2}\right)+1\right] \mathrm{e}^{-i \omega t}\\
&\approx \sum_{k=1}^M g_k^2 \mathrm{e}^{-i\omega t},\;\omega_k,g_k\in\mathbb{R}\backslash\{0\}
\end{aligned}$$
The code allows for the estimation of frequencies and coefficients in the system plus bosonic bath model using Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS). 



## Usage

1. Set the required parameters in an input file `input.txt` (see below).  
2. Run the main script:
   ```
   python ./src/main.py input.txt
   ```
3. The output will include the estimated frequencies and coefficients (saved as `omega_g.txt`), along with a plot of the resultant BCF (saved as `bcf.png`).


## Parameter Configuration

### Setting Parameters

To customize the simulation, you need to adjust certain parameters in the following files:

- **`input.txt`**: This file contains important global parameters such as:
  - `method`: Discretization method, such as:
    - `ID`: ID approach
    - `BSDO`: BSDO method
  - `temperature`: Specifies the temperature of the system in [$\mathrm{K}$].
  - `Tc` (double): Cutoff time in [$\mathrm{fs}$].
  - `Omegac` (double): Cutoff frequency in [$\mathrm{cm}^{-1}$].
  - `N_t` (integer): Number of sample points in the time domain.
  - `N_w` (integer): Number of sample points in the frequency domain.
  - `wmax`: The maximum frequency cutoff used in the numerical integration of the FDT to calculate the BCF.
  - `eps` (double): Threshold for the ID.
  - `frank` (integer): Rank for the ID.  When frank is set to a value larger than -1 ($(\text{frank})>-1$), ID is performed based on the rank.
  - `stype`: The type of spectral density $J(\omega)$ (`"PWR"`, `"TMn"`, `"BOn"`).  The program supports several types of spectral density, such as:
    - Power-law with exponential cutoff (`PWR`) 
      $$J(\omega)=\pi\alpha\omega_c^{1-s}\omega^s\mathrm{e}^{-\omega/\omega_c}$$
    - Sum of Tannor-Meyer type spectral densities (`TMn`)
      $$J(\omega)=\sum_{j=1}^n \frac{4\Gamma_j\lambda_j(\Omega_j^2+\Gamma_j^2)\omega}{\left[(\omega+\Omega_j)^2+\Gamma_j^2\right]\left[(\omega-\Omega_j)^2+\Gamma_j^2\right]}$$
    - Sum of Brownian spectral densities (`BOn`)
      $$J(\omega)=\sum_{j=1}^n 2\lambda_j\frac{\zeta_j \Omega_j^2\omega}{(\omega^2-\Omega_j^2)^2+\zeta_j^2\Omega_j^2}$$
  - Parameters for specific spectral density types, such as:
    - `s`, `alpha`, `gamc` for Power-law Exponential (`PWR`).
    - `Omg`, `Gam`, `Lam` for Tannor-Meyer (`TMn`) and Brownian Oscillator (`BOn`).

### Example Parameter Settings


## Cite `EDR-ID`
If you find the framework useful in your research, we would be grateful if you could cite our publications:
- H. Takahashi and R. Borrelli, J. Chem. Phys. 161, 151101 (2024). (https://doi.org/10.1063/5.0232232) 
- H. Takahashi and R. Borrelli, submitted to J. Chem. Phys.

Here are the bibtex entries:
```bib
@article{TakahashiBorrelli2024JCP,
  title = {Effective Modeling of Open Quantum Systems by Low-Rank Discretization of Structured Environments},
  author = {Takahashi, Hideaki and Borrelli, Raffaele},
  year = {2024},
  month = oct,
  journal = {The Journal of Chemical Physics},
  volume = {161},
  number = {15},
  pages = {151101},
  issn = {0021-9606},
  doi = {10.1063/5.0232232}
}
```


## Authors

Hideaki Takahashi (takahashi.hideaki.w33@kyoto-u.jp)


## License

This project is distributed under the [BSD 3-clause License](./LICENSE.md).

## References


