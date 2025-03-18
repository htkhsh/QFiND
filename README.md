# QFiND

This repository contains Python codes for constructing an effective discrete representation of a system-bath model. In other words, this codes provides an approximation of the bath correlation function $C(t)$ for a given spectral density $J(\omega)$

$$
\begin{aligned}
C(t)&=\frac{1}{2\pi} \int_{-\infty}^{\infty} \mathrm{d}\omega J(\omega)\left[\mathrm{coth}\left(\frac{\beta \omega}{2}\right)+1\right] \mathrm{e}^{-i \omega t}\\
&\approx \sum_{k=1}^M g_k^2 \mathrm{e}^{-i\omega_k t}
\end{aligned}
$$

where $`\omega_k,g_k \in ℝ \backslash \{0\}`$.
The code allows for the estimation of frequencies and coefficients in the system plus bosonic bath model using Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS). 



## Usage

1. Set the required parameters in an input file `input.txt` (see below).  
2. Run the main script:
   ```
   python ./src/qfind.py input.txt
   ```
3. The output will include the estimated frequencies and coefficients (saved as `omega_g.txt`), along with a plot of the resultant BCF (saved as `bcf.png`).


## Parameter Configuration

### Setting Parameters

To customize the simulation, you need to adjust certain parameters in the following files:

- **`input.txt`**: This file contains important global parameters such as:
  - `method`: Discretization method, such as:
    - `BSDO`: BSDO method
    - `ID`: ID approach
    - `LOG`: Logarithmic discretization
    - `MDM`: Mode density method
  - `temperature`: Specifies the temperature of the system in $[\mathrm{K}]$.
  - `Tc` (double): Cutoff time in $[\mathrm{fs}]$.
  - `Omega_min` (double): Minimum cutoff frequency in $[\mathrm{cm}^{-1}]$.
  - `Omega_max` (double): Maximum cutoff frequency in $[\mathrm{cm}^{-1}]$.
  - `N_t` (integer): Number of sample points in the time domain.
  - `N_w` (integer): Number of sample points in the frequency domain.
  - `wmax_quad`: The maximum frequency cutoff used in the numerical integration.
  - `eps` (double): Threshold for the ID.
  - `frank` (integer): Rank for the ID.  When frank is set to a value larger than 0 (`frank`>0), ID is performed based on the rank.
  - `stype`: The type of spectral density $J(\omega)$ (`PWR`, `TM`, `BO`).  The program supports several types of spectral density, such as:
    - Power-law with exponential cutoff (`PWR`) 
      $$J(\omega)=\pi\alpha\omega_c^{1-s}\omega^s\mathrm{e}^{-\omega/\omega_c}$$
    - Sum of Tannor-Meyer type spectral densities (`TM`)
      $$J(\omega)=\sum_{j=1}^n \frac{4\Gamma_j\lambda_j(\Omega_j^2+\Gamma_j^2)\omega}{\left[(\omega+\Omega_j)^2+\Gamma_j^2\right]\left[(\omega-\Omega_j)^2+\Gamma_j^2\right]}$$
    - Sum of Brownian spectral densities (`BO`)
      $$J(\omega)=\sum_{j=1}^n 2\lambda_j\frac{\Gamma_j \Omega_j^2\omega}{(\omega^2-\Omega_j^2)^2+\Gamma_j^2\Omega_j^2}$$
  - Parameters for specific spectral density types, such as:
    - `s`, `alpha`, `gamc` for Power-law Exponential (`PWR`).
    - `Omg`, `Gam`, `Lam` for Tannor-Meyer type (`TM`) and Brownian Oscillator (`BO`).

### Example Parameter Settings
See the directory `./examples`.

## Cite `QFiND`
If you find the framework useful in your research, we would be grateful if you could cite our publications:
- H. Takahashi and R. Borrelli, J. Chem. Phys. 161, 151101 (2024) (https://doi.org/10.1063/5.0232232) 
- H. Takahashi and R. Borrelli,  J. Chem. Theory Comput. 21, 2206–2218 (2025) (https://doi.org/10.1021/acs.jctc.4c01728)

Here are the bibtex entries:
```bib
@article{TakahashiBorrelli2024JCP,
  title = {Effective modeling of open quantum systems by low-rank discretization of structured environments},
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

@article{TakahashiBorrelli2025JCTC,
  title = {Discretization of {{Structured Bosonic Environments}} at {{Finite Temperature}} by {{Interpolative Decomposition}}: {{Theory}} and {{Application}}},
  shorttitle = {Discretization of {{Structured Bosonic Environments}} at {{Finite Temperature}} by {{Interpolative Decomposition}}},
  author = {Takahashi, Hideaki and Borrelli, Raffaele},
  year = {2025},
  month = mar,
  journal = {Journal of Chemical Theory and Computation},
  volume = {21},
  number = {5},
  pages = {2206--2218},
  publisher = {American Chemical Society},
  issn = {1549-9618},
  doi = {10.1021/acs.jctc.4c01728},
}

```


## Authors

Hideaki Takahashi (hideaki.takahashi@unito.it)


## License

This project is distributed under the [BSD 3-clause License](./LICENSE.md).


