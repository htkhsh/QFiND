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


# Effective Discrete Renpresentation of Spectral Density Based on Interpolative Decomposition

This repository contains Python codes for constructing an effective system-bath model. The code allows for the estimation of frequencies and coefficients in the system plus bosonic bath model using Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS).

$$C(t)=\frac{1}{2\pi} \int_{-\infty}^{\infty} \mathrm{d}\omega J(\omega)\left[\mathrm{coth}\left(\frac{\beta \omega}{2}\right)+1\right] \mathrm{e}^{-i \omega t}$$

## Usage

1. Set the required parameters in `my_function.py` (see below).  
2. Run the main script:
   ```
   python ./src/main.py
   ```
3. The output will include the estimated frequencies and coefficients (saved as `omega_g.txt`), along with a plot of the resultant BCF (saved as `bcf.png`).


## Parameter Configuration

### Setting Parameters

To customize the simulation, you need to adjust certain parameters in the following files:

- **`global_value.py`**: This file contains important global parameters such as:
  - `temperature`: Specifies the temperature of the system in Kelvin.
  - `tc` (double): Cutoff time in [fs].
  - `omegac` (double): Cutoff frequency in [cm^{-1}].
  - `M` (integer): Number of sample points in the time domain.
  - `N` (integer): Number of sample points in the frequency domain.
  - `eps` (double): Threshold for numerical calculations.
  - `stype`: The type of spectral density model (`"PWR"`, `"TMn"`, `"BOn"`).
  - `wmax`: The maximum frequency cutoff for the spectral density calculation.
  - Parameters for specific spectral density types, such as:
    - `s`, `alpha`, `gamc` for Power-law Exponential (`PWR`).
    - `Omg`, `Gam`, `Lam` for Tannor-Meyer (`TMn`) and Brownian Oscillator (`BOn`).

### Example Parameter Settings

1. **Temperature Setting**: Set the system temperature in `global_value.py`:
   ```python
   temperature = 300.0  # Temperature in Kelvin
   ```

2. **Spectral Density Type**: Choose the spectral density type based on the physical model:
   ```python
   stype = "PWR"  # Power-law exponential decay model
   ```

3. **Spectral Density Parameters**: For the Power-law model, adjust parameters like:
   ```python
   s = 1.0
   alpha = 0.5
   gamc = 1.0
   ```


## Functionality

### Spectral Density Types

The program supports several types of spectral density, such as:
- Power-law with exponential cutoff (`PWR`)
- Tannor-Meyer type (`TMn`)
- Brownian oscillator (`BOn`)

These are specified in `global_value.py`.

### Frequency Estimation

The program uses Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS) for estimating frequencies from the time-domain correlation data. The method is implemented in `edr.py`.

### Plotting

The Bath Correlation Function (BCF) is visualized using `matplotlib`. The plotting function in `plot.py` generates a plot of the BCF and the approximation error and saves it as `bcf.png`.



## Cite `EDR-ID`
If you like `EDR-ID`, we would appreciate it if you starred the repository in order to help us increase its visibility. Furthermore, if you find the framework useful in your research, we would be grateful if you could cite our publications:
- H. Takahashi and R. Borrelli, J. Chem. Phys. 161, 151101 (2024). (https://doi.org/10.1063/5.0232232) 
- H. Takahashi and R. Borrelli, submitted to J. Chem. Phys.

Here are the bibtex entries:
```bib
@article{HierarchicalEOM-jl2023,
  doi = {10.1038/s42005-023-01427-2},
  url = {https://doi.org/10.1038/s42005-023-01427-2},
  year = {2023},
  month = {Oct},
  publisher = {Nature Portfolio},
  volume = {6},
  number = {1},
  pages = {313},
  author = {Huang, Yi-Te and Kuo, Po-Chen and Lambert, Neill and Cirio, Mauro and Cross, Simon and Yang, Shen-Liang and Nori, Franco and Chen, Yueh-Nan},
  title = {An efficient {J}ulia framework for hierarchical equations of motion in open quantum systems},
  journal = {Communications Physics}
}
```


## Authors

Hideaki Takahashi (takahashi.hideaki.w33@kyoto-u.jp)


## License

[![license](https://img.shields.io/badge/license-New%20BSD-blue.svg)](http://en.wikipedia.org/wiki/BSD_licenses#3-clause_license_.28.22Revised_BSD_License.22.2C_.22New_BSD_License.22.2C_or_.22Modified_BSD_License.22.29)

This project is distributed under the [BSD 3-clause License](./LICENSE.md).
