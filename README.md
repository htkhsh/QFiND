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

## Code Structure

### Main Files

- **`corrfunc.py`**: 
  - Defines functions to compute the real (`S_exact`) and imaginary (`A_exact`) parts of the BCF. It numerically integrates the fluctuation-dissipation theorem (FDT) to evaluate these components.

- **`edr.py`**: 
  - Implements frequency estimation using Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS). This allows for decomposing time-domain data to estimate its frequency content.

- **`eval.py`**: 
  - Contains functions to evaluate the frequency spectrum based on time correlation functions, using the methods implemented in `edr.py`.

- **`global_value.py`**: 
  - Contains global configuration parameters (e.g., temperature, spectral density type, and max frequency). These are used across the other modules to ensure consistent simulation settings.

- **`main.py`**: 
  - The main script to run the simulation. It initializes the required parameters, runs frequency estimations, and generates plots of the correlation functions.

- **`my_function.py`**: 
  - Utility functions for setting up parameters for the simulation.

- **`plot.py`**: 
  - Defines the plotting logic to visualize the real and imaginary parts of the BCF.

- **`specdens.py`**: 
  - Contains functions for computing different types of spectral density (power law, Tannor-Meyer, Brownian). These are used in the correlation function calculations.

## Usage

1. Set the required parameters in `my_function.py` (see below).  
2. Run the main script:
   ```
   python main.py
   ```
3. The output will include the estimated frequencies and coefficients (saved as `omega_g.txt`), along with a plot of the resultant bath correlation function (saved as `bcf.png`).

## Parameter Configuration

### Setting Parameters

To customize the simulation, you need to adjust certain parameters in the following files:

- **`global_value.py`**: This file contains important global parameters such as:
  - `temperature`: Specifies the temperature of the system in Kelvin.
  - `stype`: The type of spectral density model (`"PWR"`, `"TMn"`, `"BOn"`).
  - `wmax`: The maximum frequency cutoff for the spectral density calculation.
  - Parameters for specific spectral density types, such as:
    - `s`, `alpha`, `gamc` for Power-law Exponential (`PWR`).
    - `Omg`, `Gam`, `Lam` for Tannor-Meyer (`TMn`) and Brownian Oscillator (`BOn`).

- **`my_function.py`**: Use the `setpara()` function to set up parameters needed by the simulation. Customize the default values within this function for your use case.

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

## Authors
Hideaki Takahashi (takahashi.hideaki.w33@kyoto-u.jp)

## Cite `EDR-ID`
If you like `EDR-ID`, we would appreciate it if you starred the repository in order to help us increase its visibility. Furthermore, if you find the framework useful in your research, we would be grateful if you could cite our publications
 [ [Commun. Phys. 6, 313 (2023)](https://doi.org/10.1038/s42005-023-01427-2)  ] 
using the following bibtex entry:
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

## License

This project is licensed under the MIT License.
