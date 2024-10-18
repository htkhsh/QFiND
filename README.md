
# Bath Correlation Function (BCF) Simulation

This repository contains Python code for simulating the real and imaginary parts of the Bath Correlation Function (BCF) using various spectral density models. The code allows for the estimation of frequency spectra based on the time-correlation data using Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS).

## Installation

1. Clone the repository.
2. Ensure you have the required Python packages installed:
   ```bash
   pip install numpy scipy matplotlib
   ```

## Code Structure

### Main Files

- **`corrfunc.py`**: 
  - Defines functions to compute the real (`S_exact`) and imaginary (`A_exact`) parts of the BCF. It integrates the spectral density over frequency to evaluate these components.

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

1. Set the required parameters in `global_value.py` and `my_function.py`.
2. Run the main script:
   ```bash
   python main.py
   ```
3. The output will include the estimated frequencies and amplitudes, along with a plot of the Bath Correlation Function (saved as `bcf.png`).

## Functionality

### Spectral Density Types

The program supports several types of spectral density, such as:
- Power law exponential decay (`PWR`)
- Tannor-Meyer oscillator model (`TMn`)
- Brownian oscillator model (`BOn`)

These are specified in `global_value.py` and are used to compute the BCF in `corrfunc.py`.

### Frequency Estimation

The program uses Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS) for estimating frequencies from the time-domain correlation data. The method is implemented in `edr.py`.

### Plotting

The Bath Correlation Function (BCF) is visualized using `matplotlib`. The plotting function in `plot.py` generates a plot of the real and imaginary parts of the BCF and saves it as `bcf.png`.

## Example

You can modify the parameters in `global_value.py` to test different spectral density types and temperature settings. After running `main.py`, the results will be printed to the console, and a plot of the correlation function will be saved.

## License

This project is licensed under the MIT License.
