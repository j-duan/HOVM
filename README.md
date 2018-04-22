# HOVM
Higher order variational models (HOVM) for image processing

## What is HOVM
HOVM includes up to 8 (+1 first order total variaiton model) higher order variational models for image denoising. It is written in matlab and super easy to run (download scripts and hit the run bottom directly in matlab). All of the concerned variational models are implemented by the fast alternating direction method of multipliers (i.e. split Bregman) with the finite difference discretistion. It avoids directly tackling the resulting higher order partial differential equations, which can be difficult to discretise to solve computationally. The main idea of the split Bregman is to break down the original problem into several subproblems, each of which can be solved analytically by using fast Fourier transform (FFT), soft-thresholding equations, etc. Therefore, the overall computational cost is low and the convergence rate is fast. The code has been made as straightforwardly as possible, so they shall be easy to understand by referring to the following corresponding literature.


## System requirements

The code uses only basic matlab built-in functions, which should be working across multiple versions of matlab (2013, 2015, 2017, etc).


## How to cite
If you find the code or a certain part of it useful, please consider giving appropriate credit to it by citing the following relevant papers. The code is easier to read by refering to [1] and [2]. The matlab scripts have same names as those in Table 1 in [1]. Thank you for your interest.

[1] Lu, W., Duan, J., Qiu, Z., Pan, Z., Liu, R. W., & Bai, L. (2016). Implementation of high‐order variational models made easy for image processing. Mathematical Methods in the Applied Sciences, 39 (14), 4208-4233. [[doi]](https://doi.org/10.1002/mma.3858)

[2] Duan, J., Qiu, Z., Lu, W., Wang, G., Pan, Z., & Bai, L. (2016). An edge-weighted second order variational model for image decomposition. Digital Signal Processing, 49, 162-181. [[doi]](https://doi.org/10.1016/j.dsp.2015.10.010)
