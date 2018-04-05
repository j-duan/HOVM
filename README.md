# HOVM
Higher order variational models (HOVM) for image processing

## What is HOVMID
HOVM includes up to 8 (+1 first order total variaiton model) higher order variational models for image denoising. It was written in matlab and super easy to run. These variational models are implemented by the fast alternating direction method of multipliers (i.e. split Bregman) with the finite difference discretistion. It avoids directly tackling the resulting higher order partial differential equations, which can be difficult to discretise to solve computationally. The main idea of the split Bregman is to break down the original problem into several subproblems, each of which can be solved analytically by using fast Fourier transform (FFT), soft-thresholding equations, etc. Therefore, the overall computational cost is low and the convergence rate is fast. The code has been made as straightforwardly as possible, so they shall be easy to understand by referring to the corresponding literature.


## System requirements

The code uses only basic matlab built-in functions, which should be working across multiple versions of matlab (2013, 2015, 2017, etc).


## How to cite
In the event you find the code or a certain part of it useful, please consider giving appropriate credit to it by citing one or some of the following relevant papers. The code is easier to read by refering to [1] and [2]. Thank you.

[1] Lu, W., Duan, J., Qiu, Z., Pan, Z., Liu, R. W., & Bai, L. (2016). Implementation of high‐order variational models made easy for image processing. Mathematical Methods in the Applied Sciences, 39(14), 4208-4233.

[2] Duan, J., Qiu, Z., Lu, W., Wang, G., Pan, Z., & Bai, L. (2016). An edge-weighted second order variational model for image decomposition. Digital Signal Processing, 49, 162-181.

