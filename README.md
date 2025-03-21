[![Matlab](https://img.shields.io/badge/MATLAB-R2023a-BLUE.svg)](https://nl.mathworks.com/products/new_products/release2023a.html)

# Insights into the explainability of Lasso-based DeePC for nonlinear systems

Matlab implementation of data-driven control strategies as described in the article
"Insights into the explainability of Lasso-based DeePC for nonlinear systems."

## Prerequisites

You will need:

- `Matlab`
- `Matlab statistics_toolbox` (see https://nl.mathworks.com/products/statistics.html)
- `CVX` (see https://cvxr.com/)
- `Gurobi` (see https://www.gurobi.com/)

## Installation

To clone this repository, see https://nl.mathworks.com/help/simulink/ug/clone-git-repository.html

## Benchmark case study: the unbalanced disk

The first-principle model of the DC motor connected to an unbalanced disk system is described as follows

```math
\begin{equation}
    \ddot{\theta}_t = \frac{mgl}{J}\text{cos}(\theta_t)-\frac{1}{\tau}\dot{\theta}_t + \frac{K_m}{\tau}v_t,
\end{equation}
```
being the angular position $\theta_t$ and the input voltage $v_t$ the measured output $y_t =\theta_t$ and applied input $u_t = v_t$, respectively. 

Toward an explainable Lasso-based DeePC, we construct a two-block grouped Hankel matrix by utilizing the (raw) input/output data gathered around the system's operating points (OPs) and stored in [`data`](data).

<p align="center">
  <img src="imgs/grouped_Hankel.png" width="60%" alt='A two-block grouped Hankel data structure'>
</p>

The following gifs compare the reference tracking performed by a Lasso-DeePC formulation with the above data structure.

<p align="center">
     <img src="gifs/IO_tr_C12.gif" alt="iLasso-DeePC trajectory tracking, data selection, and BPIs">
</p> 

## License
This project is licensed under the terms of the `CC-BY-4.0` license.
See [LICENSE](LICENSE) for more details.


## References
-
