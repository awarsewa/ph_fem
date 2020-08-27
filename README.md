# Port-Hamiltonian Modeling of 3D Truss Structures and Frames
The code in this repository can be used to assemble 3D linear elastic mechanical port-Hamiltonian systems composed of beam and rod elements. 
A nonlinear model of a hydraulic piston actuator is also included that can be coupled with mechanical structures.

A detailed explanation of the methods and algorithms is available in the corresponding [paper](https://doi.org/10.1016/j.apm.2020.07.038.) (accepted manuscript available on [arxiv](http://arxiv.org/abs/2008.07985)).
Also, please refer to the provided examples to learn about the usage of the framework. 

## Dependencies
The function [lgwt](https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes) written by Greg von Winckel is used to compute the Gauss-Legendre quadrature weights and nodes.
It is included in the repository so there is no need to download it.

A [simplectic integrator package](https://www.mathworks.com/matlabcentral/fileexchange/7686-symplectic-integrators) written by Francisco Beron-Vera is used in some of the example.
We use the supplied symplectic s-stage Gauss-Legendre solver *gls* and a modified version called *linear_gls* that works faster for linear system.
Both are included in the repository so there is also no need to download the whole package.

For the simulation of nonlinear systems (fast calculation of Jacobians), [CasADi](https://web.casadi.org/) is used. 
Please download a recent version of the tool and add it to the matlab path to run the nonlinear examples. 

## License 
See the [LICENSE](LICENSE) file for license rights and limitations (GNU GPL3).