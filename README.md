# Port-Hamiltonian Modeling of 3D Truss Structures and Frames
The code in this repository can be used to assemble 3D linear elastic mechanical port-Hamiltonian systems composed of beam and rod elements. 
A nonlinear model of a hydraulic piston actuator is also included that can be coupled with mechanical structures.
It was developed in the course of the [Collaborative Research Centre 1244](https://www.sfb1244.uni-stuttgart.de/en/) (CRC 1244) "*Adaptive skins and structures for the built environment of tomorrow*" funded by the German Resarch Foundation. Please also check the project website and learn about our reserach if you're interested. As of now, the framework in this repository is not actively futher developed by me, but you can still get in touch. Since the CRC 1244 is now in its second funding period, updates to this framework might follow sooner or later and will be linked to here (this might be delayed though). 

The objects [PH_LinearSystem](basic_objects/PH_LinearSystem.m) and [PH_NonlinearSystem](basic_objects/PH_NonlinearSystem.m) provide the core functionality of the framework.
Code for the spatial discretization of varios mechanical elements is contained in the folder [discretization](discretization).
A custom tool for calculating the mass and stiffness matrices of systems composed of truss elements and beams is included in [fem_tool](fem_tool).
It is used for comparing the spatial discretization approach used in the PH formulation to conventional FEM.

A detailed explanation of the methods and algorithms is available in the corresponding [paper](https://doi.org/10.1016/j.apm.2020.07.038) (accepted manuscript available on [arxiv](http://arxiv.org/abs/2008.07985)).
Also, please refer to the provided examples to learn about the usage of the framework. 

## Dependencies
The function [lgwt](https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes) written by Greg von Winckel is used to compute the Gauss-Legendre quadrature weights and nodes.
It is included in the repository so there is no need to download it.

A [simplectic integrator package](https://www.mathworks.com/matlabcentral/fileexchange/7686-symplectic-integrators) written by Francisco Beron-Vera is used in some of the examples.
We use the supplied symplectic s-stage Gauss-Legendre solver *gls* and a modified version called *linear_gls* that works faster for linear system.
Both are included in the repository so there is also no need to download the whole package.

For the simulation of nonlinear systems (fast calculation of Jacobians), [CasADi](https://web.casadi.org/) is used. 
Please download a recent version of the tool and add it to the matlab path to run the nonlinear examples. 

## Contributing
If you'd like to further develop the **ph_fem** code, add examples of its usage or otherwise contribute, feel free to [contact me](mailto:alexander.warsewa@posteo.net). I can add you as a contributer or get you in touch with other people working in the CRC 1244.

## License 
See the [LICENSE](LICENSE) file for license rights and limitations (GNU GPL3).
