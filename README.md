# ONSAS.jl

---
## *An Open Nonlinear Structural Analysis Solver*

| **Documentation** |**Status** |**License** |  **Community** |
|:-----------------:|:----------------------------:|:-----------:| :-----------:|
| [![docs-dev][dev-img]][dev-url] | [![CI][ci-img]][ci-url] [![codecov][cov-img]][cov-url] |  [![license][lic-img]][lic-url] | [![community][com-img]][com-url] |


[dev-img]: https://img.shields.io/badge/docs-stable-blue.svg
[dev-url]: https://ONSAS.github.io/ONSAS.jl/dev/
[ci-img]: https://github.com/ONSAS/ONSAS.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/ONSAS/ONSAS.jl/actions/workflows/CI.yml?query=branch%3Amain
[cov-img]: https://codecov.io/gh/ONSAS/ONSAS.jl/branch/main/graph/badge.svg
[cov-url]: https://codecov.io/gh/ONSAS/ONSAS.jl
[doi-img]: https://zenodo.org/badge/170120732.svg
[doi-url]: https://zenodo.org/badge/latestdoi/
[lic-img]: https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000
[lic-url]: https://github.com/mvanzulli/ONSAS.jl/blob/main/LICENSE
[com-url]: https://app.slack.com/client/T04QWNG5T2Q/C04R6TMDV0R
[com-img]: https://user-images.githubusercontent.com/50339940/224228295-4ff7c201-975f-4e87-af51-86036a10f289.svg



---
## 💾 [Installing](https://julialang.github.io/Pkg.jl/v1/managing-packages/)

`ONSAS.jl` is a registered Julia package and [as such you can install it](https://julialang.github.io/Pkg.jl/v1/managing-packages/) by activating the `pkg` mode (type `]`, and to leave it, type `<backspace>`),
followed by

```julia
pkg> add ONSAS
```

---
## 💻 Contributing
All forms of contributing are kindly welcome (bug reports, docs, features, suggestions, questions and discussions). See CONTRIBUTING or  for more details.

---
## 📓 [Manual](http://onsas.github.io/ONSAS.jl/dev/)
To help you get started with using ONSAS.jl, we have provided resources in the [docs](http://onsas.github.io/ONSAS.jl/dev/). This includes internal documentation for interfaces, algorithms, and examples

---
## ⚙️ [Software architecture](http://juliareach.github.io/ONSAS.jl/dev/)
Project has been partitioned into three major groups for simultaneous development.The input group for defining the structure(mesh, materials, elements and boundary conditions), the fea package for the Finite Element Analysis and the output to post-process and visualize the simulation results. For a better understanding of the system architecture, please refer to the [UML diagram](https://github.com/ONSAS/ONSAS.jl/blob/main/UML.drawio) provided in this repository. 

