# Benchmark Report for *ONSAS*

## Job Properties
* Time of benchmark: 29 Oct 2023 - 9:20
* Package commit: dirty
* Julia commit: bed2cd
* Julia command flags: None
* Environment variables: None

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                                                                                     | time            | GC time   | memory          | allocations |
|------------------------------------------------------------------------------------------------------------------------|----------------:|----------:|----------------:|------------:|
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 0.5, nelems = 2321, nnodes = 816, npoints = 1000"]`   | 147.156 ms (5%) |           |  59.58 MiB (1%) |     1882863 |
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 0.5, nelems = 2321, nnodes = 816, npoints = 125"]`    |  20.103 ms (5%) |           |   8.22 MiB (1%) |      262401 |
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 1.0, nelems = 10804, nnodes = 3092, npoints = 1000"]` | 675.662 ms (5%) | 10.921 ms | 261.83 MiB (1%) |     8510117 |
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 1.0, nelems = 10804, nnodes = 3092, npoints = 125"]`  |  94.399 ms (5%) |           |  35.34 MiB (1%) |     1151102 |
| `["linear_cylinder_internal_pressure", "solve!, ms = 0.5, nelems = 2321, nnodes = 816"]`                               | 111.429 ms (5%) |  6.336 ms | 174.21 MiB (1%) |      711754 |
| `["linear_cylinder_internal_pressure", "solve!, ms = 1.0, nelems = 10804, nnodes = 3092"]`                             |    1.937 s (5%) | 57.378 ms | 809.43 MiB (1%) |     3300405 |
| `["linear_cylinder_internal_pressure", "structure, ms = 0.5, nelems = 2321, nnodes = 816"]`                            |  66.666 ms (5%) |           |  12.59 MiB (1%) |      167123 |
| `["linear_cylinder_internal_pressure", "structure, ms = 1.0, nelems = 10804, nnodes = 3092"]`                          | 263.812 ms (5%) |           |  49.96 MiB (1%) |      668519 |
| `["linear_extension", "solve, ms = 0.1, nelems = 9239, nnodes = 2137"]`                                                | 462.767 ms (5%) | 43.934 ms | 688.37 MiB (1%) |     2790054 |
| `["linear_extension", "solve, ms = 0.2, nelems = 1316, nnodes = 398"]`                                                 |  56.211 ms (5%) |  6.260 ms |  98.15 MiB (1%) |      398717 |
| `["linear_extension", "solve, ms = 0.3, nelems = 579, nnodes = 201"]`                                                  |  22.264 ms (5%) |           |  43.33 MiB (1%) |      175760 |
| `["linear_extension", "solve, ms = 0.4, nelems = 256, nnodes = 107"]`                                                  |   9.352 ms (5%) |           |  19.29 MiB (1%) |       78027 |
| `["linear_extension", "solve, ms = 0.5, nelems = 144, nnodes = 62"]`                                                   |   5.276 ms (5%) |           |  10.84 MiB (1%) |       43971 |
| `["uniaxial_compression", "solve, ms = 0.1, nelems = 9239, nnodes = 2137"]`                                            | 933.792 ms (5%) | 39.755 ms | 539.15 MiB (1%) |     1811290 |
| `["uniaxial_compression", "solve, ms = 0.2, nelems = 1316, nnodes = 398"]`                                             | 122.751 ms (5%) |           |  76.80 MiB (1%) |      259376 |
| `["uniaxial_compression", "solve, ms = 0.3, nelems = 579, nnodes = 201"]`                                              |  52.257 ms (5%) |           |  33.88 MiB (1%) |      114479 |
| `["uniaxial_compression", "solve, ms = 0.4, nelems = 256, nnodes = 107"]`                                              |  23.003 ms (5%) |           |  15.11 MiB (1%) |       50940 |
| `["uniaxial_compression", "solve, ms = 0.5, nelems = 144, nnodes = 62"]`                                               |  12.902 ms (5%) |           |   8.47 MiB (1%) |       28722 |
| `["uniaxial_extension", "solve, ms = 0.1, nelems = 9239, nnodes = 2137"]`                                              | 864.477 ms (5%) | 26.018 ms | 378.16 MiB (1%) |     1598793 |
| `["uniaxial_extension", "solve, ms = 0.2, nelems = 1316, nnodes = 398"]`                                               | 112.273 ms (5%) |           |  53.86 MiB (1%) |      229108 |
| `["uniaxial_extension", "solve, ms = 0.3, nelems = 579, nnodes = 201"]`                                                |  49.050 ms (5%) |           |  23.79 MiB (1%) |      101162 |
| `["uniaxial_extension", "solve, ms = 0.4, nelems = 256, nnodes = 107"]`                                                |  21.291 ms (5%) |           |  10.65 MiB (1%) |       45052 |
| `["uniaxial_extension", "solve, ms = 0.5, nelems = 144, nnodes = 62"]`                                                 |  12.143 ms (5%) |           |   5.96 MiB (1%) |       25410 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["linear_cylinder_internal_pressure"]`
- `["linear_extension"]`
- `["uniaxial_compression"]`
- `["uniaxial_extension"]`

## Julia versioninfo
```
Julia Version 1.9.3
Commit bed2cd540a1 (2023-08-24 14:43 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin22.4.0)
  uname: Darwin 22.6.0 Darwin Kernel Version 22.6.0: Wed Jul  5 22:22:05 PDT 2023; root:xnu-8796.141.3~6/RELEASE_ARM64_T6000 arm64 arm
  CPU: Apple M1 Pro: 
                 speed         user         nice          sys         idle          irq
       #1-10  2400 MHz      14592 s          0 s       5587 s     136212 s          0 s
  Memory: 16.0 GB (722.40625 MB free)
  Uptime: 1560.0 sec
  Load Avg:  3.6640625  3.4326171875  3.17822265625
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-14.0.6 (ORCJIT, apple-m1)
  Threads: 1 on 8 virtual cores
```