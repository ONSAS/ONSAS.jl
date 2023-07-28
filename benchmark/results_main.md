# Benchmark Report for *ONSAS*

## Job Properties
* Time of benchmark: 28 Jul 2023 - 15:14
* Package commit: fa3c09
* Julia commit: e4ee48
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
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 0.5, nelems = 2318, nnodes = 815, npoints = 1000"]`   | 197.971 ms (5%) |  3.904 ms |  59.75 MiB (1%) |     1888246 |
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 0.5, nelems = 2318, nnodes = 815, npoints = 125"]`    |  25.150 ms (5%) |           |   7.99 MiB (1%) |      254965 |
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 1.0, nelems = 10738, nnodes = 3083, npoints = 1000"]` | 868.462 ms (5%) | 20.091 ms | 260.36 MiB (1%) |     8461935 |
| `["linear_cylinder_internal_pressure", "point_eval_handler, ms = 1.0, nelems = 10738, nnodes = 3083, npoints = 125"]`  | 114.946 ms (5%) |           |  34.23 MiB (1%) |     1114856 |
| `["linear_cylinder_internal_pressure", "solve!, ms = 0.5, nelems = 2318, nnodes = 815"]`                               | 476.596 ms (5%) | 15.553 ms | 176.75 MiB (1%) |      685352 |
| `["linear_cylinder_internal_pressure", "solve!, ms = 1.0, nelems = 10738, nnodes = 3083"]`                             |    6.789 s (5%) | 70.361 ms | 817.30 MiB (1%) |     3162421 |
| `["linear_cylinder_internal_pressure", "structure, ms = 0.5, nelems = 2318, nnodes = 815"]`                            |  85.243 ms (5%) |           |  12.59 MiB (1%) |      167008 |
| `["linear_cylinder_internal_pressure", "structure, ms = 1.0, nelems = 10738, nnodes = 3083"]`                          | 336.675 ms (5%) |  8.143 ms |  49.80 MiB (1%) |      666040 |
| `["linear_extension", "solve, ms = 0.1, nelems = 9239, nnodes = 2137"]`                                                | 594.309 ms (5%) | 57.288 ms | 699.26 MiB (1%) |     2688424 |
| `["linear_extension", "solve, ms = 0.2, nelems = 1316, nnodes = 398"]`                                                 |  71.715 ms (5%) |  7.744 ms |  99.72 MiB (1%) |      384240 |
| `["linear_extension", "solve, ms = 0.3, nelems = 579, nnodes = 201"]`                                                  |  32.142 ms (5%) |  4.087 ms |  44.03 MiB (1%) |      169390 |
| `["linear_extension", "solve, ms = 0.4, nelems = 256, nnodes = 107"]`                                                  |  11.736 ms (5%) |           |  19.59 MiB (1%) |       75210 |
| `["linear_extension", "solve, ms = 0.5, nelems = 144, nnodes = 62"]`                                                   |   6.249 ms (5%) |           |  11.01 MiB (1%) |       42386 |
| `["uniaxial_compression", "solve, ms = 0.1, nelems = 9239, nnodes = 2137"]`                                            | 536.751 ms (5%) | 46.366 ms | 556.67 MiB (1%) |     1654227 |
| `["uniaxial_compression", "solve, ms = 0.2, nelems = 1316, nnodes = 398"]`                                             |  61.741 ms (5%) |  3.979 ms |  79.32 MiB (1%) |      237004 |
| `["uniaxial_compression", "solve, ms = 0.3, nelems = 579, nnodes = 201"]`                                              |  23.839 ms (5%) |           |  34.99 MiB (1%) |      104636 |
| `["uniaxial_compression", "solve, ms = 0.4, nelems = 256, nnodes = 107"]`                                              |   9.797 ms (5%) |           |  15.60 MiB (1%) |       46588 |
| `["uniaxial_compression", "solve, ms = 0.5, nelems = 144, nnodes = 62"]`                                               |   4.966 ms (5%) |           |   8.75 MiB (1%) |       26274 |
| `["uniaxial_extension", "solve, ms = 0.1, nelems = 9239, nnodes = 2137"]`                                              | 645.027 ms (5%) | 34.516 ms | 395.39 MiB (1%) |     1450969 |
| `["uniaxial_extension", "solve, ms = 0.2, nelems = 1316, nnodes = 398"]`                                               |  86.997 ms (5%) |  4.100 ms |  56.34 MiB (1%) |      208052 |
| `["uniaxial_extension", "solve, ms = 0.3, nelems = 579, nnodes = 201"]`                                                |  34.952 ms (5%) |  3.945 ms |  24.88 MiB (1%) |       91898 |
| `["uniaxial_extension", "solve, ms = 0.4, nelems = 256, nnodes = 107"]`                                                |  11.495 ms (5%) |           |  11.13 MiB (1%) |       40956 |
| `["uniaxial_extension", "solve, ms = 0.5, nelems = 144, nnodes = 62"]`                                                 |   6.740 ms (5%) |           |   6.23 MiB (1%) |       23106 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["linear_cylinder_internal_pressure"]`
- `["linear_extension"]`
- `["uniaxial_compression"]`
- `["uniaxial_extension"]`

## Julia versioninfo
```
Julia Version 1.9.2
Commit e4ee485e909 (2023-07-05 09:39 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin22.4.0)
  uname: Darwin 21.4.0 Darwin Kernel Version 21.4.0: Mon Feb 21 20:35:58 PST 2022; root:xnu-8020.101.4~2/RELEASE_ARM64_T6000 x86_64 i386
  CPU: Apple M1 Pro: 
              speed         user         nice          sys         idle          irq
       #1  2400 MHz     188397 s          0 s     183670 s    1177879 s          0 s
       #2  2400 MHz     184700 s          0 s     165940 s    1199286 s          0 s
       #3  2400 MHz     220330 s          0 s      65208 s    1264389 s          0 s
       #4  2400 MHz     153453 s          0 s      43030 s    1353443 s          0 s
       #5  2400 MHz     103146 s          0 s      25207 s    1421574 s          0 s
       #6  2400 MHz      78635 s          0 s      16516 s    1454776 s          0 s
       #7  2400 MHz      42743 s          0 s       8396 s    1498787 s          0 s
       #8  2400 MHz      25366 s          0 s       4160 s    1520401 s          0 s
  Memory: 32.0 GB (171.8359375 MB free)
  Uptime: 286638.0 sec
  Load Avg:  4.10400390625  3.66796875  3.25
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-14.0.6 (ORCJIT, westmere)
  Threads: 1 on 8 virtual cores
```