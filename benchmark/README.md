## How to run the benchmarks

The files in this folder define a benchmark suite with the tools provided by
[PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl) and
[BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl).

To run the benchmarks, execute:

```julia
julia> using PkgBenchmark

julia> results = benchmarkpkg("ONSAS", retune=true)
```

A shell script, `runbenchmarks.sh`, is also provided. In that case, ensure to run
the script from the directory containing your ONSAS development environment, or otherwise define the environment variable `ENV["ONSAS"]` containing such environment. Moreover, `julia_command` needs to point to your Julia binary.

## How to compare benchmarks

To compare current version to another tagged version, commit or branch:

```julia
julia> results = judge("ONSAS", <tagged-version-or-branch>)
```

## Exporting results

To export the benchmark results to a JSON file:

```julia
julia> writeresults("results.json", results)
```

To export the benchmark results to a Markdown file:

```julia
julia> export_markdown("results.md", results)
```
