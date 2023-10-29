#!/bin/bash

julia_command="julia"

$julia_command -e 'import Pkg; Pkg.add("PkgBenchmark"); Pkg.add("BenchmarkTools")'
time $julia_command -e 'import Pkg; Pkg.activate(get(ENV, "ONSAS", @__DIR__)); using PkgBenchmark; results = benchmarkpkg("ONSAS", retune=true); export_markdown("results.md", results)'

find . -name "*.msh" -type f -delete
