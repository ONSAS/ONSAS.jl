#!/bin/sh

julia -e 'import Pkg; Pkg.activate(get(ENV, "ONSAS", @__DIR__)); using PkgBenchmark; results = benchmarkpkg("ONSAS", retune=true); export_markdown("benchmark/results.md", results)'
