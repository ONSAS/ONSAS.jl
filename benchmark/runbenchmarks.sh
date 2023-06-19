#!/bin/sh

time $JULIA_EXEC -e 'import Pkg; Pkg.activate(get(ENV, "ONSAS", @__DIR__)); using PkgBenchmark; results = benchmarkpkg("ONSAS", retune=true); export_markdown("results.md", results)'

find . -name "*.msh" -type f -delete
