"Return the as variables prints from the REPL."
function capture_print(output)
    # Extract the number of nodes and elements as integers.
    nnodes = parse(Int, match(r"(\d+) nodes (\d+) elements", output).captures[1])
    nelems = parse(Int, match(r"(\d+) nodes (\d+) elements", output).captures[2])
    return nnodes, nelems
end

"Delete files with `extension` in a folder with `example_folder` path."
function delete_files(example_folder::String, extension::String)
    return foreach(rm, filter(endswith(extension), readdir(example_folder; join=true)))
end

"Creates example folders paths from the example name. This assumes that the example folder 
is named as the example, and the example script is bench_<example_name>.jl."
function joinpath_example_folder(example_name::String)
    example_folder = joinpath(pkgdir(ONSAS), "benchmark", example_name)
    bench_path = joinpath(example_folder, "bench_" * example_name * ".jl")
    return example_folder, bench_path
end
