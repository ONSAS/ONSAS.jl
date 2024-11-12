using Aqua
using ONSAS

@testset "Aqua.jl" begin
    Aqua.test_all(ONSAS;
                  ambiguities=false,
                  unbound_args=true,
                  undefined_exports=true,
                  project_extras=true,
                  stale_deps=true,
                  deps_compat=true,
                  piracies=true,
                  persistent_tasks=true)
end
