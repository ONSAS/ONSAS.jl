mutable struct MaterialsData{T}
    youngModulus::T
    poissonRatio::T
end

struct ElementsData
    type::String
    crossSectionArea::Float64
end

mutable struct BoundaryCondsData
    DirichletNodalDOFs::Vector{Int}
    DirichletNodalVals::Vector{Float64}
    NeumannNodalDOFs::Vector{Int}
    NeumannNodalVals::Vector{Float64}
end

struct BCsIndexes
    diriDofs::Vector{Int}
    neumDofs::Vector{Int}
end

struct SystemMatrix{T}
    matrix::T
end
