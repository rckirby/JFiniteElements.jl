using JFiniteElements
using Base.Test

# write your own tests here
@test 1 == 1

include("testTriangulation.jl")

include("testAssembly.jl")
