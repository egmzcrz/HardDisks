push!(LOAD_PATH, "/Users/enrique/Documents/HardDisks/src/")

using Documenter, HardDisks

makedocs(sitename="Hard Disks Project Documentation")

deploydocs(repo = "github.com/egmzcrz/HardDisks.jl.git")
