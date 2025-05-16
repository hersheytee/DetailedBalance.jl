using Pkg

# List of dev-only packages to ensure are present
devtools = ["Revise", "GLMakie"]

# Install any missing dev packages
for pkg in devtools
    if !haskey(Pkg.installed(), Symbol(pkg))
        println("📦 Installing missing package: $pkg")
        Pkg.add(pkg)
    end
end

# Activate the dev environment in the `dev/` folder
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load development tools
using Revise
using GLMakie

# Load your package (assumed to be one directory up)
Pkg.develop(path = joinpath(@__DIR__, ".."))  # ensures the dev version is loaded
using DetailedBalance

println("✅ Dev environment loaded! Revise + GLMakie + DetailedBalance are ready.")