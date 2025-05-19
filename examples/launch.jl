# This script is used to set up the environment for running examples in the package.
# It activates the package environment, loads necessary packages
using Pkg

# List of packages to ensure are present for running examples
devtools = ["GLMakie"]

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
using GLMakie

# Load your package (assumed to be one directory up)
Pkg.develop(path = joinpath(@__DIR__, ".."))  # ensures the dev version is loaded
using DetailedBalance

loaded = Base.loaded_modules_array()
loaded_names = [string(Base.PkgId(mod).name) for mod in loaded]

println("✅ Example environment loaded! Ready to go.")
println("🔍 Packages currently loaded:")

for pkg in devtools
    if pkg in loaded_names
        println("  ✔️  $pkg is loaded")
    else
        println("  ❌  $pkg is NOT loaded")
    end
end