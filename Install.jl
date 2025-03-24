import Pkg as p 

packages = ["CSV", "DataFrames", "Interpolations", "Plots", "Debugger", 
            "HTTP", "JSON", "SpecialFunctions", "LinearAlgebra"]

for pkg in packages
    p.add(pkg)
end