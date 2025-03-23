import Pkg as p 

packages = ["CSV", "DataFrames", "Interpolations", "Plots", "Debugger"]

for pkg in packages
    p.add(pkg)
end