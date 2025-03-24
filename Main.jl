using CSV, DataFrames, Interpolations, Plots
include("Constants.jl")
# Placeholder indexes
const iH2O = 0
const iCH4 = 2
const iSO2 = 3
const iCO2 = 1
# === Load Atmospheric Data ===
df_atm = CSV.read("Data/atmosphere_data.csv", DataFrame, types=Float64, stripwhitespace=true, delim=',')  # Load atmospheric data from a CSV file


# Extract necessary columns
z = df_atm.altitude   # Altitude in km
P = df_atm.pressure   # Pressure in Pa
T = df_atm.temperature  # Temperature in K
species = [:H2O, :CO2, :CH4, :SO2]  # List of species (update based on your dataset)
# n_species = Dict(sp => df_species[!, sp] for sp in species)  # Number densities per species (molecules/cm³)
df_species = CSV.read("Data/Species_densities_test.csv", DataFrame, types=Float64, stripwhitespace=true, delim=',')  # Load atmospheric data from a CSV file

# === Load Opacity Data (Cross-Sections in cm²/molecule) ===
function load_opacity(species::Symbol)
    # Placeholder: Replace with real cross-section data
    λ = range(0.3, 5.0, length=100)  # Wavelength range in microns
    σ = exp.(-λ) * 1e-22  # Fake cross-section 
    return λ, σ
end

# Interpolate opacities for each species
opacity_data = Dict(sp => linear_interpolation(load_opacity(sp)...) for sp in species) 

# === Compute Optical Depth ===
function optical_depth(λ::Float64, b::Float64)
    τ = 0.0  # Initialize optical depth
    for i in 1:length(z)-1
        dz = abs(z[i+1] - z[i])  # Layer thickness
        if b >= (z[i] + dz)
            return 0.0  # No contribution if the impact parameter is larger than the layer
        end
        ds = sqrt((z[i] + dz)^2 - b^2) - sqrt(z[i]^2 - b^2)  # Path length through layer
        if isnan(ds) || ds < 0
            continue
        end
        
        for sp in species
            σ = opacity_data[sp](λ)  # Interpolated cross-section
            for i in range(1,length(df_species[!,sp]))
                n = df_species[!,sp][i] * 1e6  # Convert cm³ to m³
                τ += σ * n * ds  # Summing contributions from all species
            end
        end
    end
    return τ
end

# === Compute Effective Transit Radius ===
function transit_radius(λ::Float64, R_p::Float64)
    b_vals = range(0, R_p + maximum(z), length=100)  # Impact parameters
    A = 0.0  # Initialize projected area
    for b in b_vals
        τ = optical_depth(λ, b)
        A += (1 - exp(-τ)) * 2π * b * (b_vals[2] - b_vals[1])  # Integral over impact parameter
    end
    return sqrt(A / π)  # Effective radius
end

# === Compute Transit Spectrum ===
function transit_spectrum(R_p::Float64, R_star::Float64)
    λ_range = range(0.3, 5.0, length=100)  # Wavelength range in microns
    R_eff = [transit_radius(λ, R_p) for λ in λ_range]  # Compute R_p(λ)
    print(R_eff)
    transit_depth = (R_eff ./ R_star) .^ 2  # Fractional transit depth
    print(transit_depth)

    # Plot results
    plot(λ_range, transit_depth, xlabel="Wavelength (µm)", ylabel="Transit Depth", lw=2, label="Simulated Spectrum", show = true)
    savefig("Plots/test.png")
    return λ_range, transit_depth
end

# === Run Simulation ===
R_p = 317*R_EARTH  # Planet radius in meters
R_star = R_SUN  # Stellar radius in meters
λ_vals, transit_depth_vals = transit_spectrum(R_p, R_star)
