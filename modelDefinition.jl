# vessel name
const vessel_name = "straight_vessel"
const ℓ = 0.02    # vessel length

# numerical constants
const N_nodes = 11
const iterations_limit = 1  # iterations limit
const results_tollerance = 1e-5

# physical constants
const E = 5.2e5  # Young's modulus
const h  = zeros(Float64, N_nodes) + 1.e-3  # h(x) thickness
const R0 = zeros(Float64, N_nodes) + 1.e-2  # R0(x) diastolic lumen radius
const P_ext = 2*133.32

# boundary conditions
#=
1: half sin function
2: from .dat file
=#
const inlet_function = 1
const cardiac_period = 1.   # [s]

#= !!!
Always put a value for flow_amplitude. In the case of
inlet_function == 2, put 0.
Otherwise, initialiseModel() returns error because it cannot
initilise QTypes.heartData()
=#
const flow_amplitude = 3.e-4 # [l/s]

# outlet 3wk parameters
const Rc = 1.e7
const Rd = 1.e8
const Cd = 1.e-8

# blood parameters
const ρ = 1056.
const μ = 2.5e-2
