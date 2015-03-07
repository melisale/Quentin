module QTypes

export vesselData, bloodData, heartData, capillariesData

#=================
 TYPES DEFINITION
==================#

type vesselData
  label :: String   # vessel reference name

  # numerical constants
  N  :: Int64     # number of nodes
  dx :: Float64   # Δx

  # physical constants
  beta :: Array{Float64, 1}   # β
  R0   :: Array{Float64, 1}   # unstressed lumen radius

  # iterative solution
  A :: Array{Float64, 1}   # cross-sectional area
  u :: Array{Float64, 1}   # axial velocity
  Q :: Array{Float64, 1}   # flow
  P :: Array{Float64, 1}   # pressure
  c :: Array{Float64, 1}   # wave velocity
end

type bloodData
  mu  :: Float64   # μ dynamic viscosity
  rho :: Float64   # ρ density
end

type heartData
  T  :: Float64  # cardiac period
  Aq :: Float64  # flow wave amplitude

  #=
  1: half sin function of half-period 0.3T
  2: from .dat
  =#
  i_func :: Int16 # inlet punping function
end

type capillariesData
  Rc :: Float64
  Rd :: Float64
  Cd :: Float64
end

end
