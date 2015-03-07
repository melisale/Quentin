module QLib
using QTypes

include("modelDefinition.jl")

export initialiseModel, setInletBC

#=====================
 FUNCTIONS DEFINITION
=====================#

function initialiseModel()
  #1- initialise vessel
  # physical constant
  β = E.*h.*sqrt(π)./(0.75*π.*R0.*R0)

  # numerical constant
  Δx = ℓ/N_nodes

  # arrays to store iterative solution at current time step
  A = zeros(Float64, N_nodes)
  u = zeros(Float64, N_nodes)
  Q = zeros(Float64, N_nodes)
  P = zeros(Float64, N_nodes)
  c = zeros(Float64, N_nodes)

  v = QTypes.vesselData(vessel_name,
                    N_nodes, Δx,
                    β, R0,
                    A, u, Q, P, c)

  #2- initialise blood
  b = QTypes.bloodData(μ, ρ)

  #3- initialise heart
  hrt = QTypes.heartData(cardiac_period, flow_amplitude,
                         inlet_function)

  #4- initialise capillaries
  cap = QTypes.capillariesData(Rc, Rd, Cd)

  return v, b, hrt, cap
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#==========
 Inlet BC
==========#

# Set inlet flow at each time accordingly to inlet_function
function setInletBC(time, heart_data, vessel_data)
  hrt = heart_data
  v   = vessel_data

  if hrt.i_func == 1 # half sin
    v.Q[1] = getHalfSin(time, hrt.T, hrt.Aq)

  elseif hrt.i_func == 2 # .dat
      println("ERROR: inlet data not defined!")
  end
end

#= Inlet functions =#
function getHalfSin(time, period, q_amplitude)
  if time <= period*0.3
    return q_amplitude * sin( time*π / (0.3*period) )
  else
    return 0.
  end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===========
 Outlet BC
===========#

function setOutletBC(time, vessel_data, capillaries_data)
  v   = vessel_data
  cap = capillaries_data



end

function threeElementWindkessel(mt, vessel, delta_t, Rc, Rd, Cd)
  m = vessel
  DT= delta_t

#   Rc = vessel.c[end] * vessel.ρ / vessel.a[end]

	m.q[m.N] = m.q[m.N-1]

  tv = DT/(Rd*Cd)

	pt = m.p[m.N] - Rc*m.q[m.N]
	pt = (pt + tv*m.q[m.N]*Rd) / (1+tv)

	m.p[m.N] = pt + m.q[m.N]*Rc

	ta = m.sa0 + (m.p[m.N] - P_ext)/m.G0

  m.a[m.N] = ta*ta

	m.u[m.N] = m.q[m.N] / m.a[m.N]
	m.c[m.N] = sqrt(m.G0 * sqrt(m.a[m.N]) / (2.*m.ρ))

	m.λ1[m.N] = (m.p[m.N] + m.ρ * m.c[m.N] * m.u[m.N])/2.
	m.λ2[m.N] = (m.p[m.N] - m.ρ * m.c[m.N] * m.u[m.N])/2.
end















end
