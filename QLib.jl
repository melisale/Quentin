module QLib
using QTypes

include("modelDefinition.jl")

export initialiseModel, setInletBC, setOutletBC, finiteDifferencesSolver

#=====================
      FUNCTIONS
=====================#

function initialiseModel()
  #1- initialise vessel
  # physical constant
  β = E.*h.*sqrt(π)./(0.75*π.*R0.*R0)

  # numerical constant
  Δx = ℓ/N_nodes

  # arrays to store iterative solution at current time step
  A = zeros(Float64, N_nodes) + π.*R0.*R0
  u = zeros(Float64, N_nodes)
  Q = zeros(Float64, N_nodes)
  P = zeros(Float64, N_nodes)
  c = zeros(Float64, N_nodes) + sqrt(sqrt(A).*β./ρ)

  v = QTypes.vesselData(vessel_name,
                    N_nodes, Δx,
                    β, R0, P_ext,
                    A, u, Q, P, c)

  #2- initialise blood
  b = QTypes.bloodData(μ, ρ)

  #3- initialise heart
  hrt = QTypes.heartData(cardiac_period, flow_amplitude,
                         inlet_function)

  #4- initialise capillaries
  cap = QTypes.capillariesData(Rc, Rd, Cd)

  #5- compute Δt from CFL condition; assume c >> u
  Δt = minimum(Δx./c)

  num_model = QTypes.numericalModel(Δt, iterations_limit,
                                    results_tollerance)

  return v, b, hrt, cap, num_model
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

function setOutletBC(delta_t, vessel_data, capillaries_data, blood_data)
  threeElementWindkessel(delta_t, vessel_data, capillaries_data, blood_data)
end

function threeElementWindkessel(delta_t, vessel_data, capillaries_data, blood_data)
  v   = vessel_data
  b   = blood_data
  wk  = capillaries_data
  Δt  = delta_t
  #----#

  v.Q[end] = v.Q[end-1]

  dtRC = Δt/(wk.Rd * wk.Cd)
  Pt = v.P[end] - wk.Rc * v.Q[end]
  Pt = (Pt + dtRC * v.Q[end] * wk.Rd) / (1 + dtRC)

  v.P[end] = Pt + v.Q[end]*wk.Rc

  At = (v.P[end] - v.P_ext)/v.beta[end]

  v.A[end] = At*At

  v.u[end] = v.Q[end] / v.A[end]
  v.c[end] = sqrt(sqrt(v.A[end])*v.beta[end]/b.rho)

# 	m.λ1[m.N] = (m.p[m.N] + m.ρ * m.c[m.N] * m.u[m.N])/2.
# 	m.λ2[m.N] = (m.p[m.N] - m.ρ * m.c[m.N] * m.u[m.N])/2.
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#========
 Solver
========#

function finiteDifferencesSolver(idx, vessel_data, blood_data, delta_t)
  MacCormackSolver(idx, vessel_data, blood_data, delta_t)
end

# without convective term
function MacCormackSolver(idx, vessel_data, blood_data, delta_t)
  v = vessel_data
  b = blood_data
  Δt = delta_t
  Δx = v.dx
  #-------------#

  # barred quantities
  Ab = zeros(Float64, v.N) # ̄A
  ub = zeros(Float64, v.N) # ̄u

  if idx == 1 # first iteration - inviscid flow and flat velocity profile

    # predictor step
    for i in 2:v.N-1

      dA_dt = -( v.u[i]*(v.A[i] - v.A[i-1]) + v.A[i]*(v.u[i] - v.u[i-1]) )/Δx
      du_dt = ( (v.u[i]*v.u[i]/v.A[i] -
                   v.beta[i]/(b.rho * sqrt(v.A[i])))*(v.A[i] - v.A[i-1]) +
                 v.u[i]*(v.u[i] * v.u[i-1]) )/Δx

      Ab[i] = v.A[i] + dA_dt*Δt
      ub[i] = v.u[i] + du_dt*Δt
    end

    Ab[end] = v.A[end]
    ub[end] = v.u[end]

    # corrector step
    for i in reverse([2:v.N-1])

      # barred derivatives
      dA_dt = -( ub[i]*(Ab[i+1] - Ab[i]) + Ab[i]*(ub[i+1] - ub[i]) )/Δx
      du_dt = ( (ub[i]*ub[i]/Ab[i] -
                   v.beta[i]/(b.rho * sqrt(Ab[i])))*(Ab[i+1] - Ab[i]) +
                 ub[i]*(ub[i+1] - ub[i]) )/Δx

      v.A[i] = 0.5*( v.A[i] + Ab[i] + dA_dt*Δt )
      v.u[i] = 0.5*( v.u[i] + ub[i] + dA_dt*Δt )
    end


  else # viscid flow and velocity profile from W-W theory
      return
  end


end



















#!!
end
