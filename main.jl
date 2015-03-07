# load library
using QLib

# initialise model
vessel, blood, heart, capillaries, num_model = initialiseModel()

# Main loop

idx = 1
while true # solve until convergence

  for i in 0:(heart.T/num_model.delta_t - 1) # solve along one time step

    time = i*num_model.delta_t  # calculate current simulation time

    # Boundary Conditions
    setInletBC(time, heart, vessel)
    setOutletBC(num_model.delta_t, vessel, capillaries, blood)

    # Solver
    finiteDifferencesSolver(idx, vessel, blood, num_model.delta_t)
  end

  # Check convergence
  monitor = 1.
  if monitor <= num_model.toll
    println("Reached convergence!")
    break
  elseif idx >= num_model.it_lim
    println("$idx - Iterations limit reached without convergence D:")
    break
  else
    idx += 1
  end

end
