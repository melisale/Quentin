using QLib

vessel, blood, heart, capillaries = initialiseModel()

setInletBC(0., heart, vessel)
