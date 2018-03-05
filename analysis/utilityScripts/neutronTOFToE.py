import sys

TOF = float(sys.argv[1]) # in ns

if(TOF<=0):
    print("Error: TOF = " + str(TOF) + " ns is not allowed")
    sys.exit()

distance = 2709.4 # in cm, for Ni experiment
neutronMass = 939.565 # in MeV/c^2
C = 2.9979*(10**8) # in m/s

velocity = ((10**7)*distance)/TOF # in m/s

if(velocity>C):
    print("Error: velocity cannot exceed speed of light")
    sys.exit()

energy = (((1-((velocity/C)**2))**-0.5)-1)*neutronMass
print "Energy = " + str(energy) + " MeV"
