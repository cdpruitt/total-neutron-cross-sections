import sys

energy = float(sys.argv[1]) # in MeV

if(energy<=0):
    print("Error: energy = " + str(energy) + " MeV is not allowed")
    sys.exit()

distance = 2709.4 # in cm, for Ni experiment
neutronMass = 939.565 # in MeV/c^2
C = 2.9979*(10**8) # in m/s

velocity = ((1-(1/((energy/neutronMass)+1))**2)**0.5)*C # in m/s

if(velocity>C):
    print("Error: velocity cannot exceed speed of light")
    sys.exit()

TOF = ((10**7)*distance)/velocity

print "TOF = " + str(TOF) + " ns"
