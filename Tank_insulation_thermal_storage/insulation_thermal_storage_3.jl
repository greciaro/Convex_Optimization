
using JuMP
using AmplNLWriter


###### PARAMETERS and DATA ###################################
#Salt
T0Initial = 520 + 273.15 #[K]
T0Min = 131 + 273.15 #[K]
#T0 at the frist hour with exponential decay
T0 = 789.8331 #[K]
DensitySalt = 1870 #[kg/m^3]
HeatCapacitySalt = 1.6 #[kJ/(kg-K)]
CostSalt = 0.5 #[$/kg]
#Steel
DensitySteel = 7.8 #[tonnes/m^3]
ThermalConductivitySTeel = 40 #[W/(m-K)]
CostSteel = 4000 #[$/tonne]
#Insultion
ThermalConductivityInsulation = 0.4 #[W/(m-K)]
ConvectiveHeatTransferInsulation = 15 #[W/(m^2-K)]
CostInsulation = 800 #[$/m^3]
#Physical Constants
StefanBoltzmanCosntant = 5.67*10^-8 #[W/(m^2-K^4)]
#Tank
HeightTank = 10 #[m]
NTanksMax = 50 #[Number of Tanks]
#Air
Tinf =12.5 + 273.15 #[K]
#Cost Penalty for radius
CostPenalty = 1500 #[$/m^3]
# Time
T = 12 #Total time steps
# Aveage Heat Loss per hour
#AverageHeatLoss = 10*10^6/12



############CALCULATING PARAMETERS#####################################
#Mass considering the amount of usable energy
MassSalt = EnergyStorage/(HeatCapacitySalt*(T0Initial-T0Min)) #[kg]
InitialHeatStorage = MassSalt*HeatCapacitySalt*(T0Initial)#[kJ]
EnergyLossSalt = HeatCapacitySalt*MassSalt*(T0Initial-T0)

#### INITIALIZE MODEL ######################################################################
# using the knitro solver
ampl_solver = "snopt"
m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,ampl_solver), ["outlev=2"]))

###### DECISION VARIABLES #########################
# Number of tanks
@variable(m, NTanks >= 1,start=50) #[#]
# Internal redius of tanks
@variable(m, rTankIn >= 0,start=5)#[m]
# Thickness of Steel
@variable(m, LSteel >= 0,start=.1) #[m]
# Thickness of Insulation
@variable(m, LInsulation >=0,start=.3) #[m]
# Salt temperature
#@variable(m, T0 >=131+273.15) #[K]
# Steel temperature
@variable(m, T1 >=0) #[K]
# Insulator temperature
@variable(m, T2 >=0) #[K]

###### NONLINEAR EXPRESSIONS USING DECISON VARIABLES ##################################################################################################################################################################
# Steel Volume [m3]
@NLexpression(m, VolumeSteel, NTanks*pi*((2*LSteel+HeightTank)*rTankIn^2+HeightTank*(rTankIn+LSteel)^2))
# Insulation Volume [m^3]
@NLexpression(m, VolumeInsulation, NTanks*pi*(2*LInsulation*rTankIn^2+HeightTank*((rTankIn+LSteel+LInsulation)^2-(rTankIn+LSteel)^2)))

@NLexpression(m,VolumeTanks, NTanks*pi*(rTankIn^2)*HeightTank)
#Salt Energy Loss
#@NLexpression(m, EnergyLossSalt, HeatCapacitySalt*MassSalt*(T0Initial-T0))
#Steel Conductive Energy Loss
@NLexpression(m, EnergyLossSteel, NTanks*2*pi*3600*ThermalConductivitySTeel*(T0-T1)*((rTankIn^2)/LSteel + HeightTank/(log((rTankIn + LSteel)/rTankIn))))
#Insulation Conductive Energy Loss
@NLexpression(m, ConductiveLossInsulation, NTanks*(T1-T2)*2*pi*3600*ThermalConductivityInsulation*((rTankIn^2)/LInsulation + HeightTank/(log((rTankIn + LSteel + LInsulation)/(rTankIn + LSteel)))))

#Insulation Convective Energy Loss
@NLexpression(m, ConvectiveLossInsulation, (T2-Tinf)*NTanks*2*pi*(rTankIn^2+ HeightTank*(rTankIn+LSteel+LInsulation))*ConvectiveHeatTransferInsulation*3600)
#Insulation Radiative Energy Loss
@NLexpression(m, RadiativeLossInsulation, (T2^4-Tinf^4)*NTanks*2*pi*(rTankIn^2+ HeightTank*(rTankIn+LSteel + LInsulation))*StefanBoltzmanCosntant*3600)

##### CONSTRAINTS: GENERAL ########################
# Salt temperature
#@constraint(m, T0 <= T0Initial)
#@constraint(m, T0 >= T0Min)
# Steel Thickness
@constraint(m, LSteel >= 0.02) #[m]
@constraint(m, LSteel <= 0.1) #[m]
# Insulation Thickness
@constraint(m, LInsulation >= 0) #[m]
@constraint(m, LInsulation <= 0.3) #[m]
# Internal Tank radius
@constraint(m, rTankIn >= 0.125) #[m]
@constraint(m, rTankIn <= 5) #[m]
# Number of Tanks
@constraint(m, NTanks <= NTanksMax) #Number of Tanks


##### CONSTRAINTS: ENERGY BALANCE ###############################################
#Salt Energy Loss = Conductive Energy through Steel
@NLconstraint(m, EnergyLossSalt == EnergyLossSteel)
# Conductive Energy through Steel  = Conductive Energy through Insulation
@NLconstraint(m, EnergyLossSteel == ConductiveLossInsulation)
# Insulation Energy Loss = Steel Energy Loss
@NLconstraint(m, ConductiveLossInsulation == ConvectiveLossInsulation + RadiativeLossInsulation)
#Constraint of Tank Volume
@NLconstraint(m, VolumeTanks == MassSalt/DensitySalt)

##### CONSTRAINTS: SALT ENERGY LOSS ################
# Salt Energy Loss <= # Maximum Affordable Energy Loss per hour
#@NLconstraint(m, EnergyLossSalt*12 <= MaxEnergyLoss)

#### Objective ###########################################################################################################
# Minimize the total costs of Thermal Energy Storage
@NLobjective(m, Min, CostSteel*DensitySteel*VolumeSteel + CostInsulation*VolumeInsulation + CostPenalty*rTankIn^3)

# Solve the model
solve(m)

#DeltaT0 = T0Initial-getvalue(T0)
#DeltaT2 = getvalue(T1)-getvalue(T2)

#### PRINTING RESULTS ####
println("\n")
println("\n")
println("##############################################\n")
println("OPTIMIZATION OF THERMAL ENERGY STORAGE\n")
println("##############################################\n")
println("\n")
println("Total Cost: ", getobjectivevalue(m))
println("\n")
println("Number of tanks: ", getvalue(NTanks))
println("Internal Tank Radius[m]: ", getvalue(rTankIn))
println("Steel Thickness[m]: ", getvalue(LSteel))
println("Insulation Thickness[m]: ", getvalue(LInsulation))
println("\n")
println("Temperature Profile\n")
println("Time ","-","       T0       \n")
for i=1:12
  print("hour= ",i,"   -   ",T0Initial*((1-0.004182)^i)," \n")
end
println("---")
