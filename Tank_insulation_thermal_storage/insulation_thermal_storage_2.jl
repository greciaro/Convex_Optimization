
using JuMP
using AmplNLWriter


###### PARAMETERS and DATA ###################################
#Salt
T0Initial = 520 + 273.15 #[K]
T0Min = 131 + 273.15 #[K]
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
HeightTank = 3 #[m]
NTanksMax = 50 #[Number of Tanks]
#Air
Tinf =12.5 + 273.15 #[K]
#Cost Penalty for radius
CostPenalty = 10000 #[$/m^3]
#Requirements
EnergyStorage = 100*10^6 #kJ
# Time
T = 12 #Total time steps

############CALCULATING PARAMETERS#####################################
#Heat Storage
MassSalt = EnergyStorage/(HeatCapacitySalt*(T0Initial-T0Min)) #[kg]
InitialHeatStorage = MassSalt*HeatCapacitySalt*(T0Initial)#[kJ]
HeatStorageLimit = 0.90*EnergyStorage #[kJ]


#### INITIALIZE MODEL ######################################################################
# using the knitro solver
ampl_solver = "snopt"
m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,ampl_solver), ["outlev=2"]))

###### DECISION VARIABLES #########################
# Number of tanks
@variable(m, NTanks >=1,start=50) #[#]
# Internal redius of tanks
@variable(m, rTankIn >=0, start=2.5) #[m]
# Thickness of Steel
@variable(m, LSteel >=0,start=0.5) #[m]
# Thickness of Insulation
@variable(m, LInsulation >=0,start=0.15) #[m]
# Salt temperature
@variable(m, T0[1:T+1] >=0) #[K]
# Steel temperature
@variable(m, T1[1:T] >=0) #[K]
# Insulator temperature
@variable(m, T2[1:T] >=0) #[K]

###### NONLINEAR EXPRESSIONS USING DECISON VARIABLES ##################################################################################################################################################################
# Steel Volume [m3]
@NLexpression(m, VolumeSteel, NTanks*pi*((2*LSteel+HeightTank)*rTankIn^2+HeightTank*(rTankIn+LSteel)^2))
# Insulation Volume [m^3]
@NLexpression(m, VolumeInsulation, NTanks*pi*(2*LInsulation*rTankIn^2+HeightTank*((rTankIn+LSteel+LInsulation)^2-(rTankIn+LSteel)^2)))

#Salt Energy Loss
@NLexpression(m, EnergyLossSalt[t=1:T], HeatCapacitySalt*MassSalt*(T0[t]-T0[t+1]))
#Steel Conductive Energy Loss
@NLexpression(m, EnergyLossSteel[t=1:T], NTanks*2*pi*3600*ThermalConductivitySTeel*(T0[t+1]-T1[t])*((rTankIn^2)/LSteel + HeightTank/(log((rTankIn + LSteel)/rTankIn))))
#Insulation Conductive Energy Loss
@NLexpression(m, ConductiveLossInsulation[t=1:T], NTanks*(T1[t]-T2[t])*2*pi*3600*ThermalConductivityInsulation*((rTankIn^2)/LInsulation + HeightTank/(log((rTankIn + LSteel + LInsulation)/(rTankIn + LSteel)))))

#Insulation Convective Energy Loss
@NLexpression(m, ConvectiveLossInsulation[t=1:T], (T2[t]-Tinf)*NTanks*2*pi*(rTankIn^2+ HeightTank*(rTankIn+LSteel+LInsulation))*ConvectiveHeatTransferInsulation*3600)
#Insulation Radiative Energy Loss
@NLexpression(m, RadiativeLossInsulation[t=1:T], (T2[t]^4-Tinf^4)*NTanks*2*pi*(rTankIn^2+ HeightTank*(rTankIn+LSteel+LInsulation))*StefanBoltzmanCosntant*3600)

#Heat Storage Balance
@NLexpression(m, HeatinSalt[t=1:T+1], HeatCapacitySalt*MassSalt*T0[t])

##### CONSTRAINTS: GENERAL ########################
#Heat Balance
@NLconstraint(m,[t=1:T], HeatinSalt[t+1] == HeatinSalt[t]-EnergyLossSalt[t])
#@NLconstraint(m,HeatinSalt[13] >= HeatStorageLimit)
# Salt temperature
@constraint(m, T0[1] == T0Initial)
@constraint(m, [t=1:T], T0[t+1] <= T0Initial)
#@constraint(m, [t=1:T], T0[t+1] >= T0Min)
# Steel Thickness
@constraint(m, LSteel >= 0.02) #[m]
@constraint(m, LSteel <= 0.1) #[m]
# Insulation Thickness
@constraint(m, LInsulation >= 0) #[m]
@constraint(m, LInsulation <= 0.3) #[m]
# Internal Tank radius
@constraint(m, rTankIn >= 0.125) #[m]
@constraint(m, rTankIn <= 2.5) #[m]
# Number of Tanks
@constraint(m, NTanks <= NTanksMax) #Number of Tanks


##### CONSTRAINTS: ENERGY BALANCE ###############################################
# Salt Energy Loss = Conductive Energy through Steel
@NLconstraint(m, [t=1:T], EnergyLossSalt[t] <= EnergyLossSteel[t])
# Conductive Energy through Steel  = Conductive Energy through Insulation
@NLconstraint(m, [t=1:T], EnergyLossSteel[t] <= ConductiveLossInsulation[t])
# Insulation Energy Loss = Steel Energy Loss
@NLconstraint(m, [t=1:T], ConductiveLossInsulation[t] <= ConvectiveLossInsulation[t] + RadiativeLossInsulation[t])

##### CONSTRAINTS: SALT ENERGY LOSS ################
# Salt Energy Loss <= # Maximum Affordable Energy Loss per hour
#@NLconstraint(m, EnergyLossSalt <= MaxEnergyLoss)

#### Objective ###########################################################################################################
# Minimize the total costs of Thermal Energy Storage
@NLobjective(m, Min, CostSteel*DensitySteel*VolumeSteel + CostInsulation*VolumeInsulation + CostPenalty*rTankIn^3)

# Solve the model
solve(m)


#### PRINTING RESULTS ####
println("##############################################\n")
println("OPTIMIZATION OF THERMAL ENERGY STORAGE\n")
println("##############################################\n")
println("\n")
println("Total Cost: ", getobjectivevalue(m))
println("\n")
println("Number of tanks: ", getvalue(NTanks))
println("Internal Tank Radius[m]: ", getvalue(rTankIn))
println("Steel Thicness[m]: ", getvalue(LSteel))
println("Insulation Thicness[m]: ", getvalue(LInsulation))
println("\n")
println("Emperature Profiles\n")
println("Time ","-","       T0       ","-","       T1       ","-","      T2      ","\n")
for i=1:12
  print("hour= ",i,"   -   ",getvalue(T0[i]),"   -   ",getvalue(T1[i]),"   -   ",getvalue(T2[i]),"\n")
end
println("---")
