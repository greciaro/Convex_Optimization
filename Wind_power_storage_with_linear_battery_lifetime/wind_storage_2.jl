#To add any package
#Pkg.add("CSV")
# Initialize JuMP to allow mathematical programming models
using JuMP
# Initialize MILP solver Clp
using Clp
using Plots
using CSV
using DataFrames



# Set of time steps
T = 432
TIME = 1:T

############ Define parameters and data ###########

AvailablePower = CSV.read("Power.csv")[1:end,1]./6 # [MWh per time step]
PriceofPower = CSV.read("PriceofPower.csv")[1:end, 1] # [$ per time step]

Ccap = 250000 #[$/MWhcapacity]
Lc = 3000 #[MWh lifetime stored in battery/MWh of capacity]
StorageLimit = 200 #[MWh]
InitialStorage = 0 #[MWh]
RoundtripLosses = 0.0645 #[MWh/MWh stored]

#ConstantCost = ((1-RoundtripLosses)*Ccap)/(Lc*StorageLimit)

println("Amount of timesteps = ",size(AvailablePower))
#println("Constant Cost= ",ConstantCost)
#println("Production Tirn",ProductionTrin)

########## Declare model  ##########
# Define the model name and solver. In this case, model name is "m"
m = Model(solver=ClpSolver())

######## Decision variables ########
@variable(m, StorePowerIn[1:T] >= 0)
@variable(m, StorePowerOut[1:T] >= 0)
@variable(m, InStorage[1:T] >= 0)

######## Objective Function #########

@objective(m, Max, sum((AvailablePower[t]-StorePowerIn[t])*PriceofPower[t] for t=1:T) + sum(StorePowerOut[t]*(1-RoundtripLosses)*PriceofPower[t] for t=1:T) - (StorageLimit*Ccap*sum(StorePowerOut[t] for t=1:T)/StorageLimit/Lc))

############# Constraints ############
@constraint(m, InStorage[1] == InitialStorage)
@constraint(m, [t=1:T-1], InStorage[t+1] == InStorage[t] + StorePowerIn[t+1] - StorePowerOut[t+1])
@constraint(m, [t=1:T], InStorage[t] <= StorageLimit)
@constraint(m, [t=1:T], StorePowerIn[t] <= AvailablePower[t])
#@constraint(m, [t=1:T], StorePowerOut[t] <= 10)

########### Print and solve ##########
status = solve(m)
# Print more detailed results to screen
println("Mazimum profit: ", getobjectivevalue(m))
InStorage_optimal = getvalue(InStorage)
println("In Storage = ", InStorage_optimal)
#p1 = plot(TIME, InStorage_optimal)
#p2 = plot(TIME, PriceofPower)
plot(TIME, InStorage_optimal,seriestype=:line, ylabel="Power Stored in the Battery [MWh]",xlabel="Time step", label="")
plot(TIME, PriceofPower,seriestype=:line, ylabel="Price of Power per MWh",xlabel="Time step", label="")
