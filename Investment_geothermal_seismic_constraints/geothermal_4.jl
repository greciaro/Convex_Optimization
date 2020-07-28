# Initialize MILP solver Clp
using Cbc
using CSV
using Plots
using DataFrames

Population = CSV.read("Population.csv")[1:end,1] # []people per block|
FaultDistance = CSV.read("FaultDistance.csv")[1:end, 1] # [km]
GeothermalResource = CSV.read("GeothermalResource.csv")[1:end, 1] # [MW/km^2]
TransmissionDistance = CSV.read("TransmissionDistance.csv")[1:end, 1] # [km]
DrillingRightsCost = CSV.read("DrillingRightsCost.csv")[1:end, 1] # [$/block]
EGSInstalledCost = 3*10^6 #[$/MW]
TransmissionCost = 1*10^6 #[$/km]
Area = 100 #km^2 per block
#####################################
#Population Matrix
################################
GeothermalResourceMatrix = zeros(112,44)
FaultDistanceMatrix = zeros(112,44)
TransmissionDistanceMatrix = zeros(112,44)
DrillingRightsCostMatrix = zeros(112,44)
k = 1
for j = 1:44
    for i=1:112
      GeothermalResourceMatrix[i,j] = GeothermalResource[k]
      FaultDistanceMatrix[i,j] = FaultDistance[k]
      TransmissionDistanceMatrix[i,j] = TransmissionDistance[k]
      DrillingRightsCostMatrix[i,j] = DrillingRightsCost[k]
    global k = k+1
    end
end

k=1
ExtraMatrix = zeros(112+20,44+20)
for j= 11:54
for i = 11:122
ExtraMatrix[i,j] = Population[k]
global k = k+1
end
end

SumationPopulation = zeros(112,44)
for j = 1:44
for i = 1:112
SumationPopulation[i,j] =sum(sum(ExtraMatrix[a,b] for a =i:i+20) for b =j:j+20)
end
end
###########################################################

########## Declare model  ##########
# Define the model name and solver. In this case, model name is "m"
m = Model(solver=CbcSolver())

######## Decision variables ########
@variable(m, EGS[1:112,1:44], Bin)

######## Objective Function #########

@objective(m, Min, sum(sum(EGS[i,j]*(EGSInstalledCost*GeothermalResourceMatrix[i,j]*Area + TransmissionDistanceMatrix[i,j]*TransmissionCost + DrillingRightsCostMatrix[i,j]) for i =1:112) for j=1:44))

############# Constraints ############
@constraint(m, [i=1:112,j=1:44], sum(sum(EGS[i,j]*GeothermalResourceMatrix[i,j]  for i =1:112) for j=1:44) >= 100)

for j = 1:44
for i = 1:112
if FaultDistanceMatrix[i,j] <= 100
   @constraint(m, EGS[i,j] == 0)
end
end
end
for j = 1:44
for i = 1:112
if SumationPopulation[i,j] >= 100000
   @constraint(m, EGS[i,j] == 0)
end
end
end

########### Print and solve ##########
status = solve(m)
# Print more detailed results to screen

#println("Energy storage size: ",StorageLimit[s])
println("Minimum cost: ", getobjectivevalue(m))
#println(s)
#Optimal_profit[1,s] = getobjectivevalue(m)
#getvalue(EGS)
x = findall(x-> x>=1, getvalue(EGS))

heatmap(getvalue(EGS))
#println("Yes")

#for j = 1:44
#for i = 1:112
#totalEGS
#end
#end
