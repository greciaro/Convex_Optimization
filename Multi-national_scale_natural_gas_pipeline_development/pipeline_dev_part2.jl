
####################################
######### Initialize tools #########
####################################


# Initialize JuMP to allow mathematical programming models
using JuMP

# Initialize MILP solver Clp
using Clp


###############################
######### Define Sets #########
###############################

# Set of countries
COUNTRIES = ["Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador","Peru","Trinidad&Tobago","Uruguay","Venezuela"]
nCOUNTRIES = length(COUNTRIES)


# Set of time steps
T = 17 #From 2009 to 2025
TIME = 1:T

###################################################
############ Define parameters and data ###########
###################################################

# Natural Gas production per country [x10^9 m^3]
InitialProduction = [41.37, 12.63, 10.28, 1.36, 10.48, 0.28, 3.48, 40.61, 0, 18.43]

# Natural Gas consumption per country [x10^9 m^3]]
InitalConsumption = [43.13, 2.83, 18.72, 2.83, 8.69, 0.28, 3.48, 20.87, 0.04, 20.22]

# Natural Gas production growth per country [%/yr]
ProductionGrowth = [2, 3, 3, 2, -1, 3, 5, -3, 0, 7]

# Natural Gas consumption growth per country [%/yr]
ConsumptionGrowth = [3, 2, 2, 3, 3, 3, 2, 4, 2, 2]

# Distances between South American capitol cities. Indexed n for columns, m for rows.
Distance =  [0      1891.6 2486.1 1111.2 4643.6 4219.2 3116.4 3129.0 230.4 5084.6;
			  1891.6 0      1489.6 1598.4 1197.8 1175.8 775.5  1505.6 1775.7 1466.6;
			  2486.1 1489.6 0      3249.5 3995.9 4180.5 3541.3 2403.3 2378.3 3862.9;
			  1111.2 1598.4 3249.5 0      4243.2 3606.3 2464.6 3103.8 1340.9 4899.6;
			  4643.6 1197.8 3995.9 4243.2 0       991.1 1878.9  955.7 4767.8 1025.9;
			  4219.2 1175.8 4180.5 3606.3  991.1 0      1142.1 1542.2 4382.6 2010.6;
			  3116.4  775.5 3541.3 2464.6 1878.9 1142.1 0      1521.5 3293.9 2744.0;
			  3129.0 1505.6 2403.3 3103.8  955.7 1542.2 1521.5 0      3161.5  590.1;
			  230.4  1775.7 2378.3 1340.9 4767.8 4382.6 3293.9 3161.5 0      5165.3;
			  5084.6 1466.6 3862.9 4899.6 1025.9 2010.6 2744.0  590.1 5165.3 0     ;]

Production = zeros(T,nCOUNTRIES)
Consumption = zeros(T,nCOUNTRIES)
AmountGasTransported = zeros(T,nCOUNTRIES,nCOUNTRIES)
CostTransport = zeros(T,nCOUNTRIES,nCOUNTRIES)
n = 0
nn = 0


for t=1:T
	for n=1:nCOUNTRIES
			Production[t,n] = InitialProduction[n]*exp(ProductionGrowth[n]*(t-(t-1)))
			Consumption[t,n] = InitalConsumption[n]*exp(ConsumptionGrowth[n]*(t-(t-1)))
			for nn =1:nCOUNTRIES
		    CostTransport[t,n,nn] = (0.02/1000)*Distance[n,nn]
		end
end
end

#for t=1:T
t=1

			  # Define the model name and solver. In this case, model name is "nn"
			  m = Model(solver=ClpSolver())

			  ######## Decision variable ########
			  @variable(m, AmountGasTransported[1:T,1:nCOUNTRIES,1:nCOUNTRIES] >= 0)


			  ######## Objective Functions #########
			  # Single objective for minimizing cost
			  @objective(m, Min, sum(sum(CostTransport[t,n,nn]*AmountGasTransported[t,n,nn] for nn=1:nCOUNTRIES) for n=1:nCOUNTRIES))


			  ############# Constraints ############
			  # Production in each country constraint
			  @constraint(m, sum(Production[t,n] for n=1:nCOUNTRIES) == sum(AmountGasTransported[t,n,nn] for nn=1:nCOUNTRIES))
			  @constraint(m, sum(Consumption[t,nn] for nn=1:nCOUNTRIES) == sum(AmountGasTransported[t,n,nn] for n=1:nCOUNTRIES))
#end
solve(m)
print(m)

println("Objective value: ", getobjectivevalue(m))
println("Transport = ",getvalue(AmountGasTransported))
