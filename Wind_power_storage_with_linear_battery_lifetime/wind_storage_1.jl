
#REVISED SIMPLEX METHOD

# Initialize JuMP to allow mathematical programming models
using JuMP

# Initialize MILP solver Clp
using Clp

###################################################################
#ENTERING EQUATION AND RESTRICTIONS
#################################################################
#Objective Function
c = [-3 -4 -5]
#Subject to:
A = [1 2 2.5;
     1 2 1.33;
     4 2 5;
     2 3 1]
#RHS
b = [14; 12; 14; 6]

#Variables name
Xnames = ["x1" "x2" "x3" "s1" "s2" "s3" "s4"]
#############

############## Completing the matrices of Simplex method
#Saving number of restrictions
In = size(A,1)

#Augmented matrix
Eye = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1]

B = zeros(size(A))
BV = zeros(In)
Xsol = [A Eye b
        c zeros(1,In) 0]
########################################################


using LinearAlgebra

# Input Parameters
A=[1 2 2.5;1 2 1.33;4 2 5;2 3 1]
c=[-3 -4 -5]
c_original=[-3 -4 -5]
b=[14 12 14 6]'
c_B=[0 0 0 0]
B=Matrix{Float64}(I, 4, 4)

#While loop until NBV coefficients are non-negaive
while min(c_original-c_B*inv(B)*A...)< -10^-8

# Finding new entering basic variable from updates top row of tableau
index_Col=findmin([c_original-c_B*inv(B)*A]...)[2][2]

# Finding leaving variable
index_Row=findmin(b./A[:,index_Col])[2][1]

# Replacing leaving variable column with entering variable column from A Matrix
B[:,index_Row]=A[:,index_Col]

# Computing x_B: Current Solution
x_B=inv(B)*b

# Updating c and c_B vectors as variables enter and leave basis
c_B[index_Row]=c_original[index_Col]

#Computing objective function value
minus_f=-c_B*inv(B)*b

# Printing first B, x_B and -f
println("B Matrix: "); display(B)
println("x_B Vector: "); display(x_B)
println("- f(x): "); display(minus_f)
display([c_original-c_B*inv(B)*A -c_B*inv(B)])
end
