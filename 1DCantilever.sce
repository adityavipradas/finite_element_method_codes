//FEM Solver for a static cantilever subjected to axial tensile forces using linear 1D elements
//Displacement at node 1 is always 0

function[] = plot_format()
    //Get the handle of current axes
    g = gca()
    //Give labels and set label properties
    g.labels_font_color=5
    g.font_size=3
    g.grid=[1,1]
    g.box="off"
endfunction

clc
A = input("Enter the uniform cross-sectional area of the cantilever(mm2): ")
E = input("Enter the Youngs Modulus of the material used(N/mm2): ")
L = input("Enter the length of the cantilever(mm): ")
n = input("Enter the number of 1D linear elements to be used: ")
F = zeros(n+1, 1)
K = zeros(n+1, n+1)
Q = zeros(1, 6*n)
Z = zeros(1, 6*n)
l = 1
m = 1
stiff = (A*E)/(L/n)
Kele = [stiff, (-1 * stiff); (-1 * stiff), stiff]

//Form global force vector
for i = 1:1:(n+1)
    printf("\nFor node %d:\n",i)
    F(i, 1) = input("Enter the force on node in newton: ")
end

//Direct Stiffness Method
for i = 1:1:n
    for j = i:1:(i+1)
        for k = i:1:(i+1)
            K(j, k) = K(j, k) + Kele(l, m)
            m = 2
        end
        m = 1
        l = 2
    end
    l = 1
    m = 1
end
//Using Penalty Approach to solve the equation
Kmax = max(K)
C = Kmax *(10^4)
Kp = K
Kp(1, 1) = K(1, 1) + C
U = inv(Kp) * F
R = -1 * C * U(1, 1)
//Displacements at the intermediate points
j = 1
k = 0
for i = 1:1:n
    for zeta = 0:0.2:1
        q = (zeta * U((i+1), 1)) + ((1 - zeta) * U(i, 1))
        strain = (U((i+1), 1) - U(i, 1))/(L/n)
        stress = E * strain
        Z(1, j) = (k*(L/n)) + zeta*(L/n)
        Q(1, j) = q
        Strain(1, j) = strain
        Stress(1, j) = stress
        j = j + 1
    end
    k = k + 1
end
funcprot(0)
subplot(2,2,1)
plot(Z, Q, '-*')
title('Displacement against Length')
xlabel('Cantilever Length')
ylabel('Displacement')
plot_format()
subplot(2,2,2)
plot(Z, Strain, '-*')
title('Strain against Length')
xlabel('Cantilever Length')
ylabel('Strain')
plot_format()
subplot(2,2,3)
plot(Z, Stress, '-*')
title('Stress against Length')
xlabel('Cantilever Length')
ylabel('Stress')
plot_format()
//Print Results
printf("Global stiffness matrix:\n")
disp(K)
printf("\nUsing penalty approach\n")
printf("\nDisplacement vector(mm):\n")
disp(U)
printf("Reaction at the rigid support(N):\n")
disp(R)
