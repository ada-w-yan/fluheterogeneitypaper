# the all-or-nothing susceptibility model with two compartments

# initial conditions
initial(T1[1]) <- alpha * T_0
initial(T1[2]) <- (1 - alpha) * T_0
initial(R[1]) <- 0
initial(R[2]) <- 0
initial(I1[1]) <- 0
initial(I1[2]) <- 0
initial(V) <- V_0
initial(F1) <- 0
initial(A) <- A_0
initial(AUC) <- 0
initial(AUC_log10) <- 0
initial(dose_response1) <- 0
initial(dose_response2) <- 0

# equations
deriv(T1[1:2]) <- -max(0, beta1[i] * T1[i] * V) - max(0, phi[i] * T1[i] * F1) + max(0, rho[i] * R[i])
deriv(R[1:2]) <- max(0, phi[i] * T1[i] * F1) - max(0, rho[i] * R[i])
deriv(I1[1:2]) <- max(0, beta1[i] * T1[i] * V) - max(0, delta * I1[i])
deriv(V) <- max(0, sum(pI)) - max(0, (c1 + k1 * A) * V)
deriv(F1) <- max(0, sum(qI)) - max(0, d * F1)
deriv(A) <- m * A * (1 - A / A_max)
deriv(AUC) <- V
deriv(AUC_log10) <- max(0, log10(V))
deriv(dose_response1) <- max(0, 1 - exp(-gamma1 * V))
deriv(dose_response2) <- max(0, 1 - exp(-gamma2 * V))

dim(T1) <- 2
dim(I1) <- 2
dim(R) <- 2
dim(beta1) <- 2
dim(p) <- 2
dim(q) <- 2
dim(phi) <- 2
dim(rho) <- 2
dim(pI) <- 2
dim(qI) <- 2

# parameter values
beta1[] <- user()
delta <- user()
p[] <- user()
q[] <- user()
phi[] <- user()
rho[] <- user()
c1 <- user()
d <- user()
m <- user()
k1 <- user()
alpha <- user()
T_0 <- user()
V_0 <- user()
A_0 <- user()
A_max <- user()
gamma1 <- user()
gamma2 <- user()
pI[] <- p[i] * I1[i]
qI[] <- q[i] * I1[i]