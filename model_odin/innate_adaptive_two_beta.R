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

# equations
deriv(T1[1:2]) <- -max(0, beta1[i] * T1[i] * V) - max(0, phi * T1[i] * F1) + max(0, rho * R[i])
deriv(R[1:2]) <- max(0, phi * T1[i] * F1) - max(0, rho * R[i])
deriv(I1[1:2]) <- max(0, beta1[i] * T1[i] * V) - max(0, delta * I1[i])
deriv(V) <- max(0, sum(pI)) - max(0, (c1 + k1 * A) * V)
deriv(F1) <- max(0, sum(qI)) - max(0, d * F1)
deriv(A) <- m * A * (1 - A / A_max)

dim(T1) <- 2
dim(I1) <- 2
dim(R) <- 2
dim(beta1) <- 2
dim(pI) <- 2
dim(qI) <- 2

# parameter values
beta1[] <- user()
delta <- user()
p <- user()
q <- user()
phi <- user()
rho <- user()
c1 <- user()
d <- user()
m <- user()
k1 <- user()
alpha <- user()
T_0 <- user()
V_0 <- user()
A_0 <- user()
A_max <- user()
pI[] <- p * I1[i]
qI[] <- q * I1[i]