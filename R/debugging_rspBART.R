# Just running the default values so I can go through the function and debug all
#the things

mlbench.friedman1.nointercation <- function (n, sd = 1)
{
  x <- matrix(runif(4 * n), ncol = 4)
  y <- 10 * sin(pi * x[, 1])
  y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1 <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1.break <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)
  y[x<0] <- y[x<0] + 5
  y[x>=0] <- y[x>=0] - 5

  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

library(mlbench)
n_ <- 101
# sim_train <- mlbench.friedman1.nointercation(n = n_,sd = 1)  |> as.data.frame()
# x_test <- mlbench.friedman1.nointercation(n = n_,sd = 1) |> as.data.frame() |> dplyr::select(dplyr::starts_with("x"))

sim_train <- mlbench.d1.break(n = n_,sd = 0.5)  |> as.data.frame()
x_test <- mlbench.d1.break(n = n_,sd = 0.5) |> as.data.frame() |> dplyr::select(dplyr::starts_with("x"))

x_train <- sim_train |> dplyr::select(dplyr::starts_with("x"))
y_train <- sim_train$y

# x_train <- x_train[,1:5]
# x_test <- x_test[,1:5]
n_tree <- 10
node_min_size = 1
n_mcmc = 2000
n_burn = 500
alpha = 0.95
beta = 2
df = 3
sigquant = 0.9
kappa = 2
tau = 100
scale_bool = TRUE
stump = FALSE
no_rotation_bool = FALSE
numcut = 100L # Defining the grid of split rules
usequants = TRUE
delta <- 1
source("R/other_functions.R")

# Splines parameters
nIknots = 5
dif_order = 0

