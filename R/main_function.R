rm(list=ls())
source("R/debugging_rspBART.R")
source("R/tree_functions.R")
set.seed(42)

# Creating the main function from the rspBART
rspBART <- function(x_train,
                    y_train,
                    x_test,
                    node_min_size = 15,
                    n_mcmc = 2000,
                    n_burn = 500,
                    alpha = 0.95,
                    beta = 2,
                    df = 3,
                    sigquant = 0.9,
                    kappa = 2,
                    # Splines parameters
                    nIknots = 3,
                    dif_order = 2,
                    tau = 100,
                    scale_bool = TRUE,
                    stump = FALSE,
                    no_rotation_bool = FALSE,
                    numcut = 100L, # Defining the grid of split rules
                    usequants = FALSE
) {

  # Verifying if x_train and x_test are matrices
  if(!is.data.frame(x_train) || !is.data.frame(x_test)){
    stop("Insert valid data.frame for both data and xnew.")
  }


  # Getting the valid
  dummy_x <- base_dummyVars(x_train)

  # Create a list
  if(length(dummy_x$facVars)!=0){
    for(i in 1:length(dummy_x$facVars)){
      # See if the levels of the test and train matches
      if(!all(levels(x_train[[dummy_x$facVars[i]]])==levels(x_test[[dummy_x$facVars[i]]]))){
        levels(x_test[[dummy_x$facVars[[i]]]]) <- levels(x_train[[dummy_x$facVars[[i]]]])
      }
      df_aux <- data.frame( x = x_train[,dummy_x$facVars[i]],y)
      formula_aux <- stats::aggregate(y~x,df_aux,mean)
      formula_aux$y <- rank(formula_aux$y)
      x_train[[dummy_x$facVars[i]]] <- as.numeric(factor(x_train[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

      # Doing the same for the test set
      x_test[[dummy_x$facVars[i]]] <- as.numeric(factor(x_test[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

    }
  }

  # Getting the train and test set
  x_train_scale <- as.matrix(x_train)
  x_test_scale <- as.matrix(x_test)

  # Scaling x
  x_min <- apply(as.matrix(x_train_scale),2,min)
  x_max <- apply(as.matrix(x_train_scale),2,max)

  # Storing the original
  x_train_original <- x_train
  x_test_original <- x_test


  # Normalising all the columns
  for(i in 1:ncol(x_train)){
    x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
    x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
  }



  # Creating the numcuts matrix of splitting rules
  xcut_m <- matrix(NA,nrow = numcut,ncol = ncol(x_train_scale))
  for(i in 1:ncol(x_train_scale)){

    if(nrow(x_train_scale)<numcut){
      xcut_m[,i] <- sort(x_train_scale[,i])
    } else {
      xcut_m[,i] <- seq(min(x_train_scale[,i]),
                        max(x_train_scale[,i]),
                        length.out = numcut+2)[-c(1,numcut+2)]
    }
  }




  # =========================================================================================================
  # Getting the Splines Basis functions
  # =========================================================================================================

  knots <- apply(x_train_scale,
                 2,
                 function(x){seq(min(x),max(x),length.out = nIknots)})


  D_train <- matrix(NA,
                       nrow = nrow(x_train_scale),
                       ncol = (nrow(knots)+3)*length(dummy_x$continuousVars))

  D_test <- matrix(NA,
                    nrow = nrow(x_test_scale),
                    ncol = (nrow(knots)+3)*length(dummy_x$continuousVars))

  # Setting new parameters for the spline
  ndx <- nIknots
  dx <- 1/ndx
  # New_knots
  new_knots <- matrix()
  new_knots <- matrix(mapply(x_min,x_max, FUN = function(MIN,MAX){seq(from = 0-3*dx, to = 1+3*dx, by = dx)}), ncol = length(dummy_x$continuousVars)) # MIN and MAX are 0 and 1 respectively, because of the scale
  colnames(new_knots) <- dummy_x$continuousVars

  basis_size <- (nrow(knots)+3) # Change this value to the desired size of each sublist
  D_seq <- 1:ncol(D_train)  # Replace this with the columns of D

  # Creating a vector
  basis_subindex <- split(D_seq, rep(1:(length(D_seq) %/% basis_size), each = basis_size, length.out = length(D_seq)))

  # Creating the natural B-spline for each predictor
  for(i in 1:length(basis_subindex)){
    B_train_obj <- splines::spline.des(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                       knots = new_knots[,dummy_x$continuousVars[i]],
                                       ord = 4,
                                       derivs = 0*x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = FALSE)$design
    # Returning to MOTR-BART
    D_train[,basis_subindex[[i]]] <- B_train_obj


    D_test[,basis_subindex[[i]]] <- splines::spline.des(x = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                           knots = new_knots[,dummy_x$continuousVars[i]],
                                           ord = 4,
                                           derivs = 0*x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design
  }

  # R-th difference order matrix
  # if(dif_order!=0){
  #   D <- D_gen(p = ncol(D_train),n_dif = dif_order)
  # } else {
  #   D <- diag(nrow = ncol(D_train))
  # }

  # Scaling the y
  min_y <- min(y_train)
  max_y <- max(y_train)

  # Getting the min and max for each column
  min_x <- apply(x_train_scale,2,min)
  max_x <- apply(x_train_scale, 2, max)

  # Scaling "y"
  if(scale_bool){
    y_scale <- normalize_bart(y = y_train,a = min_y,b = max_y)
    tau_gamma <- tau_mu <- (4*n_tree*(kappa^2))

  } else {
    y_scale <- y_train

    tau_gamma <- tau_mu <- (4*n_tree*(kappa^2))/((max_y-min_y)^2)
  }


  # Calculating the penalty matrix
  # diag(D) <- diag(D) + dif_order*tau_gamma
  # P <- crossprod(D)
  # P[1,1] <- P[1,1] + tau_gamma
  # P[nrow(P),ncol(P)] <- P[nrow(P),ncol(P)] + tau_gamma
  # P <- P + diag(tau_gamma,nrow = nrow(P))

  # Getting the naive sigma value
  nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

  # Calculating tau hyperparam
  a_tau <- df/2

  # Calculating lambda
  qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
  lambda <- (nsigma*nsigma*qchi)/df
  d_tau <- (lambda*df)/2


  # Call the bart function
  tau_init <- nsigma^(-2)

  mu_init <- mean(y_scale)

  # Creating the vector that stores all trees
  all_tree_post <- vector("list",length = round(n_mcmc-n_burn))


  # =====================================================================
  # ========= From here I gonna initialise the BART function itself =====
  # =====================================================================
  n_post <- (n_mcmc-n_burn)
  all_trees <- vector("list", n_mcmc)
  all_betas <- vector("list",n_mcmc)
  tau_beta <- n_tree # In this first scenario we are going to work with a single value of \tau
  all_tau_beta <- numeric(n_mcmc)
  all_tau_gamma <- numeric(n_mcmc)
  all_delta <- numeric(n_mcmc)
  all_tau <- numeric(n_mcmc)

  all_y_hat <- matrix(NA,nrow = n_mcmc,ncol = nrow(x_train_scale))
  all_trees_fit <- vector("list",n_mcmc)
  all_trees <- vector("list",n_mcmc)
  forest <- vector("list",n_tree)

  # Partial component pieces
  partial_train_fits <- vector("list", n_tree)

  proposal_outcomes <- setNames(data.frame(matrix(nrow = 0, ncol = 6)),
                                c("tree_number" , "proposal", "status","mcmc_iter", "new_tree_loglike", "old_tree_loglike"))
  all_train_indexes <- data.frame(matrix(data = NA,nrow = nrow(xcut_m),ncol = ncol(xcut_m)))

  # Gonna create a list of lists to store all the indexes for all split rules and cutpoints
  all_var_splits <- vector("list",ncol(x_train_scale))
  names(all_var_splits) <- colnames(x_train_scale)

  # Iterating over all possible x.columns
  for(i in 1:length(all_var_splits)){

    # Creating the dummy for a list of index to store all numeric split values
    all_cut_points <- vector("list", nrow(xcut_m))


    for(j in 1:length(all_cut_points)){

      # Getting the node indexes object
      left_train_list <- vector("list",length = 1L)
      names(left_train_list) <- "left_train"
      right_train_list <- vector("list",length = 1L)
      names(right_train_list) <- "right_train"
      left_test_list <- vector("list",length = 1L)
      names(left_test_list) <- "left_test"
      right_test_list <- vector("list",length = 1L)
      names(right_test_list) <- "right_test"

      node_index <- append(left_train_list, right_train_list) |>
                    append(left_test_list) |> append(right_test_list)

      all_cut_points[[j]]$left_train <-  which(x_train_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_train <-  which(x_train_scale[,i] >= xcut_m[j,i])
      all_cut_points[[j]]$left_test <-  which(x_train_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_test <-  which(x_train_scale[,i] >= xcut_m[j,i])

    }

    all_var_splits[[i]] <- all_cut_points

  }

  # Creating the "data" list which contain all elements necessary to run
  #most of the functions
  data <- list(x_train = x_train_scale,
               x_test = x_test_scale,
               y_train = y_scale,
               D_train = D_train,
               D_test = D_test,
               basis_subindex = basis_subindex,
               all_var_splits = all_var_splits,
               n_tree = n_tree,
               tau_mu = tau_mu,
               tau_gamma = tau_gamma,
               tau = tau,
               a_tau = a_tau,
               d_tau = d_tau,
               tau_beta = tau_beta,
               delta = delta,
               # P = P,
               node_min_size = node_min_size)

  #   So to simply interepret the element all_var_splits each element correspond
  #to each variable. Afterwards each element corresponds to a cutpoint; Finally,
  #inside that level we would have the index for the the left and right nodes;

  # Initialing for storing post samples
  post <- 0

  # Initialising all the stumps
  for(k in 1:data$n_tree){
    forest[[k]] <- stump(data = data)
  }

  # and tree predictions
  trees_fit <- matrix(0,nrow = n_tree,ncol = nrow(x_train_scale))

  for(i in 1:n_tree){
    trees_fit[i,] <- y_scale/n_tree
  }

  # Initialsing the loop
  for(i in 1:n_mcmc){

    # Initialising the partial train tree fits
    partial_train_fits[[i]] <- vector("list",data$n_tree)
    names(partial_train_fits[[i]]) <- paste0("tree",1:data$n_tree)

    # Initialising orogress bar
    progress <- i / n_mcmc * 100

    for(t in 1:data$n_tree){

        # Calculating the partial residuals
        if(n_tree>1){
          partial_residuals <- y_scale-colSums(trees_fit[-t,,drop = FALSE])
        } else {
          partial_residuals <- y_scale
        }


      # Sample a verb
      verb <- sample(c("grow","prune", "change"), prob = c(0.3,0.3,0.4),size = 1)

      # Forcing to grow when only have a stump
      if(length(forest[[t]])==1){
        verb <- "grow"
      }


      # Sampling a verb
      if(verb == "grow"){
        forest[[t]] <- grow(tree = forest[[t]],
                            curr_part_res = partial_residuals,
                            data = data)
      } else if (verb == "prune"){
        forest[[t]] <- prune(tree = forest[[t]],
                             curr_part_res = partial_residuals,
                             data = data)
      } else if (verb == "change"){
        forest[[t]] <- change(tree = forest[[t]],
                              curr_part_res = partial_residuals,
                              data = data)
      }

      # cat("Forest size:", (length(forest)),"\n")
      # cat("Forest verb:", verb,"\n")

      # Updating the intercept
      forest[[t]] <- updateGamma(tree = forest[[t]],
                                 curr_part_res = partial_residuals,
                                 data = data)

      # Updating the betas
      forest[[t]] <- updateBetas(tree = forest[[t]],
                                 curr_part_res = partial_residuals,
                                 data = data)

      # Getting the predictions
      tree_predictions <- getPredictions(tree = forest[[t]],
                                         data = data)

      trees_fit[t,] <- rowSums(tree_predictions$y_train_hat)

      partial_train_fits[[t]] <- tree_predictions$y_train_hat


    }
    # Getting final predcition
    y_hat <- colSums(trees_fit)

    # Seeing the results for the unidimensional cases.
    plot(x_train_scale,y_scale)
    for(plot_i in 1:n_tree){
      points(x_train_scale,trees_fit[plot_i,],pch=20,col = ggplot2::alpha(plot_i,0.2))
    }
    points(x_train_scale,y_hat,col = "blue",pch=20)



    # Updating all other parameters
    data$tau_beta <- update_tau_betas(forest = forest,data = data)
    # data$tau_gamma <- update_tau_gamma(forest = forest,data = data)

    # Updating delta
    # data$delta <- update_delta(data = data)


    # Getting tau
    data$tau <- update_tau(y_train_hat = y_hat,
                           data = data)


    # Storing all predictions
    all_trees[[i]] <- forest
    all_tau[[i]] <- data$tau
    all_tau_gamma[[i]] <- data$tau_gamma
    all_trees_fit[[i]] <- partial_train_fits
    all_y_hat[i,] <- y_hat
    all_tau_beta[i] <- data$tau_beta
    all_delta[i] <- data$delta


    # Print progress bar
    cat("\rProgress: [", paste(rep("=", floor(progress / 5)), collapse = ""),
        paste(rep(" ", floor((100 - progress) / 5)), collapse = ""),
        "] ", sprintf("%.1f%%", progress))

    # Flush the output
    flush.console()

    # Simulate some work
    Sys.sleep(0.1)

  }

  # ====== Few analyses from the results ======
  plot(all_tau,type = "l")
  y1_hat <- matrix(0,nrow = n_post,ncol = nrow(data$x_train))

  for(i in 1:nrow(y1_hat)){
    for(t in 1:n_tree){
      y1_hat[i,] <- all_trees_fit[[i]][[t]][,1]
    }
  }

  plot(x_train_scale[,1],colMeans(y1_hat))
  plot(all_tau_beta, type = "l")
  plot(all_tau_gamma,type = "l")
  plot(all_delta, type = "l")
  plot(x_train_scale,y_scale)
  points(x_train_scale,colMeans(all_y_hat[1001:2000,]),pch= 20, col = "blue")
  # ============================================

  curr <- 0
  tree_lengths <- numeric()

  for(i in 1:length(all_trees)){

      curr <- curr + 1

      for(j in 1:n_tree){
          tree_lengths[curr] <- length(all_trees[[i]][[j]])
      }

  }

  # tree_lengths |> table()
  #
  # gamma_tree <- matrix(0,nrow = 200, ncol = n_tree)
  #
  # for(i in 1:200){
  #
  #     for(k in 1:10){
  #
  #       gamma_tree[i,k] <- all_trees[[i]][[k]]$node0$betas_vec[1,1]
  #     }
  # }


}




