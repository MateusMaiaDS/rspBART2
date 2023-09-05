library(MOTRbart)
model <- motr_bart(x = cbind(x_train,x_train),y = y_train,sparse = FALSE)
plot(x_train$x,y_train)
points(x_train$x,model$y_hat |> colMeans(),pch = 20)
model$vars_betas
