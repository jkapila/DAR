#########################################################################
#
#         Implementation of Distribution Assertive Regression
#
#########################################################################


# Definign Mean Average Percentage Error
mape <- function(actual,predicted,asmean=FALSE){
  ape <- (abs(actual-predicted) * 100)/ actual
  ape <- if(asmean)
    mean(ape)
  else
    round(ape,5)
  ape
}

decileBinner <- function(data,target_var,splitname,breaks = 10){
  
  if(missing(splitname)){
    splitname <- "splits"
  }
  data <- data %>%  mutate(s = as.numeric(cut(data[,target_var],
                                              breaks = breaks,rigth = T)))
  name_ <- colnames(data)
  name_ [name_=="s"] <- splitname
  colnames(data) <- name_
  return(data)
}

# Normalising the data
normalise <- function(x){
  min_ <- min(x,na.rm = T)
  max_ <- max(x,na.rm = T)
  y <- (x-min_)/(max_-min_)
  return(y)
}

dummyVar <- function(data,name,keepVar=FALSE){
  data_ <- data.frame(data)
  colnames(data_) <- name
  for(t in unique(data_[,name])) {
    new_col <- paste(name,t,sep="_")
    data_[,new_col] <- ifelse(data_[,name]==t,1,0)
  }
  if(keepVar){
    data <- cbind(data,data_)
    return(data)
  }else{
    return(data_)
  }
}

# Segregating data
# Aproach 1: KNN
## this is to be impleted from c for calculating distance
## and learning knn classifications

## simple knn function
knn <- function(mat, k){
  cat('\nGot Data :',nrow(mat),'\n')
  if(!is.matrix(mat)){
    mat <- as.matrix(mat)
  }
  n <- nrow(mat)
  if (n <= k) stop(" kNN's k can not be more than nrow(data)-1! 
                   Reduce k and/or increase samples!  ")
  neigh <- matrix(0, nrow = n, ncol = k)
  # library(fields)
  ## This sholud be looped in chunks or implemented in C++
  dist.mat <- fields::rdist(mat, mat)
  # print(dist.mat)
  for(i in 1:n) {
    euc.dist <- dist.mat[i, ]
    # print(euc.dist)
    neigh[i, ] <- order(euc.dist)[2:(k + 1)]
  }
  
  return(neigh)
}

knn(iris[,1:4],3)

dafr <- function(formula,df,model,family,dec.front=c(2),dec.back=c(2),
                 knn.neighbours=5,...){
  
  dname <- paste(deparse(substitute(formula)))
  # definign missing values
  if (missing(model)){
    model <- glm
  }
  if (missing(family)){
    family <- "gaussian"
  }
  
  # breaking the df set for initial split
  # y_actual <- df[,formula$y]
  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y_actual <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
    # Z <- (rownames(data) %in% cut(y_actual,breaks = 10,right = TRUE))
  }
  else {
    mf <- model.frame(formula, data = df)
    y_actual <- model.response(mf)
    X <- model.matrix(formula, data = df)
    # Z <- (rownames(data) %in% cut(y_actual,breaks = 10,right = TRUE))
  }
  
  # making original base model
  mod_orig <- model(formula = formula ,data = df, family = family,...=...)
  y_orig <- predict(mod_orig,df)
  
  # geting mape curve
  results <- data.frame(actuals=y_actual,original=y_orig)
  results <- decileBinner(results,"actuals",splitname = "splits")
  rownames(results) <- rownames(df)
  cat('\nActual and Prediction by Single Model :\n')
  print(results)
  curve_ape <- summarize(group_by(results,splits),mape(actuals,original,asmean=TRUE))
  curve_ape <- data.frame(curve_ape)
  colnames(curve_ape) <- c("splits","mape")
  cat('\nSplits MAPE distributions:\n')
  print(curve_ape)
  plot(curve_ape,type="b",main =" Plot of Unsplitted Absolute Percentage Error",
       ylab = "Mean Absolute Percentage Error",xlab="Split Index")
  
  # looking at split distirbution
  hist(results[,"splits"],main = "Split Distribution",xlab = "Split Index",breaks = 10)
  cat('\n Early Failure Region : ',dec.front," Wear Out Failure Region :",dec.back,'\n')
  
  # Breaking the dataset by deciles and remodelling
  # vectorising the deciles
  if(NROW(dec.front)==1&!is.null(dec.front)){
    dec.front <- seq(1,dec.front)
    cat('\nFront Splits:',dec.front)
    # print(results[results$splits %in% dec.front,])
    front_idx <- row.names(results[results$splits %in% dec.front,])
    # data_front <- df[,]
    # print(knn(data_front,knn.neighbours))
    cat("\nFront Data has: ",
        NROW(front_idx)," rows")
  }else{
    dec.front <- c()
  }
  if(NROW(dec.back)==1&!is.null(dec.back)){
    dec.back <- seq(10,(11-dec.back))
    cat('\nBack Splits:',dec.back)
    # data_back <- df[rownames(results[results$splits %in% dec.back,]),]
    back_idx <- row.names(results[results$splits %in% dec.back,])
    cat("\nBack Data has: ",NROW(back_idx)," rows")
  }else{
    dec.back <- c()
  }
  
  # data_mid <- df[rownames(results[!results$splits %in% c(dec.back,dec.front),]),]
  mid_idx <- rownames(results[!results$splits %in% c(dec.back,dec.front),])
  cat("\nMid Data has: ",
      NROW(mid_idx)," rows\nData has: ",nrow(df)," rows")
  
  # generating splitted models and calculating mapes
  pred_dec <- c()
  models <- list()
  if(length(dec.front)>0){
    
    mod_front <- model(formula = formula ,data = df[front_idx,], family = family,...=...)
    models$Front_Model <- mod_front
    # cat(nrow(dff <- df[front_idx,]))
    models$knn.front <- knn(df[front_idx,],knn.neighbours)
    pred_front <- predict(mod_front,df[front_idx,])
    pred_dec <- pred_front
  }
  if(NROW(mid_idx)>0){
    mod_mid <- model(formula = formula ,data = df[mid_idx,], family = family,...=...)
    models$Mid_Model <- mod_mid
    models$knn.mid <- knn(df[mid_idx,],knn.neighbours)
    pred_mid <- predict(mod_mid,df[mid_idx,])
    pred_dec <- c(pred_dec,pred_mid)
  }
  if(length(dec.back)>0){
    mod_back <- model(formula = formula ,data = df[back_idx,], family = family,...=...)
    models$Back_Model <- mod_back
    models$knn.back <- knn(df[back_idx,],knn.neighbours)
    pred_back <- predict(mod_back,df[back_idx,])
    pred_dec <- c(pred_dec,pred_back)
  }
  cat("\nDeciled Prediction has:",NROW(pred_dec)," value \n")
  # Replotting curve of mape
  results[,"dec_pred"] <- pred_dec
  print(results)
  curve_ape_dec <- summarize(group_by(results,splits),mape(actuals,dec_pred,asmean=TRUE))
  curve_ape_dec <- data.frame(curve_ape_dec)
  colnames(curve_ape_dec) <- c("split","mape_dec")
  curve_ape <- cbind(curve_ape,round(curve_ape_dec[["mape_dec"]],2))
  colnames(curve_ape)[3]  <- "mape_dec"
  print(curve_ape)
  plot(curve_ape_dec,type="b",main =" Plot of Splitted Absolute Percentage Error",
       ylab = "Mean Absolute Percentage Error",xlab="Split Index")
  
  
  # models <- list(exists(mod_front),exists(mod_mid),exists(mod_back))
  dafr <- list(formula = dname,models=models,results= results[,c(3,1,2,4)],
               mapes=curve_ape,split.freq=table(results[,"splits"]))
  return(dafr)
  
}

# Using available data mtcars and LifeCycleSavings
data("mtcars")
df2 <- mtcars[order(mtcars[,"mpg"]),]
data("LifeCycleSavings")
df <- LifeCycleSavings[order(LifeCycleSavings[,"sr"]),]

knn(df,1)

library(mlbench)
data("BostonHousing")
df <- BostonHousing[order(BostonHousing[,'medv']),]

# Using our t
mod_dfar <- dafr(medv ~ ., df = df #, dec.back = 3,knn.neighbours = 3
)
summary(mod_dfar)
mod_dfar$call
mod_dfar$models
mod_dfar$results
mod_dfar$mapes
mod_dfar$split.freq

mod_dfar <- dafr(mpg ~. , data = df2)
mod_dfar <- dafr(mpg ~. , data = df2,dec.front = 3)
mod_dfar <- dafr(mpg ~. , data = df2,dec.front = 3,dec.back = 3)
mod_dfar <- dafr(mpg ~. , data = df2,dec.front = 4,dec.back = 3)

mod_dfar <- dafr(sr ~ ., data = df,model = lm)
mod_dfar <- dafr(sr ~ ., data = df,dec.front = 3,model = lm)
mod_dfar <- dafr(sr ~ ., data = df,dec.front = 3,dec.back = NULL,model = lm)
mod_dfar <- dafr(sr ~ ., data = df,dec.front = 3,dec.back = 3,model = lm)

mod_dfar <- dafr(medv ~. , data = df)
mod_dfar <- dafr(medv ~. , data = df,dec.front = 3)
mod_dfar <- dafr(medv ~. , data = df,dec.front = 3,dec.back = 3)
mod_dfar <- dafr(medv  ~. , df = df, dec.front = 2, dec.back = 3)


# Making Prediction on data
predict.dafr <- function(dafr,data){
  if (!is.null(r <- get0(myVarName, envir = myEnvir))) {
    ## ... deal with r ...
  }
  
}
