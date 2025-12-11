setwd("C:\\Users\\thchen7\\Downloads\\coding")

#
df3 <- readRDS("Dnase-seq_data.rds") 

set.seed(11111)
training_count = 750 #10% of peak counts (7464)
all_count = 7464 #all peak counts
peak_length = 100
trimmed_length = 1 #trimmed features length

SampleCount = training_count 
SampleLength = peak_length
sampledTrainingDataIndices = sample(0:all_count-1, SampleCount)
sampledTrainingRowStartIndices = (sampledTrainingDataIndices*(SampleLength-trimmed_length)+1)
sampledTrainingRowMask = rep(FALSE, nrow(df3))
for(i in 1:length(sampledTrainingRowStartIndices)){
  sampledTrainingRowMask[sampledTrainingRowStartIndices[i]:(sampledTrainingRowStartIndices[i]+SampleLength-trimmed_length-1)] = TRUE
}

SampledDF = df3[sampledTrainingRowMask,] #training_count*(peak_length-trimmed_length)
testingDF = df3[!sampledTrainingRowMask,] 

library(depmixS4)
library(dplyr)

# EM 
control_em <- function(maxit, tol, crit){
  em.control(maxit = maxit, tol = tol, crit = crit)
}

fit_and_evaluate <- function(formulas, data, nstates, family_list, ntimes,
                             transition_formula, replicates, em_control) {
  
  best_model <- NULL
  best_bic <- .Machine$double.xmax
  loglik_list <- c()
  aic_list <- c()
  bic_list <- c()
  
  for (i in 1:replicates) {
    model <- depmix(response = formulas, 
                    data = data, 
                    nstates = nstates, 
                    family = family_list, 
                    ntimes = ntimes, 
                    transition = transition_formula)
    
    fitted_model <- tryCatch({
      fit(model, emcontrol = em_control, verbose = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(fitted_model)) {
      current_bic <- BIC(fitted_model)
      current_aic <- AIC(fitted_model)
      current_ll <- logLik(fitted_model)
      
      bic_list <- c(bic_list, current_bic)
      aic_list <- c(aic_list, current_aic)
      loglik_list <- c(loglik_list, current_ll)
      
      if (current_bic < best_bic) {
        best_model <- fitted_model
        best_bic <- current_bic
      }
    }
  }
  
  return(list(
    best_model = best_model,
    best_bic = best_bic,
    best_aic = ifelse(length(aic_list) > 0, min(aic_list), NA),
    best_logLik = ifelse(length(loglik_list) > 0, max(loglik_list), NA)
  ))
}

# grid_search
grid_search_depmix <- function(param_grid, data, family_list, ntimes, 
                               replicates = 3) {
  results <- list()
  summary_df <- data.frame()
  
  for (i in 1:nrow(param_grid)) {
    row <- param_grid[i, ]
    
    message(sprintf("Trying config: R=%s, T=%s, nstates = %d, maxit = %d, tol = %.1e, crit = %s",
                    row$formulas, row$transition_formula, row$nstates, row$maxit, row$tol, row$crit))
    
    em_ctrl <- control_em(row$maxit, row$tol, as.character(row$crit))
    #print(flatten(row$formulas))
    #print(rep(list(gaussian()), length(flatten(row$formulas))))
    #print(row$transition_formula[[1]])
    res <- fit_and_evaluate(flatten(row$formulas), data, row$nstates, rep(list(gaussian()),length(flatten(row$formulas))), ntimes,
                            as.formula(row$transition_formula[[1]]), replicates, em_ctrl)

    summary_df <- rbind(summary_df, data.frame(
      response = as.character(row$formulas),
      transitions = as.character(row$transition_formula),
      nstates = row$nstates,
      maxit = row$maxit,
      tol = row$tol,
      crit = row$crit,
      BIC = res$best_bic,
      AIC = res$best_aic,
      logLik = res$best_logLik
    ))
    
    results[[i]] <- list(
      config = row,
      model = res$best_model
    )
  }
  
  return(list(
    summary = summary_df,
    models = results
  ))
}


# HMM
formulas_1 = list(
  MGW~1, HelT~1, ProT~1, Roll~1,
  EP~1
)
#DSDM
formulas_2 = list(
  MGW~1+MGWlag1+MGWlag2+MGWlag3,
  HelT~1+HelTlag1,
  ProT~1+ProTlag1,
  Roll~1,
  EP~1+EPlag1
)

transition_formula = ~1

param_grid <- expand.grid(
  formulas = list(formulas_1,formulas_2),
  transition_formula = c(transition_formula),
  nstates = c(4,5,6,7,8,9),         
  maxit = c(500),         
  tol = c(1e-08),        
  crit = c("relative")  
)

family_list = rep(list(gaussian()), 5)
ntimes = rep(99,training_count)
grid_data = SampledDF
replicates_number = 10

# 执行网格搜索
start_time <- proc.time()
grid_results <- grid_search_depmix(param_grid = param_grid,
                                   data = grid_data,
                                   family_list = family_list,
                                   ntimes = ntimes,
                                   replicates = replicates_number)
end_time <- proc.time()
print(paste("耗时（时）:", (end_time - start_time)["elapsed"]/3600))

# print summary
print(grid_results$summary)
summary_data <- as.data.frame(grid_results$summary)
# found BIC min
best_row <- grid_results$summary[which.min(grid_results$summary$BIC), ]
print(best_row)
# best model
best_model <- grid_results$models[[which.min(grid_results$summary$BIC)]]$model
best_model



all_train_performanceDF <- data.frame(
      Method = character(),
      Shape = character(),
      Model = character(), 
      ME = numeric(),
      RMSE = numeric(),
      MAE = numeric(),
      MPE = numeric(),
      MAPE = numeric(),
      stringsAsFactors = FALSE
    )
all_test_performanceDF <- data.frame(
    Method = character(),
    Shape = character(),
    Model = character(), 
    ME = numeric(),
    RMSE = numeric(),
    MAE = numeric(),
    MPE = numeric(),
    MAPE = numeric(),
    stringsAsFactors = FALSE
  )

for(model_idx in 1:length(grid_results$models)) {
  # Evaluate on testing data
  cat("\n\nEvaluating on testing data...")
  
  test_predicted_values <- list()
  for(param in shape_params) {
    test_predicted_values[[param]] <- rep(0, nrow(testingDF))
  }
  
  # Process testing data in chunks
  i <- 1
  while(i <= nrow(testingDF)) {
    if(i %% 1000 == 0) {
      cat(sprintf("\nProcessing testing row %d of %d", i, nrow(testingDF)))
    }
    
    chunk_end <- min(i + SampleLength - 4 - 1, nrow(testingDF))
    onePeak <- testingDF[i:chunk_end,]
    
    # Create test model with same structure as training
    model <- depmix(response = grid_results$models[[model_idx]]$config$formulas, 
                    data = onePeak, 
                    nstates = grid_results$models[[model_idx]]$config$nstates, 
                    family = rep(list(gaussian()), length(grid_results$models[[model_idx]]$config$formulas)), 
                    ntimes = rep(99,nrow(onePeak)), 
                    transition = as.formula(grid_results$models[[model_idx]]$config$transition_formula[[1]]))

    # Set parameters from trained model
    modTest <- setpars(modTest, getpars(current_model))
    viterbiResults <- viterbi(modTest)
    
    # Generate predictions
    for(j in 1:(chunk_end-i+1)) {
      currentState <- viterbiResults$state[j]
      currentResponseModels <- current_model@response[[currentState]]
      
      for(param_idx in 1:length(shape_params)) {
        param <- shape_params[param_idx]
        coefficients <- currentResponseModels[[param_idx]]@parameters$coefficients
        
        # Get lagged values based on formula
        if(param == "MGW") {
          valueVec <- c(1, onePeak[j,c('MGWlag1','MGWlag2','MGWlag3')])
        } else if(param %in% c("HelT","ProT","Roll","EP")) {
          valueVec <- c(1, onePeak[j,paste0(param,'lag1')])
        } else {
          valueVec <- 1
        }
        
        test_predicted_values[[param]][i+j-1] <- sum(unlist(valueVec) * coefficients)
      }
    }
    
    i <- i + SampleLength - 4
  }
  
  # Calculate and save testing performance metrics
  cat("\n\nCalculating testing performance metrics...")
  test_performanceDF <- data.frame(
    Method = character(),
    Shape = character(),
    Model = character(), 
    ME = numeric(),
    RMSE = numeric(),
    MAE = numeric(),
    MPE = numeric(),
    MAPE = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(param in shape_params) {
    metrics <- accuracy(test_predicted_values[[param]], testingDF[,param])
    test_performanceDF <- rbind(test_performanceDF,
                               c(grid_results$summary[model_idx,]$response, param, sprintf("Model_%d", model_idx), metrics))
  }
  
  colnames(test_performanceDF) <- c("Method", "Shape", "Model", "ME", "RMSE", "MAE", "MPE", "MAPE")
  
  write.csv(test_performanceDF,
            sprintf("model_%d_testing_performance.csv", model_idx),
            row.names = FALSE)
  
  cat(sprintf("\nModel %d evaluation complete - results saved to CSV files", model_idx))
}




MGWmodel = auto.arima(SampledDF[,'MGW'])
HelTmodel = auto.arima(SampledDF[,'HelT'])
ProTmodel = auto.arima(SampledDF[,'ProT'])
Rollmodel = auto.arima(SampledDF[,'Roll'])
EPmodel = auto.arima(SampledDF[,'EP'])

performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'MGW', 'Training', accuracy(MGWmodel)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'HelT', 'Training', accuracy(HelTmodel)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'ProT', 'Training', accuracy(ProTmodel)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'Roll', 'Training', accuracy(Rollmodel)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'EP', 'Training', accuracy(EPmodel)[1,1:5]))

performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'MGW', 'Testing', accuracy(testingDF[,'MGW'],(Arima(testingDF[,'MGW'], model=MGWmodel))$fitted)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'HelT', 'Testing', accuracy(testingDF[,'HelT'],(Arima(testingDF[,'HelT'], model=HelTmodel))$fitted)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'ProT', 'Testing', accuracy(testingDF[,'ProT'],(Arima(testingDF[,'ProT'], model=ProTmodel))$fitted)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'Roll', 'Testing', accuracy(testingDF[,'Roll'],(Arima(testingDF[,'Roll'], model=Rollmodel))$fitted)[1,1:5]))
performanceDF = rbind(performanceDF,c('ARIMA(auto)', 'EP', 'Testing', accuracy(testingDF[,'EP'],(Arima(testingDF[,'EP'], model=EPmodel))$fitted)[1,1:5]))


var.1 <- lineVar(SampledDF[,1:5], lag=1)
fittedValues = fitted(var.1)
fittedValues = rbind(colMeans(SampledDF[,1:5]),fittedValues)
performanceDF = rbind(performanceDF,c('VAR(1)', 'MGW', 'Training', accuracy(fittedValues[,'MGW'], SampledDF[,'MGW'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'HelT', 'Training', accuracy(fittedValues[,'HelT'], SampledDF[,'HelT'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'ProT', 'Training', accuracy(fittedValues[,'ProT'], SampledDF[,'ProT'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'Roll', 'Training', accuracy(fittedValues[,'Roll'], SampledDF[,'Roll'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'EP', 'Training',accuracy(fittedValues[,'EP'], SampledDF[,'EP'])))
fittedValues2 = matrix(0, nrow = nrow(testingDF), ncol = 5)
fittedValues2[1,] = colMeans(SampledDF[,1:5])
for(i in 1:(nrow(testingDF)-1)){
  fittedValues2[i+1,] =  predict(var.1, newdata=testingDF[i,1:5], n.ahead=1)
}
colnames(fittedValues2) = c('MGW','HelT','ProT','Roll','EP')
performanceDF = rbind(performanceDF,c('VAR(1)', 'MGW', 'Testing', accuracy(fittedValues2[,'MGW'], testingDF[,'MGW'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'HelT', 'Testing', accuracy(fittedValues2[,'HelT'], testingDF[,'HelT'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'ProT', 'Testing', accuracy(fittedValues2[,'ProT'], testingDF[,'ProT'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'Roll', 'Testing', accuracy(fittedValues2[,'Roll'], testingDF[,'Roll'])))
performanceDF = rbind(performanceDF,c('VAR(1)', 'EP', 'Testing',accuracy(fittedValues2[,'EP'], testingDF[,'EP'])))



var.2 <- lineVar(SampledDF[,1:5], lag=2)
fittedValues = fitted(var.2)
fittedValues = rbind(colMeans(SampledDF[,1:5]),fittedValues)
fittedValues = rbind(colMeans(SampledDF[,1:5]),fittedValues)
performanceDF = rbind(performanceDF,c('VAR(2)', 'MGW', 'Training', accuracy(fittedValues[,'MGW'], SampledDF[,'MGW'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'HelT', 'Training', accuracy(fittedValues[,'HelT'], SampledDF[,'HelT'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'ProT', 'Training', accuracy(fittedValues[,'ProT'], SampledDF[,'ProT'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'Roll', 'Training', accuracy(fittedValues[,'Roll'], SampledDF[,'Roll'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'EP', 'Training',accuracy(fittedValues[,'EP'], SampledDF[,'EP'])))
fittedValues2 = matrix(0, nrow = nrow(testingDF), ncol = 5)
fittedValues2[1,] = colMeans(SampledDF[,1:5])
fittedValues2[2,] = colMeans(SampledDF[,1:5])
for(i in 1:(nrow(testingDF)-2)){
  fittedValues2[i+2,] = predict(var.2, newdata=testingDF[i:(i+1),1:5], n.ahead=1)
}
colnames(fittedValues2) = c('MGW','HelT','ProT','Roll','EP')
performanceDF = rbind(performanceDF,c('VAR(2)', 'MGW', 'Testing', accuracy(fittedValues2[,'MGW'], testingDF[,'MGW'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'HelT', 'Testing', accuracy(fittedValues2[,'HelT'], testingDF[,'HelT'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'ProT', 'Testing', accuracy(fittedValues2[,'ProT'], testingDF[,'ProT'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'Roll', 'Testing', accuracy(fittedValues2[,'Roll'], testingDF[,'Roll'])))
performanceDF = rbind(performanceDF,c('VAR(2)', 'EP', 'Testing',accuracy(fittedValues2[,'EP'], testingDF[,'EP'])))



var.3 <- lineVar(SampledDF[,1:5], lag=3)
fittedValues = fitted(var.3)
fittedValues = rbind(colMeans(SampledDF[,1:5]),fittedValues)
fittedValues = rbind(colMeans(SampledDF[,1:5]),fittedValues)
fittedValues = rbind(colMeans(SampledDF[,1:5]),fittedValues)
performanceDF = rbind(performanceDF,c('VAR(3)', 'MGW', 'Training', accuracy(fittedValues[,'MGW'], SampledDF[,'MGW'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'HelT', 'Training', accuracy(fittedValues[,'HelT'], SampledDF[,'HelT'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'ProT', 'Training', accuracy(fittedValues[,'ProT'], SampledDF[,'ProT'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'Roll', 'Training', accuracy(fittedValues[,'Roll'], SampledDF[,'Roll'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'EP', 'Training',accuracy(fittedValues[,'EP'], SampledDF[,'EP'])))
fittedValues2 = matrix(0, nrow = nrow(testingDF), ncol = 5)
fittedValues2[1,] = colMeans(SampledDF[,1:5])
fittedValues2[2,] = colMeans(SampledDF[,1:5])
fittedValues2[3,] = colMeans(SampledDF[,1:5])
for(i in 1:(nrow(testingDF)-3)){
  fittedValues2[i+3,] = predict(var.3, newdata=testingDF[i:(i+2),1:5], n.ahead=1)
}
colnames(fittedValues2) = c('MGW','HelT','ProT','Roll','EP')
performanceDF = rbind(performanceDF,c('VAR(3)', 'MGW', 'Testing', accuracy(fittedValues2[,'MGW'], testingDF[,'MGW'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'HelT', 'Testing', accuracy(fittedValues2[,'HelT'], testingDF[,'HelT'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'ProT', 'Testing', accuracy(fittedValues2[,'ProT'], testingDF[,'ProT'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'Roll', 'Testing', accuracy(fittedValues2[,'Roll'], testingDF[,'Roll'])))
performanceDF = rbind(performanceDF,c('VAR(3)', 'EP', 'Testing',accuracy(fittedValues2[,'EP'], testingDF[,'EP'])))


