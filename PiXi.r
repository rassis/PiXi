library(dplyr)
library(readr)
library(mvtnorm) # For multivariate Gaussian distribution
library(moments) # For generating moments


GenerateTrainingData <- function(m, Nk, background_filename, training_prefix, logTransformExpression, logThetaMin, logThetaMax, logAlphaMin, logAlphaMax, logSigmaSqMin, logSigmaSqMax) {
	
	scenarios <- c("Conserved", "Diverged")		#Two scenarios here: "Conserved" (when no expression divergence occurs) and "Diverged" (when expression divergence occurs)
	tissues <- seq(1, m)
  training_data <- matrix(0, nrow = 2*Nk, ncol = 1 + 6*m)         #m tissues with 2 expression values and 4 parameters each and class label	
  
  for(sindex in 1:length(scenarios)) {
    scenario <- scenarios[sindex]
    for(i in 1:Nk) {
      alpha <- 10^runif(m, min = logAlphaMin, max = logAlphaMax)
      sigmaSq <- 10^runif(m, min = logSigmaSqMin, max = logSigmaSqMax)

      theta1 <- c()
      theta2 <- c()
      
      if(scenario == "Conserved") {
        theta1 <- runif(m, min = logThetaMin, max = logThetaMax)
        theta2 <- theta1
      }
      else if(scenario == "Diverged") {
        theta1 <- runif(m, min = logThetaMin, max = logThetaMax)
        theta2 <- runif(m, min = logThetaMin, max = logThetaMax)
      }

      expression_vec <- c()
			mu <- rep(0, 2)
			CovMat <- matrix(0, nrow = 2, ncol = 2)
      			
			for(j in 1:m) {
				mu[1] = (1 - exp(-alpha[j])) * theta2[j] + exp(-alpha[j]) * theta1[j]
				mu[2] = theta1[j]
				CovMat[1, 1] = sigmaSq[j] / (2 * alpha[j])
				CovMat[2, 2] = CovMat[1, 1]
				CovMat[1, 2] = exp(-2 * alpha[j]) * sigmaSq[j] / (2 * alpha[j])
				CovMat[2, 1] = CovMat[1, 2]
				expression_vec <- c(expression_vec, rmvnorm(1, mean = mu, sigma = CovMat))
			}
      
			rowIndex <- (sindex - 1)*Nk + i 
			training_data[rowIndex, 1] = sindex
      			
			for(j in 1:(2*m)) {
  			training_data[rowIndex, 1 + j] = expression_vec[j]
			}
      			
			for(j in 1:m) {
    		training_data[rowIndex, 1 + 2*m + j] = theta1[j]
      }
      
			for(j in 1:m) {
			  training_data[rowIndex, 1 + 3*m + j] = theta2[j]
			}
            		
			for(j in 1:m) {
      	training_data[rowIndex, 1 + 4*m + j] = log10(alpha[j])
      }

	    for(j in 1:m) {
        	training_data[rowIndex, 1 + 5*m + j] = log10(sigmaSq[j])
      }
    }
  }
  
	column_labels <- c("Class")

	for(j in 1:m) {
		column_labels <- c(column_labels, paste("e1", j, sep = ""), paste("e2", j, sep = ""))
	}

	for(j in 1:m) {
		column_labels <- c(column_labels, paste("Theta1", j, sep = ""))
	}

	for(j in 1:m) {
	  column_labels <- c(column_labels, paste("Theta2", j, sep = ""))			
	}

	for(j in 1:m) {
  	column_labels <- c(column_labels, paste("Alpha", j, sep = ""))
	}

	for(j in 1:m) {
  	column_labels <- c(column_labels, paste("SigmaSq", j, sep = ""))
	}

	colnames(training_data) <- column_labels

	write.table(training_data, file = paste(training_prefix, ".data", sep = ""), row.names = FALSE)

	if(logTransformExpression == FALSE) {
	  GenerateFeatures(m, background_filename, paste(training_prefix, ".data", sep = ""), training_prefix, FALSE)
	}
	
	if(logTransformExpression == TRUE) {
	  GenerateFeatures(m, background_filename, paste(training_prefix, ".data", sep = ""), training_prefix, TRUE)
	}
	
	training_data <- as_tibble( as.matrix(read.table(paste(training_prefix, ".data", sep = ""), header = TRUE)) ) %>% 
	  select(Class) %>% as.matrix()
	response <- matrix(0, nrow = nrow(training_data), ncol = 2)			#ncol = 2 because we have two classes
	
	for(i in 1:nrow(training_data)) {
	  response[i, training_data[i,1]] = 1
	}
	
	colnames(response) <- c("Conserved", "Diverged")
	
	write.table(response, file = paste(training_prefix, ".classes", sep = ""), row.names = FALSE)
	
	training_data <- as_tibble(as.matrix(read.table(paste(training_prefix, ".data", sep = ""), header = TRUE)) ) %>% 
	  select(starts_with("Theta"), starts_with("Alpha"), starts_with("Sigma")) %>% as.matrix()
	
	write.table(training_data, file = paste(training_prefix, ".responses", sep = ""), row.names = FALSE)
	
	Y_means <- colMeans(training_data)
	Y_sds <- apply(training_data, 2, sd)
	std_params <- data.frame("Ymeans" = Y_means, "Ysds" = Y_sds)
	write.table(std_params, paste(training_prefix, ".Y_stdparams", sep = ""), row.names = FALSE)
	rm(std_params)
}

GenerateFeatures <- function(m, background_filename, sample_filename, feature_filename, logTransformExpression) {
  minexp = 1e-4
	errorexp = 1e-5

	if(logTransformExpression == TRUE){
	  single <- read.table(background_filename, header = TRUE)
	  features <- as_tibble( as.matrix(read.table(sample_filename, header = TRUE)) ) %>% 
	    select(starts_with("e1"), starts_with("e2"))
	}
  
	if(logTransformExpression == FALSE){
	  single <- read.table(background_filename, header = TRUE)
	  features <- as_tibble( as.matrix(read.table(sample_filename, header = TRUE)) ) %>% 
	    select(starts_with("e1"), starts_with("e2"))
	  
	  for(i in 1:nrow(single)) {
  	  for(j in 1:ncol(single)) {
    	  single[i,j] <- log10(single[i,j] + minexp + rnorm(1, 0, errorexp))
    	}
	  }
	  
	  for(i in 1:nrow(features)) {
	    for(j in 1:ncol(features)) {
	      features[i,j] <- log10(features[i,j] + minexp + rnorm(1, 0, errorexp))
	    }
	  }
	}
	
	eS1S2dist <- sqrt( rowSums((single[,1:m] - single[,(m+1):(2*m)])^2)  ) #Euclidean distance between background expression profiles
	maxS1S2dist <- max(abs(eS1S2dist))
	eS1S2cor <- c()
	for(i in 1:nrow(single)) {
  		eS1S2cor[i] = cor(as.numeric(single[i,1:m]), as.numeric(single[i,(m+1):(2*m)]), method = "pearson") 
	}
  	
	rm(single)
  
  
  eDistNvR <- sqrt( rowSums( ( select(features, starts_with("e1")) - select(features, starts_with("e2")) )^2 ) ) # Euclidean distance between N and R #remove rest      
  
  features <- mutate(features, DistNvR = eDistNvR)
  
  rm(eDistNvR)
  
	features <- features %>%
  mutate(PBS_N = (DistNvR/2)) %>%
  mutate(PBS_R = (DistNvR/2)) %>%
  rowwise() %>%
  mutate(RankDistNvR = mean(eS1S2dist < DistNvR)) %>%
  mutate(DistNvR_m1 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 1, central = FALSE, absolute = FALSE)) %>%
  mutate(DistNvR_m2 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 2, central = FALSE, absolute = FALSE)) %>%
  mutate(DistNvR_m3 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 3, central = FALSE, absolute = FALSE)) %>%
  mutate(DistNvR_m4 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 4, central = FALSE, absolute = FALSE)) %>%
  mutate(DistNvR_m5 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 5, central = FALSE, absolute = FALSE)) %>%
  mutate(DistNvR_m6 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 6, central = FALSE, absolute = FALSE)) %>%
  mutate(DistNvR_m7 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 7, central = FALSE, absolute = FALSE)) %>%
  mutate(DistNvR_m8 = moment((eS1S2dist - DistNvR) / maxS1S2dist, order = 8, central = FALSE, absolute = FALSE)) %>%
  ungroup()
  
  eCorNvR <- c()
  
  	# Pearson correlations
  for(i in 1:nrow(features)) {
  		eCorNvR[i] <- cor( as.numeric(select(features[i,], starts_with("e1"))), as.numeric(select(features[i,], starts_with("e2"))), method = "pearson" )
  } 
  
  features <- mutate(features, CorNvR = eCorNvR)
  
  rm(eCorNvR)
  
	features <- features %>%
  rowwise() %>%
	mutate(RankCorNvR = mean(eS1S2cor < CorNvR)) %>%
	mutate(CorNvR_m1 = moment(eS1S2cor - CorNvR, order = 1, central = FALSE, absolute = FALSE)) %>%
	mutate(CorNvR_m2 = moment(eS1S2cor - CorNvR, order = 2, central = FALSE, absolute = FALSE)) %>%
	mutate(CorNvR_m3 = moment(eS1S2cor - CorNvR, order = 3, central = FALSE, absolute = FALSE)) %>%
	mutate(CorNvR_m4 = moment(eS1S2cor - CorNvR, order = 4, central = FALSE, absolute = FALSE)) %>%
	mutate(CorNvR_m5 = moment(eS1S2cor - CorNvR, order = 5, central = FALSE, absolute = FALSE)) %>%
	mutate(CorNvR_m6 = moment(eS1S2cor - CorNvR, order = 6, central = FALSE, absolute = FALSE)) %>%
	mutate(CorNvR_m7 = moment(eS1S2cor - CorNvR, order = 7, central = FALSE, absolute = FALSE)) %>%
	mutate(CorNvR_m8 = moment(eS1S2cor - CorNvR, order = 8, central = FALSE, absolute = FALSE)) %>%
	ungroup() %>% as.matrix()

	write.table(features, file = paste(feature_filename, ".features", sep = ""), row.names = FALSE)

	X_means <- colMeans(features)
	X_sds <- apply(features, 2, sd)
	std_params <- data.frame("Xmeans" = X_means, "Xsds" = X_sds)
	write.table(std_params, paste(feature_filename, ".X_stdparams", sep = ""), row.names = FALSE)
	rm(std_params)
}

ClassifierCVnn <- function(num_layers, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix) {
  if(num_layers %in% c(0, 1, 2, 3)) {
  	library(keras)
  	library(tensorflow)
  
  	CV <- 5
  	lambdas <- 10^seq(log_lambda_min, log_lambda_max, length = num_lambda)
  	gammas <- seq(gamma_min, gamma_max, length = num_gamma)
  
  	X <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
  	Y <- as.matrix(read.table(paste(training_prefix, ".classes", sep = ""), header = TRUE))
  
  	# standardize the input for training
  	X_means <- colMeans(X)
  	X_sds <- apply(X, 2, sd)
    
  	for(j in 1:ncol(X)) {
    	X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  	}
      
  	# CV-fold cross validation
  	val_loss <- array(0, dim = c(length(gammas), length(lambdas)))
  
  	# Randomly choose balanced training/validation sets per fold
  	foldid_Conserved <- sample(rep(seq(CV), length = nrow(X)/5))
  	foldid_Diverged <- sample(rep(seq(CV), length = nrow(X)/5))
  	foldid <- c(foldid_Conserved, foldid_Diverged)
  	rm(foldid_Conserved)
  	rm(foldid_Diverged)
      
  	# Perform K-fold CV, where K = CV
  	for(curr_fold in 1:CV) {
  	  Xval <- X[foldid == curr_fold, ]
  		Xtrain <- X[foldid != curr_fold, ]
  		Yval <- Y[foldid == curr_fold, ]
  		Ytrain <- Y[foldid != curr_fold, ]
  
  		# standardize the input for train and val based on train
  		temp_means <- colMeans(Xtrain)
  		temp_sds <- apply(Xtrain, 2, sd)
      			
  		for(j in 1:ncol(Xtrain)) {
      	Xtrain[,j] = (Xtrain[,j] - temp_means[j]) / temp_sds[j]
      	Xval[,j] = (Xval[,j] - temp_means[j]) / temp_sds[j]
  		}
      
  		for(i in 1:length(gammas)) {
        for(j in 1:length(lambdas)) {
  				model <- keras_model_sequential()
  				if(num_layers == 0) {
  						model %>%
  						layer_dense(units = 2,
              						activation = 'softmax',
              						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
              						input_shape = c(ncol(Xtrain)))
  				}
  
  				else if(num_layers == 1) {
  					model %>%
  				  layer_dense(units = 256, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
            						input_shape = c(ncol(Xtrain))) %>%
  					layer_dense(units = 2,
            						activation = 'softmax',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))
  				}
          
  				else if(num_layers == 2) {
      			model %>%
      			layer_dense(units = 256, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
            						input_shape = c(ncol(Xtrain))) %>%
      			layer_dense(units = 128, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
      			layer_dense(units = 2,
            						activation = 'softmax',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))
  		    }
          			
  				else if(num_layers == 3) {
  					model %>%
    				layer_dense(units = 256, 
              						activation = 'relu',
              						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
              						input_shape = c(ncol(Xtrain))) %>%
  					layer_dense(units = 128, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
  					layer_dense(units = 64, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
  					layer_dense(units = 2,
            						activation = 'softmax',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))
  				}
          
          
  			  model %>% compile(loss = 'categorical_crossentropy',
                  					optimizer = optimizer_adam(),
                  					metrics = c('categorical_crossentropy'))
        
      		history <- model %>% 
          					 fit(Xtrain, Ytrain,
              				  epochs = num_epochs,
              				  batch_size = batchsize,
              					validation_data = list(Xval, Yval),
              					verbose = 0) # verbose = 0  ensures it is silent
        
        	val_loss[i,j] <- val_loss[i,j] + as.data.frame(history) %>%
          filter(data == "validation", metric == "categorical_crossentropy") %>%
          select(value) %>% min()
    	  }
		  }
    
			rm(Xval)
			rm(Xtrain)
			rm(Yval)
			rm(Ytrain)
  	}
  
  	val_loss <- val_loss / CV
  
		gamma_opt = gammas[ which(val_loss == min(val_loss), arr.ind = TRUE)[1] ]
		lambda_opt = lambdas[ which(val_loss == min(val_loss), arr.ind = TRUE)[2] ]
		cv_results <- matrix(0, 1, 3)
		cv_results[,1] = min(val_loss)
		cv_results[,2] = gamma_opt
		cv_results[,3] = lambda_opt
		colnames(cv_results) <- c("Loss", "Gamma", "Lambda")

		write.table(cv_results, file = paste(training_prefix, ".nn.", num_layers, ".classifier_cv", sep = ""), row.names = FALSE)
		
		model <- keras_model_sequential()
		
		if(num_layers == 0) {
		  model %>%
		    layer_dense(units = 2,
		                activation = 'softmax',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
		                input_shape = c(ncol(X)))
		}
		
		else if(num_layers == 1) {
		  model %>%
		    layer_dense(units = 256, 
		                activation = 'relu',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
		                input_shape = c(ncol(X))) %>%
		    layer_dense(units = 2,
		                activation = 'softmax',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
		}
		
		else if(num_layers == 2) {
		  model %>%
		    layer_dense(units = 256, 
		                activation = 'relu',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
		                input_shape = c(ncol(X))) %>%
		    layer_dense(units = 128, 
		                activation = 'relu',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
		    layer_dense(units = 2,
		                activation = 'softmax',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
		}
		
		else if(num_layers == 3) {
		  model %>%
		    layer_dense(units = 256, 
		                activation = 'relu',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
		                input_shape = c(ncol(X))) %>%
		    layer_dense(units = 128, 
		                activation = 'relu',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
		    layer_dense(units = 64, 
		                activation = 'relu',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
		    layer_dense(units = 2,
		                activation = 'softmax',
		                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
		}
		
		model %>% 
		  compile(loss = 'categorical_crossentropy',
		          optimizer = optimizer_adam(),
		          metrics = c('categorical_crossentropy', 'accuracy'))
		
		history <- model %>% 
		  fit(X, Y,
		      epochs = num_epochs,
		      batch_size = batchsize,
		      verbose = 0) # verbose = 0  ensures it is silent
		
		model %>% save_model_hdf5(paste(training_prefix, ".nn.", num_layers, ".classifier.hdf5", sep = ""))
  }
}

PredictorCVnn <- function(num_layers, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix) {
	if(num_layers %in% c(0, 1, 2, 3)) {
    library(keras)
    library(tensorflow)
    
   	CV <- 5
		lambdas <- 10^seq(log_lambda_min, log_lambda_max, length = num_lambda)
		gammas <- seq(gamma_min, gamma_max, length = num_gamma)

		X <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
		Y <- as.matrix(read.table(paste(training_prefix, ".responses", sep = ""), header = TRUE))
    
    # standardize the input for training
		X_means <- colMeans(X)
		X_sds <- apply(X, 2, sd)
		Y_means <- colMeans(Y)
		Y_sds <- apply(Y, 2, sd)
    		
		for(j in 1:ncol(X)) {
    	X[,j] = (X[,j] - X_means[j]) / X_sds[j]
    }
    	
		for(j in 1:ncol(Y)) {
		  Y[,j] = (Y[,j] - Y_means[j]) / Y_sds[j]
    }
    
  	# CV-fold cross validation
  	val_loss <- array(0, dim = c(length(gammas), length(lambdas)))
  
  	# Randomly choose balanced training/validation sets per fold
  	foldid_Conserved <- sample(rep(seq(CV), length = nrow(X)/5))	
  	foldid_Diverged <- sample(rep(seq(CV), length = nrow(X)/5))
  	foldid <- c(foldid_Conserved, foldid_Diverged)
  	rm(foldid_Conserved)
  	rm(foldid_Diverged)
  
  	# Perform K-fold CV, where K = CV
    		
		for(curr_fold in 1:CV) {
			Xval <- X[foldid == curr_fold, ]
			Xtrain <- X[foldid != curr_fold, ]
			Yval <- Y[foldid == curr_fold, ]
			Ytrain <- Y[foldid != curr_fold, ]

			# standardize the input for train and val based on train
			temp_Xmeans <- colMeans(Xtrain)
			temp_Xsds <- apply(Xtrain, 2, sd)
			temp_Ymeans <- colMeans(Ytrain)
			temp_Ysds <- apply(Ytrain, 2, sd)
      
			for(j in 1:ncol(Xtrain)) {
  			Xtrain[,j] = (Xtrain[,j] - temp_Xmeans[j]) / temp_Xsds[j]
  			Xval[,j] = (Xval[,j] - temp_Xmeans[j]) / temp_Xsds[j]
			}
      
			for(j in 1:ncol(Ytrain)) {
      	Ytrain[,j] = (Ytrain[,j] - temp_Ymeans[j]) / temp_Ysds[j]
        Yval[,j] = (Yval[,j] - temp_Ymeans[j]) / temp_Ysds[j]
      }
      
			for(i in 1:length(gammas)) {
  			for(j in 1:length(lambdas)) {
					model <- keras_model_sequential()
    			if(num_layers == 0) {
      			model %>%
        		layer_dense(units = ncol(Ytrain),
                    						activation = 'linear',
                    						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
                    						input_shape = c(ncol(Xtrain)))
    	    }
          			
					else if(num_layers == 1) {
            model %>%
						layer_dense(units = 256, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
            						input_shape = c(ncol(Xtrain))) %>%
						layer_dense(units = ncol(Ytrain),
            						activation = 'linear',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))
		    	}
          
					else if(num_layers == 2) {
						model %>%
						layer_dense(units = 256, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
            						input_shape = c(ncol(Xtrain))) %>%
						layer_dense(units = 128, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
						layer_dense(units = ncol(Ytrain),
            						activation = 'linear',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))
				  }
          
					else if(num_layers == 3) {
  					model %>%
						layer_dense(units = 256, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
            						input_shape = c(ncol(Xtrain))) %>%
						layer_dense(units = 128, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
						layer_dense(units = 64, 
            						activation = 'relu',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])) %>%
						layer_dense(units = ncol(Ytrain),
            						activation = 'linear',
            						kernel_regularizer = regularizer_l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]))
					}
          
        	model %>% 
          compile(loss = 'mean_squared_error',
                  				optimizer = optimizer_adam(),
                  				metrics = c('mean_squared_error'))
        	
        	history <- model %>% 
          fit(Xtrain, Ytrain,
	      			epochs = num_epochs,
    					batch_size = batchsize,
    					validation_data = list(Xval, Yval),
    					verbose = 0) # verbose = 0  ensures it is silent
        
       	 	val_loss[i,j] <- val_loss[i,j] + as.data.frame(history) %>%
          filter(data == "validation", metric == "mean_squared_error") %>%
          select(value) %>% min()
        }
      }
      
			rm(Xval)
			rm(Xtrain)
			rm(Yval)
			rm(Ytrain)
  	}
    
  	val_loss <- val_loss / CV
  
		gamma_opt = gammas[ which(val_loss == min(val_loss), arr.ind = TRUE)[1] ]
		lambda_opt = lambdas[ which(val_loss == min(val_loss), arr.ind = TRUE)[2] ]
		cv_results <- matrix(0, 1, 3)
		cv_results[,1] = min(val_loss)
		cv_results[,2] = gamma_opt
		cv_results[,3] = lambda_opt
		colnames(cv_results) <- c("Loss", "Gamma", "Lambda")

  	write.table(cv_results, file = paste(training_prefix, ".nn.", num_layers, ".predictor_cv", sep = ""), row.names = FALSE)
  	
  	model <- keras_model_sequential()
  	
  	if(num_layers == 0) {
  	  model %>%
  	    layer_dense(units = ncol(Y),
  	                activation = 'linear',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
  	                input_shape = c(ncol(X)))
  	}
  	
  	else if(num_layers == 1) {
  	  model %>%
  	    layer_dense(units = 256, 
  	                activation = 'relu',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
  	                input_shape = c(ncol(X))) %>%
  	    layer_dense(units = ncol(Y),
  	                activation = 'linear',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
  	}
  	
  	else if(num_layers == 2) {
  	  model %>%
  	    layer_dense(units = 256, 
  	                activation = 'relu',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
  	                input_shape = c(ncol(X))) %>%
  	    layer_dense(units = 128, 
  	                activation = 'relu',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
  	    layer_dense(units = ncol(Y),
  	                activation = 'linear',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
  	}
  	
  	else if(num_layers == 3) {
  	  model %>%
  	    layer_dense(units = 256, 
  	                activation = 'relu',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt),
  	                input_shape = c(ncol(X))) %>%
  	    layer_dense(units = 128, 
  	                activation = 'relu',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
  	    layer_dense(units = 64, 
  	                activation = 'relu',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt)) %>%
  	    layer_dense(units = ncol(Y),
  	                activation = 'linear',
  	                kernel_regularizer = regularizer_l1_l2(l1 = gamma_opt * lambda_opt, l2 = (1 - gamma_opt) * lambda_opt))
  	}
  	
  	model %>% 
  	  compile(loss = 'mean_squared_error',
  	          optimizer = optimizer_adam(),
  	          metrics = 'mean_squared_error')
  	
  	history <- model %>% 
  	  fit(X, Y,
  	      epochs = num_epochs,
  	      batch_size = batchsize,
  	      verbose = 0) # verbose = 0  ensures it is silent
  	
  	model %>% save_model_hdf5(paste(training_prefix, ".nn.", num_layers, ".predictor.hdf5", sep = ""))
  }
}

PiXiClassifyNN <- function(training_prefix, testing_prefix, num_layers) {
	if(num_layers %in% c(0, 1, 2, 3)) {
		library(keras)
		library(tensorflow)

		X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))

		# standardize the input for testing
		std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
		X_means <- c(std_params$Xmeans)
		X_sds <- c(std_params$Xsds)
		
		for(j in 1:ncol(X)) {
  		X[,j] = (X[,j] - X_means[j]) / X_sds[j]
		}
    		
		rm(std_params)
    
		model <- load_model_hdf5(paste(training_prefix, ".nn.", num_layers, ".classifier.hdf5", sep = ""))

		Yest_num <- as.array(model %>% predict(X) %>% k_argmax())
		
		Yest <- data.frame("Class" = ifelse(Yest_num == 0, "Conserved", "Diverged"))

		write.table(Yest, paste(testing_prefix, ".nn.", num_layers, ".classifications", sep = ""), row.names = FALSE)

		probs_est <- model %>% predict(X)
		colnames(probs_est) <- c("Conserved", "Diverged")

		write.table(probs_est, paste(testing_prefix, ".nn.", num_layers, ".probabilities", sep = ""), row.names = FALSE)
  }
}

PiXiClassifySVM <- function(training_prefix, testing_prefix) {
	library(liquidSVM) # For SVM

	Xtrain <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
	Ytrain <- read.table(paste(training_prefix, ".classes", sep = ""), header = TRUE) %>%
  transmute(class = as.factor(ifelse(Conserved == 1, "Conserved", "Diverged"))) 
	
	X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))

	# standardize the input for training and testing
	std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
	X_means <- c(std_params$Xmeans)
	X_sds <- c(std_params$Xsds)
  
	for(j in 1:ncol(X)) {
  	Xtrain[,j] = (Xtrain[,j] - X_means[j]) / X_sds[j]
  	X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  	
	rm(std_params)
  
	svm_trained <- mcSVM(x = Xtrain, y = Ytrain$class, mc_type = "OvA_hinge", folds = 5, scale = FALSE, do.select = TRUE)
	preds <- c(predict(svm_trained, X))

	Yest <- data.frame("Class" = ifelse(preds <= 0.0, "Conserved", "Diverged"))

	write.table(Yest, paste(testing_prefix, ".svm.classifications", sep = ""), row.names = FALSE)

	pred_probs <- (preds + 1) / 2
	probs_est <- data.frame("Conserved" = 1 - pred_probs, "Diverged" = pred_probs)

	write.table(probs_est, paste(testing_prefix, ".svm.probabilities", sep = ""), row.names = FALSE)
}

PiXiClassifyRF <- function(training_prefix, testing_prefix) {
	library(ranger) # For RF
  
  Xtrain <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
	Ytrain <- read.table(paste(training_prefix, ".classes", sep = ""), header = TRUE) %>%
	transmute(class = as.factor(ifelse(Conserved == 1, "Conserved", "Diverged"))) 
  
	X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))
  
  # standardize the input for training and testing
	std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
	X_means <- c(std_params$Xmeans)
	X_sds <- c(std_params$Xsds)
  
	for(j in 1:ncol(X)) {
  	Xtrain[,j] = (Xtrain[,j] - X_means[j]) / X_sds[j]
  	X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  
	rm(std_params)
  
	tempData <- as.data.frame( cbind(Xtrain, Ytrain) )
	rf_trained <- ranger(class ~., data = tempData, probability = TRUE)
	preds <- predict(rf_trained, X)

	probs_est <- as.data.frame(preds$predictions)
	Yest <- data.frame("Class" = ifelse(probs_est$Conserved > 0.5, "Conserved", "Diverged"))

	write.table(Yest, paste(testing_prefix, ".rf.classifications", sep = ""), row.names = FALSE)
	write.table(probs_est, paste(testing_prefix, ".rf.probabilities", sep = ""), row.names = FALSE)
}

PiXiPredictNN <- function(training_prefix, testing_prefix, num_layers) {
  if(num_layers %in% c(0, 1, 2, 3)) {
    library(keras)
    library(tensorflow)
    
    X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))
    
    # standardize the input for testing
    std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
    X_means <- c(std_params$Xmeans)
    X_sds <- c(std_params$Xsds)
    
    for(j in 1:ncol(X)) {
      X[,j] = (X[,j] - X_means[j]) / X_sds[j]
    }
    
    rm(std_params)
    
    model <- load_model_hdf5(paste(training_prefix, ".nn.", num_layers, ".predictor.hdf5", sep = ""))
    
    std_params <- read.table(paste(training_prefix, ".Y_stdparams", sep = ""), header = TRUE)
    Y_means <- c(std_params$Ymeans)
    Y_sds <- c(std_params$Ysds)
    rm(std_params)
    
    Yest_std <- model %>% predict(X)
    Yest <- Yest_std
    
    # un-standardize the responses
    for(j in 1:ncol(Yest)) {
      Yest[, j] = Yest_std[,j] * Y_sds[j] + Y_means[j]
    }
    
    Y <- as.matrix(read.table(paste(training_prefix, ".responses", sep = ""), header = TRUE))
    colnames(Yest) <- colnames(Y)
    rm(Y)
    
    write.table(Yest, paste(testing_prefix, ".nn.", num_layers, ".predictions", sep = ""), row.names = FALSE)
  }
}

PiXiPredictSVM <- function(training_prefix, testing_prefix) {
  library(liquidSVM) # For SVM
  Xtrain <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
  Ytrain <- as.matrix(read.table(paste(training_prefix, ".responses", sep = ""), header = TRUE)) 
  X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))
  
  # standardize the input for training and testing
  std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
  X_means <- c(std_params$Xmeans)
  X_sds <- c(std_params$Xsds)
  
  for(j in 1:ncol(X)) {
    Xtrain[,j] = (Xtrain[,j] - X_means[j]) / X_sds[j]
    X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  
  rm(std_params)
  
  std_params <- read.table(paste(training_prefix, ".Y_stdparams", sep = ""), header = TRUE)
  Y_means <- c(std_params$Ymeans)
  Y_sds <- c(std_params$Ysds)
  
  for(j in 1:ncol(Ytrain)) {
    Ytrain[,j] = (Ytrain[,j] - Y_means[j]) / Y_sds[j]
  }
  
  rm(std_params)
  
  Yest <- matrix(0, nrow = nrow(X), ncol = ncol(Ytrain))
  
  for(i in 1:ncol(Ytrain)) {
    svm_trained <- lsSVM(x = Xtrain, y = Ytrain[,i], folds = 5, scale = FALSE, do.select = TRUE)
    Yest[,i] <- c(predict(svm_trained, X))
  }
  
  # un-standardize the responses
  for(j in 1:ncol(Yest)) {
    Yest[, j] = Yest[,j] * Y_sds[j] + Y_means[j]
  }
  
  colnames(Yest) <- colnames(Ytrain)
  rm(Ytrain)
  
  write.table(Yest, paste(testing_prefix, ".svm.predictions", sep = ""), row.names = FALSE)
}

PiXiPredictRF <- function(training_prefix, testing_prefix) {
  
  library(ranger)
  Xtrain <- as.matrix(read.table(paste(training_prefix, ".features", sep = ""), header = TRUE))
  Ytrain <- as.matrix(read.table(paste(training_prefix, ".responses", sep = ""), header = TRUE)) 
  X <- as.matrix(read.table(paste(testing_prefix, ".features", sep = ""), header = TRUE))
  
  # standardize the input for training and testing
  std_params <- read.table(paste(training_prefix, ".X_stdparams", sep = ""), header = TRUE)
  X_means <- c(std_params$Xmeans)
  X_sds <- c(std_params$Xsds)
  
  for(j in 1:ncol(X)) {
    Xtrain[,j] = (Xtrain[,j] - X_means[j]) / X_sds[j]
    X[,j] = (X[,j] - X_means[j]) / X_sds[j]
  }
  
  rm(std_params)
  
  std_params <- read.table(paste(training_prefix, ".Y_stdparams", sep = ""), header = TRUE)
  Y_means <- c(std_params$Ymeans)
  Y_sds <- c(std_params$Ysds)
  
  for(j in 1:ncol(Ytrain)) {
    Ytrain[,j] = (Ytrain[,j] - Y_means[j]) / Y_sds[j]
  }
  
  rm(std_params)
  
  Yest <- matrix(0, nrow = nrow(X), ncol = ncol(Ytrain))
  
  for(i in 1:ncol(Ytrain)) {
    tempData <- as.data.frame( cbind(Xtrain, Ytrain[,i]) )
    colnames(tempData) <- c(colnames(Xtrain), "regResp")
    rf_trained <- ranger(regResp ~., data = tempData)
    Yest[,i] <- c(predict(rf_trained, X)$predictions)
  }
  
  # un-standardize the responses
  for(j in 1:ncol(Yest)) {
    Yest[, j] = Yest[,j] * Y_sds[j] + Y_means[j]
  }
  
  colnames(Yest) <- colnames(Ytrain)
  rm(Ytrain)
  
  write.table(Yest, paste(testing_prefix, ".rf.predictions", sep = ""), row.names = FALSE)
}