# PiXi
PiXi is an R implementation of Piya, DeGiorgio, and Assis' (2022) three machine learning algorithms for predicting expression divergence between single-copy orthologs and their expression optima in two species.

-----------
Citing PiXi
-----------

Thank you for using PiXi.

If you use this R package, then please cite the bioRxiv preprint:
    AA Piya, M DeGiorgio, R Assis (2022) Predicting expression divergence and its evolutionary parameters between single-copy genes in two species. bioRxiv doi:           https://doi.org/10.1101/2022.07.13.499803.
	
---------------	
Getting started
---------------
Before you are able to use the PiXI R package, you will need to have the dplyr, readr, mvtnorm, and moments libraries installed in your R environment. These libraries can be installed with the following commands in R:

  install.packages("dplyr")
  install.packages("readr")
  install.packages("mvtnorm")
  install.packages("moments")
  
To use the PiXi neural network, you must also install the tensorflow, keras, and miniconda libraries with the following commands:
  install.packages("tensorflow")
  library(tensorflow)
  install_tensorflow()
  install.packages("keras")
  reticulate::install_miniconda()	#needed for some operating systems
  
To use the PiXi random forest, you must also install the ranger library with the following command:

  install.packages("ranger")

To use the PiXi support vector machine, you must also install the liquidSVM library with the following command: 

  install.packages("liquidSVM", repos="http://pnp.mathematik.uni-stuttgart.de/isa/steinwart/software/R")
  
The PiXi package comes with the script PiXi.r and an ExampleFiles directory containing example files to help you get started.
  
The PiXi package can be loaded in R with the command:

  source("PiXi.r")

Note that this command assumes that you are in the directory where the PiXi.r script is located.

To run PiXi, you will need to provide an input file containing expression data for a set of genes. The format of this file is described below, and an example file called EmpiricalData is provided in the ExampleFiles directory.

--------------------
Format of input file
--------------------
A space-delimited file with N+1 rows and 2m columns, where N is the number of genes and m is the number of conditions (e.g., tissues, ages, etc.) for expression data. The first row is a header.

Each row (rows 2 to N+1) is a gene, and each column is an absolute log10(x+1) expression value in a particular condition for that gene. The first m columns are the gene expression values for the m conditions in the first species, and the next m columns are the gene expression values for the m conditions in the second species. Specifically, columns j and m+j represent the gene expression values at condition j (j = 1, 2, ..., m) in species 1 and species 2, respectively.

The header takes the form:

  "e11" "e12" ... "e1m" "e21" "e22" ... "e2m"

where e11 to e1m denote log10-transformed gene expression values in species 1 for conditions 1 to m, and e21 to e2m denote log10-transformed gene expression values in species 2 for conditions 1 to m.
  
The ExampleFiles directory contains a file called EmpiricalData, which illustates this format for N=100 genes in m=6 tissues.

-------------------------------------------
Generating training data from an OU process
-------------------------------------------
The training data can be generated with the command:

  GenerateTrainingData(m, Nk, training_prefix, logThetaMin, logThetaMax, logAlphaMin, logAlphaMax, logSigmaSqMin, logSigmaSqMax)

where m is the number of conditions, Nk is the number of training observations for each of the two classes, training_prefix is the prefix given to all files output by this function, and the rest of the parameters are used to define the ranges that the evolutionary parameters Theta (optimal expression), Alpha (strength of selection) and Sigma-Squared (strength of phenotypic drift) of the OU process are drawn from.

The GenerateTrainingData() function also outputs the p=2m features to the file training_prefix.features, a one-hot coded matrix of classifications for all training observations to the file training_prefix.classes, a matrix of predicted expression optima for all training observations to the file training_prefix.responses, raw simulated data to the file training_prefix.data, and the means and standard deviations for each of the p features and each of the 2m model parameters in training_prefix.X_stdparams and training_prefix.Y_stdparams, respectively.

-------------
Training PiXi
-------------
Due to their speed, training of the PiXi support vector machine and random forest are performed on-the-fly while being applied to test data. Thus, this step only needs to be performed if using the PiXi neural network. Otherwise, the user can skip this section. 

The PiXi neural network is trained on the training data outputted by the GenerateTrainingData() function, as described in "Generating training data from an OU process". However, because estimating its many hyperparameters is a time-consuming process, we implement hyperparameter tuning of the PiXi neural network predictor separately from final model training. We perform five-fold cross-validation to identify optimal hyperparameters for the PiXi neural network predictor with num_layer layers (num_layer in {0, 1, 2, 3}), regularization tuning parameter lambda, and elastic net tuning parameter gamma. Conditional on the optimal lambda and gamma hyperparameters, PiXi then fits a neural network predictor with num_layer hidden layers. The optimal number of hidden layers is chosen as the one with the smallest validation loss. 

The ClassifierCVnn() function is used to train the PiXi neural network to predict classes ("Conserved" or "Diverged"), and the PredictorCVnn() function is used to train the PiXi neural network to predict expression optima (Theta1 and Theta2). We consider log(lambda) evenly distributed within [log_lambda_min, log_lambda_max] for num_lambda values, and gamma evenly distributed within [gamma_min, gamma_max] with num_gamma values, assuming a batch size of batchsize observations per epoch and trained for num_epochs epochs with the commands:

  ClassifierCVnn(num_layers, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix)
  PredictorCVnn(num_layers, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix)

where num_layers is the number of layers (0, 1, 2, or 3) of the neural network, batchsize is the number of training observations used in each epoch, num_epochs is the number of training epochs, hyperparameter lambda is drawn from log10(lambda) in interval [log_lambda_min, log_lambda_max] for num_lambda evenly spaced points, hyperparameter gamma is drawn from interval [gamma_min, gamma_max] for num_gamma evenly spaced points, and training_prefix is the prefix to all files outputted by this function.

These functions output the set of optimal hyperparameters chosen through cross-validation to the files training_prefix.nn.num_layers.classifier_cv and training_prefix.nn.num_layers.predictor_cv, and the fitted neural network models in TensorFlow format to the files training_prefix.nn.num_layers.classifier.hdf5 and training_prefix.nn.num_layers.predictor.hdf5.

---------------------------
Performing test predictions
---------------------------
PiXi can predict classes ("Conserved" or "Diverged") and expression optima (Theta1 and Theta2) for each gene in a test dataset. 

To predict classes, the following commands can be used depending on the chosen algorithm:

  PiXiClassifyNN(training_prefix, testing_prefix, num_layers)	#neural network
  PiXiClassifyRF(training_prefix, testing_prefix)				#random forest
  PiXiClassifySVM(training_prefix, testing_prefix)				#support vector machine
  
where training_prefix and testing_prefix are the prefixes to training and testing feature files that were ouputted by the GenerateTrainingData() function, and num_layers is the number of layers (0, 1, 2, or 3) of the neural network.
  
The PiXiClassifyNN(), PiXiClassifyRF(), and PiXiClassifySVM() functions output predicted classes and probabilities for each gene in the test dataset to the respective files

  test_prefix.nn.num_layers.classifications
  test_prefix.rf.classifications
  test_prefix.svm.classifications

and

  test_prefix.nn.num_layers.probabilities
  test_prefix.rf.probabilities
  test_prefix.svm.probabilities

To predict expression optima, the following commands can be used depending on the chosen algorithm:

  PiXiPredictNN(training_prefix, testing_prefix, num_layers)	#neural network
  PiXiPredictRF(training_prefix, testing_prefix)				#random forest
  PiXiPredictSVM(training_prefix, testing_prefix)				#support vector machine
  
where training_prefix and testing_prefix are the prefixes to training and testing feature files that were ouputted by the GenerateTrainingData() function, and num_layers is the number of layers (0, 1, 2, or 3) of the neural network.
  
The PiXiPredictNN(), PiXiPredictRF(), and PiXiPredictSVM() functions output the 2m predicted expression optima for each gene in the test dataset to the respective files 

  test_prefix.nn.num_layers.predictions
  test_prefix.rf.predictions
  test_prefix.svm.predictions

----------------------------
Example application of PiXi
----------------------------
Within the R environment, set the working directory to the directory containing both the PiXi.r script and the subdirectory ExampleFiles containing the example files. 

Load the functions of the PiXi package by typing the command:
  
  source("PiXi.r")

Next, generate a training dataset of 100 observations in six conditions for each of the two classes, and store the training data with prefix Training, by typing the command:

  GenerateTrainingData(6, 100, "ExampleFiles/Training", 0, 5, 0, 3, -2, 3)

The above operation will output the files

  Training.data
  Training.features
  Training.classes
  Training.responses
  Training.X_stdparams
  Training.Y_stdparams

for which log10(Theta) is drawn between 0 and 5, log10(Alpha) is drawn between 0 and 3, and log10(SigmaSq) is drawn between -2 and 3.

To train the PiXi neural network with 2 hidden layers on this training datset using five-fold cross-validation with hyperparameters log(lambda) drawn from {-3, -2, -1, 0, 1, 2, 3} and gamma drawn from {0, 0.5, 1}, assuming a batch size of 50 observations per epoch and trained for 50 epochs, type the commands:

  ClassifierCVnn(2, 50, 50, -3, 3, 7, 0, 1, 3, "ExampleFiles/Training")
  PredictorCVnn(2, 50, 50, -3, 3, 7, 0, 1, 3, "ExampleFiles/Training")

The above operations will output the files 

  Training.nn.2.classifier_cv
  Training.nn.2.predictor_cv
  Training.nn.2.classifier.hdf5
  Training.nn.2.predictor.hdf5

Note that we perform training and hyperparameter tuning of the PiXi neural network in advance of application to test data, whereas the PiXi support vector machine and random forest classifiers and predictors are trained on-the-fly during application to test data.

Finally, to predict classes and expression optima for genes in an empirical dataset with the PiXi two-layer neural network, random forest, and support vector machine, type the commands:

  PiXiClassifyNN("ExampleFiles/Training", "ExampleFiles/EmpiricalData", 2)
  PiXiPredictNN("ExampleFiles/Training", "ExampleFiles/EmpiricalData", 2)

  PiXiClassifyRF("ExampleFiles/Training", "ExampleFiles/EmpiricalData")
  PiXiPredictRF("ExampleFiles/Training", "ExampleFiles/EmpiricalData")

  PiXiClassifySVM("ExampleFiles/Training", "ExampleFiles/EmpiricalData")
  PiXiPredictSVM("ExampleFiles/Training", "ExampleFiles/EmpiricalData")

The above operations will output the files

  EmpiricalData.nn.2.classifications
  EmpiricalData.nn.2.probabilities
  EmpiricalData.nn.2.predictions
  
  EmpiricalData.rf.classifications
  EmpiricalData.rf.probabilities
  EmpiricalData.rf.predictions

  EmpiricalData.svm.classifications
  EmpiricalData.svm.probabilities
  EmpiricalData.svm.predictions
