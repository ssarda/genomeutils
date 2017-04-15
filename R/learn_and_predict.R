#' @title Learns and classifies using an SVM 
#' @description This function trains an SVM on the training set and tests its accuracy on a test set. It reports the fit, prediction results, area under curve (AUC accuracy) on test, decision values ( probability of belonging to predicted class) and original labels of test
#' @param Xtrain Predictor variables of training set
#' @param ytrain Response variables of training set
#' @param Xtest Predictor variables of test set
#' @param ytest Response variables of test set
#' @return A list with model fit, test results, AUC, probability per class, original labels
#' @export
test_svm = function(Xtrain, ytrain, Xtest, ytest)
{

	fit = svm(x=Xtrain, y=ytrain)
	results = predict(fit, Xtest, decision.values=TRUE)

	fx = attr(results, "decision.values")
	auc = -1
	if(length(fx) > 2 && length(table(ytest)) >= 2)
		auc = 0 + auc(ytest, fx)

	list(fit=fit, results=results, auc=auc, fx=fx, classes=ytest)
}



#' @title Learns and classifies using a Naive Bayes Classifier 
#' @description This function trains a NB classifier on the training set and tests its accuracy on a test set. It reports the fit, prediction results, area under curve (AUC accuracy) on test, probability of belonging to predicted class, and original labels of test
#' @param Xtrain Predictor variables of training set
#' @param ytrain Response variables of training set
#' @param Xtest Predictor variables of test set
#' @param ytest Response variables of test set
#' @return A list with model fit, test results, AUC, probability per class, original labels
#' @export
test_nb = function(Xtrain, ytrain, Xtest, ytest)
{

	fit = naiveBayes(x=Xtrain, y=ytrain)
	results = predict(fit, newdata=Xtest, type="raw")

	fx = results[, 1]
	auc = -1
	tryCatch({
		auc = 0 + auc(ytest, fx)
	}, error = function(e){ })

	fit = NULL
	results = NULL

	list(fit=fit, results=results, auc=auc, fx=fx, classes=ytest)
}




#' @title Learns and classifies using a Random Forest Classifier 
#' @description This function trains a Random Forest classifier on the training set and tests its accuracy on a test set. It reports the fit, prediction results, area under curve (AUC accuracy) on test, probability of belonging to predicted class, and original labels of test
#' @param Xtrain Predictor variables of training set
#' @param ytrain Response variables of training set
#' @param Xtest Predictor variables of test set
#' @param ytest Response variables of test set
#' @return A list with model fit, test results, AUC, probability per class, original labels
#' @export
test_rf = function(Xtrain, ytrain, Xtest, ytest)
{

	if(! is.factor(ytrain))
		ytrain = factor(ytrain)

	fit = randomForest(x=Xtrain, y=ytrain)
	results = predict(fit, Xtest, type="prob")

	fx = results[, 2]
	auc = -1
	tryCatch({
		auc = 0 + auc(ytest, fx)
	}, error = function(e){ })

	fit = NULL
	results = NULL

	list(fit=fit, results=results, auc=auc, fx=fx, classes=ytest)
}



#' @title Learns and classifies using a Linear Discriminant Analysis classifier 
#' @description This function trains a LDA classifier on the training set and tests its accuracy on a test set. It reports the fit, prediction results, area under curve (AUC accuracy) on test, probability of belonging to predicted class, and original labels of test
#' @param Xtrain Predictor variables of training set
#' @param ytrain Response variables of training set
#' @param Xtest Predictor variables of test set
#' @param ytest Response variables of test set
#' @return A list with model fit, test results, AUC, probability per class, original labels
#' @export
test_lda = function(Xtrain, ytrain, Xtest, ytest)
{

	# remove features with zero variance
	if( sum(apply(Xtrain, 2, var) == 0) > 0 )
	{
		Xtest = Xtest[, -which(apply(Xtrain, 2, var) == 0) ]
		Xtrain = Xtrain[, -which(apply(Xtrain, 2, var) == 0) ]
	}

	fit = lda(Xtrain, grouping=ytrain)
	results = predict(fit , Xtest)

	fx = results$x
	auc = -1
	if(length(fx) > 2 && length(table(ytest)) >= 2)
		auc = 0 + auc(ytest, fx)

	list(fit=fit, results=results, auc=auc, fx=fx, classes=ytest)
}


#' @title Automatic differential gene expression using limma 
#' @description This function tests for differential expression per gene using a linear model and definition contrasts that splits the input matrix into groups of samples. It outputs the fit, coefficients per gene, as well as the sum of residual squares after accounting for weights.
#' @param mat A gene expression matrix (genes x samples)
#' @param model A model of contrasts splitting samples
#' @param weights Assign differential weights either across genes or samples
#' @return It outputs a list with the model fit, coefficients per gene, as well as the sum of residual squares after accounting for differential weights
#' @export  
test_limma_lm = function(mat,mod,weights)
{
  fitln <- lmFit(log(mat),mod,weights=weights)
  b = coefficients(fitln)
  res = residuals(fitln,log(mat))
  
  sum_sq_res = sapply(seq(nrow(res)),function(i){ sum(res[i,which(weights[i,])]^2,na.rm=TRUE)/fitln$df[i]})

  rownames(fitln) = rownames(mat)
  
  list(fit=fitln, coef=b, ssr=sum_sq_res)
}
