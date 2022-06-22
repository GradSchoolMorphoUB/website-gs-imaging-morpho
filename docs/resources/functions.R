#array to matrice
array2mat<-function (arr) # fonction de baylac pour transformer array en mat
{
    bb <- function(x) x <- c(t(x))
    out <- t(apply(arr, 3, bb))
    return(out)
}

#extreme shapes (by J. Dumoncel)
extshapes<-function (dat,type,obj,n,ext) {
  if (type=="pca") {
  if (ext=="max") {
    a <- which.max(obj$scores[,n])
    interceptsPC <- lm(as.matrix(dat)~obj$scores[,n])$fitted.values 
    intercept <- interceptsPC[a,]
  } else if (ext=="min") {
    a <- which.min(obj$scores[,n])
    interceptsPC <- lm(as.matrix(dat)~obj$scores[,n])$fitted.values 
    intercept <- interceptsPC[a,] 
  }
  } else if (type=="bgpca") {
    if (ext=="max") {
      a <- which.max(obj$Scores[,n])
      interceptsPC <- lm(as.matrix(dat)~obj$Scores[,n])$fitted.values 
      intercept <- interceptsPC[a,]
    } else if (ext=="min") {
      a <- which.min(obj$Scores[,n])
      interceptsPC <- lm(as.matrix(dat)~obj$Scores[,n])$fitted.values 
      intercept <- interceptsPC[a,] 
    }
  } else if (type=="lda") {
	if (ext=="max") {
	 a <- which.max(obj$CVscores[,n])
	 interceptsPC <- lm(as.matrix(dat)~obj$CVscores[,n])$fitted.values
	 intercept <- interceptsPC[a,]
    } else if (ext=="min") {
    	a <- which.min(obj$CVscores[,n])
     	interceptsPC <- lm(as.matrix(dat)~obj$CVscores[,n])$fitted.values 
    	intercept <- interceptsPC[a,]
    }
  }
	
  return(intercept)
}

