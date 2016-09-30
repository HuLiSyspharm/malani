#' G SVM models.
#' @description Returns accuracy performance of all genes. G support vector machine (SVM) classifiers trained using G different data matrixes, are used to predict labels in test data. Models are ranked based on prediction performances.
#' @param dat.train Train data with G features and (k-1)*S/k samples. Parameter k comes from cross-validation scheme and is specified by user (default is 2).
#' @param lab.train Class labels for train data.
#' @param dat.test Test data with G features and S/k samples.
#' @param lab.test Class labels for test data.
#' @return Accuracy scores for models. Each model represents one gene.
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @export


Gsvmod <- function(dat.train,lab.train,dat.test,lab.test) {

  if(!is.null(rownames(dat.train))) rownames(dat.train) <- c()
  if(!is.null(rownames(dat.test))) rownames(dat.test) <- c()

  N <- nrow(dat.train)
  genes <- seq_len(N)
  Accs <- matrix(0,N,1)
  Accs[,1] <- genes
  for (i in genes) {
    ID <- cbind(rep(i,(N-1)),genes[!genes %in% i])
    tmp <- dat.train[ID[,1],]*dat.train[ID[,2],]
    df1 <- data.frame(t(tmp), group = lab.train)
    model <- svm(group~.,df1,type="C-classification")

    tmp <- dat.test[ID[,1],]*dat.test[ID[,2],]
    df2 <- data.frame(t(tmp), group = lab.test)
    svm.pred <- predict(model,df2)
    Accs[i] <- ((sum(svm.pred == lab.test)/length(svm.pred))*100)
  }

  return(Accs)
}
