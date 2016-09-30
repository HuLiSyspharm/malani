#' Find best performing pairs
#' @description Combine each gene in initial set with all genes in the original set. Top \code{npair}  pairs are selected to construct the \code{Q} matrix.
#' @param X Original gene expression data matrix. With G rows (number of genes) and S column (number of samples).
#' @param LX Class labels.
#' @param theta Initial gene set.
#' @param npair Given a gene in initial set, top \code{npair} best performing pairs correspond to that gene are selected (Default is 10).
#' @return Best (\code{npair*G/20}) performing pairs.
# @examples
# data(malanidata)
# int <- intGenes(dat,grp)
# Q.pairs <- pairmod(dat,grp,int)
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @export


pairmod <- function(X,LX,theta,npair = 10) {
  G <- nrow(X)
  nf <- length(theta$id.test)
  gene <- seq_len(G)
  intG <- theta$top.genes
  Q.pairs <- c()
  for (i in 1:length(intG)){
    pairs.id <- cbind(rep(intG[i],len = (G-1)),setdiff(gene,intG[i]))
    tmp.acc <- matrix(0,(G-1),nf)
    for(j in 1:nf){

      df.train <- data.frame(t(X[,-theta$id.test[[j]]]), group = LX[-theta$id.test[[j]]])
      df.test <- data.frame(t(X[,theta$id.test[[j]]]), group = LX[theta$id.test[[j]]])
      colnames(df.train) <- c(gene,"group")
      colnames(df.test) <- c(gene,"group")

      tmp.acc[,j] <- sapply(1:(G-1), function(x) {
        df1 <- df.train[,c(pairs.id[x,],'group')]
        df2 <- df.test[,c(pairs.id[x,],'group')]
        svm.model <- svm(group~.,df1,type = "C-classification")
        svm.pred <- predict(svm.model,df2)
        100*(sum((svm.pred == df2$group))/length(df2$group))
      })
    }
    pairs.id <- cbind(pairs.id,apply(tmp.acc, 1, mean))
    pairs.id <- pairs.id[order(pairs.id[,ncol(pairs.id)],decreasing = T),]
    Q.pairs <- rbind(Q.pairs,pairs.id[1:npair,1:2])
    cat('\n', i*100/length(intG),' %', ' completed.')
  }
  return(Q.pairs)
}


