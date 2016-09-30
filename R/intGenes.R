#' Select initial gene list from original data matrix.
#' @description Train G-1 SVM models in \code{k-fold} cross validation scheme to select initial genes list.
#' @param dat Original gene expression data matrix with G rows (number of genes) and S column (number of samples).
#' @param grp Class labels.
#' @param nfolds.out Outer cross validation number (default is 2).
# @param nfolds.in Inner cross validation number (default is 2).
#' @param top.per All genes are ranked based on their models performance and \code{ top.per\%} of them are selected as initial genes.
#' @return Selected initial genes.
#' @examples
#' data(malanidata)
#' int <- intGenes(dat,grp)
#' print(int$top.genes)
#' @export
#' @importFrom  e1071 svm



intGenes <- function(dat,grp,nfolds.out = 2,top.per = .05) {
  G <- nrow(dat)
  S <- ncol(dat)

  id.test <- rep(seq_len(nfolds.out),len = S)[sample(S)]
  id.test <- lapply(seq.int(1,nfolds.out), function(x) which(id.test == x))
  panels <- matrix(0,G,nfolds.out)
  for(i in 1:nfolds.out){

    dat.test <- dat[,id.test[[i]]]
    grp.test <- grp[id.test[[i]]]
    dat.train <- dat[,-id.test[[i]]]
    grp.train <- grp[-id.test[[i]]]

    panels[,i] <- Gsvmod(dat.train,grp.train,dat.test,grp.test)

  }
  panels.acc <- cbind(seq_len(G),apply(panels, 1, mean))
  panels.acc <- panels.acc[order(panels.acc[,2],decreasing = T),1]
  intG.out <- list(id.test = id.test,top.genes = panels.acc[1:floor(top.per*G)])
  return(intG.out)
}



