


#' Smooth the cell scores based on expression neighborhood
#'
#' @param cell_scores a score for each cell
#' @param NN_mat nearest neighbor matrix
#'
#' @return A smoothed cell score vector
#' @export
smoothCellScoresMatrix <- function(cell_scores, NN_mat) {
  # Check if the length of scores matches the number of vertices in the graph
  if(length(cell_scores) != dim(NN_mat)[1] ||
     length(cell_scores) != dim(NN_mat)[2]) {
    stop("The length of cell_scores must match the dim of NN_mat.")
  }

  # Convert cell_scores to a column matrix
  score_matrix <- matrix(cell_scores, nrow = length(cell_scores), ncol = 1)

  # Normalize adjacency matrix rows to sum to 1
  row_sums <- rowSums(NN_mat)
  row_sums[row_sums == 0] <- 1  # avoid division by zero for isolated cells
  NN_mat_norm <- NN_mat / row_sums

  # Matrix multiplication to smooth cell_scores
  smoothed_cell_scores <- NN_mat_norm %*% score_matrix

  # Flatten the matrix to return a vector
  return(as.vector(smoothed_cell_scores))
}


# DLG.multi.get.gene.pval<-function(cell.type,R){
#   b<-grepl("vs.",names(R))
#   pairs1<-get.strsplit(names(R),".vs.",1:2)
#   b1<-b&is.element(pairs1[,1],cell.type)
#   b2<-b&is.element(pairs1[,2],cell.type)
#   if((sum(b1)+sum(b2))==0){return(NULL)}
#   f1<-function(m1,pi = "p1"){
#     x<-m1[[pi]]
#     rownames(x)<-paste0(x$program,ifelse(x$up,".up_",".down_"),x$genes)
#     return(x)
#   }
#   m<-c(lapply(R[b1],f1),lapply(R[b2],function(m1) f1(m1,"p2")))
#   g<-unique(unlist(lapply(m,rownames)))
#   p<-get.strsplit(g,"_",1:2)
#   colnames(p)<-c("programF","genes")
#   rownames(p)<-g
#   p<-cbind.data.frame(p,program = get.strsplit(g,".",1),
#                       up = grepl("up",g))
#   names(m)<-gsub(paste0(cell.type,".vs."),"",names(m))
#   names(m)<-gsub(paste0(".vs.",cell.type),"",names(m))
#   for(i in names(m)){
#     x<-m[[i]]
#     g1<-paste0(x$program,ifelse(x$up,".up_",".down_"),x$genes)
#     idx<-match(g,g1)
#     p[,i]<-x$Z[idx]
#   }
#
#   # AD was without adjustments
#   p.up<-p.adjust.mat.per.label(get.pval.from.zscores(p[,5:ncol(p)]),p$programF)
#   p.down<-p.adjust.mat.per.label(get.pval.from.zscores(-p[,5:ncol(p)]),p$programF)
#
#   p.ub<-0.1
#   if(!is.matrix(p.up)){p.up<-as.matrix(p.up);p.down<-as.matrix(p.down)}
#   m<-cbind.data.frame(p[,c(names(m),"programF","genes")],
#                       p.up = fisher.combine(p.up),
#                       p.down = fisher.combine(p.down),
#                       n.up = rowSums(p.up<p.ub,na.rm = T),
#                       nf.up = rowMeans(p.up<p.ub,na.rm = T),
#                       n.down = rowSums(p.down<p.ub,na.rm = T),
#                       nf.down = rowMeans(p.down<p.ub,na.rm = T),
#                       p[,c("program","up")])
#   m$N<-m$n.up
#   m$N[!m$up]<-m$n.down[!m$up]
#   m$Nf<-m$nf.up
#   m$Nf[!m$up]<-m$nf.down[!m$up]
#   m$p.up[!m$up]<-1
#   m$p.down[m$up]<-1
#   return(m)
#
# }

testGeneScores <- function(object, sigmaChoice,
                           type = c("Mixed Effect", "GLM"),
                           covariates, frm){
  ## check input
  if (!is(object, "CoPro")) {
    stop("Input object must be a CoPro object")
  }

  if (!(type %in% c("Mixed Effect", "GLM"))) {
    stop("test type must be 'Mixed Effect' or 'GLM'. ")
  }

  ## write formula
  if(is.null(frm)){
    if(is.null(covariates)){
      frm = paste("y ~ x")
    }else{
      frm = paste("y ~ x +", paste(covariates, collapse = " + "))
    }
  }

  ## check if covariates exist in the meta.data
  meta_full = object@metaDataSub
  if(!all(covariates %in% colnames(meta_full))){
    stop("not all covariates present in the metadata")
  }


  cts <- object@cellTypesOfInterest
  nCC <- object@nCC

  sigmaName <- paste0("sigma_", sigmaChoice)
  ccNames <- paste0("CC_", seq_len(nCC))


  results <- setNames(vector(mode = "list", length = length(cts)), cts)
  for(i in cts){
    results[[i]] <- setNames(vector(mode = "list", length = nCC),
                                      ccNames)
  }

  if(type == "GLM"){

    for(i in cts){
      for(cci in ccNames){
        results[[i]][[cci]] <- testGeneGLM(object = object,
                                 sigmaName = sigmaName,
                                 cellTypeChoice = i,
                                 covariates = covariates,
                                 CCChoice = cci, frm = frm)
      }
    }


  }else{
    for(i in cts){
      for(cci in ccNames){
        results[[i]][[cci]] <- testGeneMixedEffect(object = object,
                                     sigmaName = sigmaName,
                                     cellTypeChoice = i,
                                     CCChoice = cci, frm = frm)
      }
    }
  }
  object@geneScoreTest = results

  return(object)
}

#' testGeneGLM
#'
#' @param object CoPro object
#' @param sigmaName sigmaName
#' @param cellTypeChoice cellTypeChoice
#' @param covariates covariates
#' @param CCChoice CCChoice
#' @param frm formula
#'
#' @returns lm results
#' @export
#' @importFrom stats lm
testGeneGLM <- function(
    object, sigmaName,cellTypeChoice,covariates,
    CCChoice,
    frm = "y ~ x + nCount_Spatial"
    ){

  meta = object@metaDataSub[object@cellTypesSub == cellTypeChoice, covariates, drop = FALSE]
  Y = getCellScores(object, sigma = as.numeric(gsub("sigma_", "", sigmaName)),
                    cellType = cellTypeChoice, ccIndex = as.numeric(gsub("CC_", "", CCChoice)),
                    verbose = FALSE)
  X = object@normalizedDataSub[object@cellTypesSub == cellTypeChoice,]

  m<-t(apply(X, MARGIN = 1,
             function(x){
               fit <- lm(formula = as.formula(frm), data = cbind(meta, x = x, y = Y))
               coef(summary(fit))["x", c("Estimate", "Pr(>|t|)")]
             }))
  return(m)
}

testGeneMixedEffect <- function(
    object, sigmaName,cellTypeChoice,
    CCChoice, frm = "y ~ (1 | samples) + x + nCount_Spatial"
    ){

        meta = object@metaDataSub[object@cellTypesSub == cellTypeChoice]
        Y = getCellScores(object, sigma = as.numeric(gsub("sigma_", "", sigmaName)), 
                          cellType = cellTypeChoice, ccIndex = as.numeric(gsub("CC_", "", CCChoice)), 
                          verbose = FALSE)
        X = object@normalizedDataSub[object@cellTypesSub == cellTypeChoice,]

        m<-t(apply(X, MARGIN = 1,
                   function(x){formula.HLM(Y,x,meta,formula = frm)}))

        colnames(m)<-c("Estimate","P")
        m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
        return(m)

}

get.onesided.p.value<-function(c,p){
  p[p==0] = min(p[p>0],na.rm = T)
  p.one.side <- p
  p.one.side[] <- NA
  b<-c>0&!is.na(c)
  p.one.side[b]=p[b]/2
  b<-c<=0&!is.na(c)
  p.one.side[b]=1-(p[b]/2)
  return(p.one.side)
}

get.p.zscores<-function(p){
  b<-p[,1]>0.5
  b[is.na(b)]<-F
  zscores<-(-log10(p[,1]))
  zscores[b]<-log10(p[b,2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}


get.cor.zscores<-function(c,p){
  v<-cbind(get.onesided.p.value(c,p),get.onesided.p.value(-c,p))
  z<-get.p.zscores(v)
  return(z)
}

apply.formula.HLM<-function(meta,X,Y,MARGIN = 1,formula = "y ~ (1 | samples) + x"){

  m<-t(apply(X,MARGIN = MARGIN,
             function(x){formula.HLM(Y,x,meta,formula = formula)}))

  colnames(m)<-c("Estimate","P")
  m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  return(m)
}

#' Regression based on a formula
#'
#' @param y y in regression
#' @param x x in regression
#' @param meta meta.data
#' @param formula formula
#' @param val value
#' @param return.all whether to return all
#'
#' @noRd
formula.HLM<-function(y,x,meta, formula = "y ~ (1 | samples) + x",
                      val = ifelse(is.numeric(x),"","TRUE"),return.all = F){
  meta$x<-x
  meta$y<-y
  f<-function(meta){
    M1 <- with(meta, lmer (formula = formula))
    if(return.all){
      c1<-summary(M1)$coef[,c("Estimate","Pr(>|t|)")]
    }else{
      c1<-summary(M1)$coef[paste0("x",val),]
      idx<-match(c("Estimate","Pr(>|t|)"),names(c1))
      c1<-c1[idx]
    }
    return(c1)
  }
  c1<-tryCatch({f(meta)},
               error = function(err){return(c(NA,NA))})
  return(c1)
}

