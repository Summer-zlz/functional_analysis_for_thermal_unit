# Fit/predict interface preserving FLCM.select.huber's objective and PGD update.
# No scaling; intercept and named load variable are unpenalized.
fhuber_prepare <- function(mydata, nbasis1=7L, nbasis2, load_name='mw') {
  stopifnot(names(mydata)[1:2]==c('id','time'), !anyDuplicated(names(mydata)),
            all(is.finite(as.matrix(mydata))))
  features <- names(mydata)[-(1:2)]
  if (!load_name %in% features) stop('Missing load covariate: ',load_name)
  if (nbasis1<4 || nbasis2<4) stop('Cubic splines need at least 4 bases')
  rng <- range(mydata$time)
  basis <- fda::create.bspline.basis(rng,nbasis2,4,
                                   seq(rng[1],rng[2],length.out=nbasis2-2))
  intercept_basis <- fda::create.bspline.basis(rng,nbasis1,4,
                                   seq(rng[1],rng[2],length.out=nbasis1-2))
  groups <- rep(seq_along(features),each=nbasis2)
  groups[rep(features==load_name,each=nbasis2)] <- 0L
  B <- fda::eval.basis(mydata$time,intercept_basis)
  S <- fda::eval.basis(mydata$time,basis)
  covariates <- as.matrix(mydata[,features])
  raw <- Matrix::Matrix(do.call(cbind,c(list(1,B),lapply(seq_along(features),
                          function(j) S*covariates[,j]))),sparse=TRUE)
  list(features=features,load_name=load_name,basis=basis,intercept_basis=intercept_basis,
       nbasis1=nbasis1,nbasis2=nbasis2,range=rng,in_mat=fda::inprod(basis,basis),
       deriv_mat=fda::inprod(basis,basis,fda::int2Lfd(2),fda::int2Lfd(2)),
       B=B,S=S,covariates=covariates,raw=raw,
       group=c(rep(0L,nbasis1+1L),groups),training_ids=sort(unique(mydata$id)))
}

fhuber_design <- function(prepared,psi,mydata=NULL,M=NULL) {
  R <- NULL
  if (is.null(M)) {
    R <- chol(prepared$in_mat+psi*prepared$deriv_mat)
    M <- solve(R)
  }
  if (is.null(mydata)) {
    B <- prepared$B; S <- prepared$S; covariates <- prepared$covariates
  } else {
    if (!identical(names(mydata)[-(1:2)],prepared$features)) stop('Feature order mismatch')
    if (!all(is.finite(as.matrix(mydata)))) stop('Nonfinite new data')
    if (any(mydata$time<prepared$range[1] | mydata$time>prepared$range[2]))
      stop('New times outside training basis domain; extrapolation refused')
    B <- fda::eval.basis(mydata$time,prepared$intercept_basis)
    S <- fda::eval.basis(mydata$time,prepared$basis)
    covariates <- as.matrix(mydata[,prepared$features])
  }
  SM <- S %*% M
  X <- do.call(cbind,c(list(1,B),lapply(seq_along(prepared$features),
                                      function(j) SM*covariates[,j])))
  if (!all(is.finite(X))) stop('Nonfinite design')
  operator <- NULL
  if (is.null(mydata) && !is.null(prepared$raw)) {
    transform <- Matrix::bdiag(c(list(Matrix::Diagonal(prepared$nbasis1+1L)),
                                 rep(list(M),length(prepared$features))))
    operator <- list(raw=prepared$raw,transform=transform,nrow=nrow(X),ncol=ncol(X))
    if(prepared$nbasis2>20L && !is.null(R)) {
      # K is banded for local B-splines. R^{-1} is dense, whereas R is sparse.
      # Solving R b=gamma and R^T g=g_raw preserves the same transformation.
      triangular <- as(Matrix::bdiag(c(list(Matrix::Diagonal(prepared$nbasis1+1L)),
                         rep(list(Matrix::Matrix(R,sparse=TRUE)),length(prepared$features)))),
                       'triangularMatrix')
      operator$triangular <- triangular
      operator$triangular_t <- Matrix::t(triangular)
    }
  }
  list(X=X,M=M,operator=operator)
}

fhuber_forward <- function(X,beta) {
  mapped <- if(is.null(X$triangular)) X$transform %*% beta else
    Matrix::solve(X$triangular,beta)
  X$raw %*% mapped
}
fhuber_reverse <- function(X,score) {
  raw_gradient <- Matrix::crossprod(X$raw,score)
  if(is.null(X$triangular)) Matrix::crossprod(X$transform,raw_gradient) else
    Matrix::solve(X$triangular_t,raw_gradient)
}

fhuber_initial_beta <- function(X,y,group,initialization='ones') {
  nc <- if(is.list(X)) X$ncol else ncol(X)
  if(initialization=='ones') return(rep(1,nc))
  stopifnot(initialization %in% c('training_unpenalized_lm','training_unpenalized_huber'))
  idx <- which(group==0L)
  unpenalized <- if(is.list(X)) as.matrix(X$raw %*% X$transform[,idx,drop=FALSE]) else X[,idx,drop=FALSE]
  coefficients <- stats::lm.fit(unpenalized,y)$coefficients
  # B-spline intercept includes a redundant constant. QR excludes dependent columns.
  coefficients[is.na(coefficients)] <- 0
  if(initialization=='training_unpenalized_huber') {
    for(it in seq_len(1000L)) {
      residual <- y-drop(unpenalized %*% coefficients)
      weights <- pmin(1,1.5/pmax(abs(residual),.Machine$double.eps))
      updated <- stats::lm.wfit(unpenalized,y,weights)$coefficients
      updated[is.na(updated)] <- 0
      change <- sqrt(sum((updated-coefficients)^2)); coefficients <- updated
      if(change<=1e-10) break
    }
  }
  beta <- numeric(nc); beta[idx] <- coefficients
  if(!all(is.finite(beta))) stop('Nonfinite training-only initialization')
  beta
}

fhuber_pgd <- function(X,y,group,alpha,lambda,delta=1.5,epsilon=1e-6,max_iter=10000L,
                       initialization='ones') {
  sparse <- is.list(X)
  nr <- if(sparse) X$nrow else nrow(X); nc <- if(sparse) X$ncol else ncol(X)
  stopifnot(length(y)==nr,length(group)==nc,all(is.finite(y)),
            is.finite(alpha),alpha>0,is.finite(lambda),lambda>=0)
  beta <- fhuber_initial_beta(X,y,group,initialization)
  blocks <- split(which(group!=0L),group[group!=0L])
  block_matrix <- do.call(cbind,blocks)
  status <- 'max_iter'; step_norm <- Inf
  for (iteration in seq_len(max_iter)) {
    fitted <- if(sparse) as.numeric(fhuber_forward(X,beta)) else drop(X %*% beta)
    residual <- y-fitted
    clipped <- pmax(-delta,pmin(delta,residual))
    gradient <- if(sparse) -as.numeric(fhuber_reverse(X,clipped))/nr else -drop(crossprod(X,clipped))/nr
    next_beta <- beta-alpha*gradient
    norms <- sqrt(colSums(matrix(next_beta[block_matrix],nrow=nrow(block_matrix))^2))
    scale <- numeric(length(norms)); positive <- norms>0
    scale[positive] <- pmax(1-lambda*alpha/norms[positive],0)
    next_beta[block_matrix] <- next_beta[block_matrix]*rep(scale,each=nrow(block_matrix))
    step_norm <- sqrt(sum((next_beta-beta)^2)); beta <- next_beta
    if (!all(is.finite(beta)) || !is.finite(step_norm)) {status <- 'nonfinite'; break}
    if (step_norm<=epsilon) {status <- 'converged'; break}
  }
  fitted <- if(sparse) as.numeric(fhuber_forward(X,beta)) else drop(X %*% beta)
  residual <- y-fitted
  objective <- mean(ifelse(abs(residual)<=delta,residual^2/2,
                          delta*abs(residual)-delta^2/2)) +
    lambda*sum(vapply(blocks,function(idx) sqrt(sum(beta[idx]^2)),numeric(1)))
  list(gamma=beta,status=status,converged=status=='converged',iterations=iteration,
       step_norm=step_norm,objective=objective,delta=delta,epsilon=epsilon,max_iter=max_iter,
       initialization=initialization)
}
fhuber_pgd <- compiler::cmpfun(fhuber_pgd)

# Columns are independent candidates; batching changes only matrix arithmetic.
fhuber_pgd_batch <- function(X,y,group,alphas,lambdas,delta=1.5,
                              epsilon=1e-6,max_iter=10000L,initialization='ones',
                              checkpoints=NULL) {
  nr <- X$nrow; nc <- X$ncol; count <- length(alphas)
  stopifnot(length(y)==nr,length(group)==nc,length(lambdas)==count)
  beta <- matrix(rep(fhuber_initial_beta(X,y,group,initialization),count),nc,count)
  blocks <- do.call(cbind,split(which(group!=0L),group[group!=0L]))
  q <- nrow(blocks); ng <- ncol(blocks)
  indices <- as.vector(blocks)
  status <- rep('max_iter',count); iterations <- rep(max_iter,count)
  step_norm <- rep(Inf,count); active <- seq_len(count)
  snapshots <- list()
  summarize <- function(iteration) {
    fitted <- as.matrix(fhuber_forward(X,beta))
    residual <- y-fitted
    loss <- ifelse(abs(residual)<=delta,residual^2/2,delta*abs(residual)-delta^2/2)
    penalty <- matrix(sqrt(colSums(matrix(beta[indices,,drop=FALSE]^2,nrow=q))),ng)
    completed <- iterations
    completed[status=='max_iter'] <- iteration
    list(gamma=beta,status=status,iterations=completed,step_norm=step_norm,
         objective=colMeans(loss)+lambdas*colSums(penalty))
  }
  if(!is.null(checkpoints)) stopifnot(all(checkpoints>0),all(checkpoints<=max_iter))
  for (it in seq_len(max_iter)) {
    current <- beta[,active,drop=FALSE]
    residual <- y-as.matrix(fhuber_forward(X,current))
    clipped <- pmax(-delta,pmin(delta,residual))
    dim(clipped) <- dim(residual)
    gradient <- -as.matrix(fhuber_reverse(X,clipped))/nr
    next_beta <- current-sweep(gradient,2,alphas[active],'*')
    norms <- matrix(sqrt(colSums(matrix(next_beta[indices,,drop=FALSE]^2,nrow=q))),ng)
    threshold <- matrix(rep(alphas[active]*lambdas[active],each=ng),ng)
    scale <- matrix(0,ng,length(active)); nonzero <- norms>0
    scale[nonzero] <- pmax(1-threshold[nonzero]/norms[nonzero],0)
    next_beta[indices,] <- next_beta[indices,,drop=FALSE]*
                           matrix(rep(scale,each=q),nrow=q*ng)
    step_norm[active] <- sqrt(colSums((next_beta-current)^2))
    beta[,active] <- next_beta
    invalid <- !is.finite(step_norm[active]) | colSums(!is.finite(next_beta))>0
    done <- !invalid & step_norm[active]<=epsilon
    status[active[invalid]] <- 'nonfinite'; status[active[done]] <- 'converged'
    iterations[active[invalid | done]] <- it
    active <- active[!invalid & !done]
    if(it %in% checkpoints) snapshots[[as.character(it)]] <- summarize(it)
    if(!length(active)) break
  }
  result <- summarize(it)
  for(checkpoint in checkpoints) {
    key <- as.character(checkpoint)
    if(is.null(snapshots[[key]])) snapshots[[key]] <- result
  }
  result$snapshots <- snapshots
  result
}
fhuber_pgd_batch <- compiler::cmpfun(fhuber_pgd_batch)

fhuber_fit <- function(y,mydata,psi,alpha,lambda,nbasis1=7L,nbasis2,
                       load_name='mw',prepared=NULL,design=NULL,
                       initialization='training_unpenalized_huber',max_iter=10000L) {
  if (is.null(prepared)) prepared <- fhuber_prepare(mydata,nbasis1,nbasis2,load_name)
  if (is.null(design)) design <- fhuber_design(prepared,psi)
  fit <- fhuber_pgd(if(is.null(design$operator)) design$X else design$operator,
                    y,prepared$group,alpha,lambda,initialization=initialization,max_iter=max_iter)
  gamma_blocks <- matrix(fit$gamma[-seq_len(prepared$nbasis1+1)],nrow=prepared$nbasis2)
  selected <- colSums(abs(gamma_blocks))>0 | prepared$features==load_name
  coefficients <- design$M %*% gamma_blocks
  colnames(coefficients) <- prepared$features
  structure(c(fit,list(psi=psi,alpha=alpha,lambda=lambda,M=design$M,
       basis=prepared$basis,intercept_basis=prepared$intercept_basis,
       nbasis1=prepared$nbasis1,nbasis2=prepared$nbasis2,range=prepared$range,
       features=prepared$features,load_name=load_name,group=prepared$group,
       selected_features=prepared$features[selected],coefficients=coefficients,
       intercept_coefficients=fit$gamma[seq_len(prepared$nbasis1+1)],
       training_ids=prepared$training_ids)),class='fhuber_model')
}

predict.fhuber_model <- function(object,newdata,allow_unconverged=FALSE,...) {
  if (!object$converged && !allow_unconverged) stop('Prediction refused: model has not converged')
  if(!object$status %in% c('converged','max_iter') || !all(is.finite(object$gamma)))
    stop('Prediction refused: numerically invalid model')
  design <- fhuber_design(object,object$psi,newdata,object$M)
  prediction <- drop(design$X %*% object$gamma)
  if (!all(is.finite(prediction))) stop('Nonfinite predictions')
  prediction
}

fhuber_metrics <- function(actual,predicted) {
  stopifnot(length(actual)==length(predicted),all(is.finite(actual)),all(is.finite(predicted)))
  error <- predicted-actual
  c(MAE=mean(abs(error)),RMSE=sqrt(mean(error^2)),
    MAPE=if(any(actual==0)) NA_real_ else 100*mean(abs(error/actual)))
}
