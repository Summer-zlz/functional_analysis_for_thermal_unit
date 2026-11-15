# Run from this folder: Rscript run_fhuber.R
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(Matrix))
source('fhuber_time_cv.R')

write_csv <- function(x,path) {
  x[] <- lapply(x,function(v) if(is.double(v)) sprintf('%.17g',v) else v)
  write.csv(x,path,row.names=FALSE,na='')
}
read_data <- function(label) {
  file <- paste0('coaldata_',tolower(label),'_unit1.csv')
  d <- if(label=='STABLE') read.csv(file,check.names=FALSE) else
    read.csv(file,colClasses=c(rep('numeric',20),'NULL'),check.names=FALSE)
  stopifnot(names(d)[1:3]==c('id','time','mw'),tail(names(d),1)=='coalseries',
            all(is.finite(as.matrix(d))),!anyDuplicated(d[,c('id','time')]))
  d[order(d$id,d$time),]
}
prepare <- function(d) fhuber_prepare(d[,-ncol(d),drop=FALSE],7L,ceiling((max(d$time)+1)/2))
fit_final <- function(d,best) fhuber_fit(d[[ncol(d)]],d[,-ncol(d),drop=FALSE],
  best$psi,1e-4,best$lambda,prepared=prepare(d),max_iter=500000L)
score <- function(model,d) {
  pred <- predict(model,d[,-ncol(d),drop=FALSE],allow_unconverged=TRUE)
  list(metrics=as.list(fhuber_metrics(d[[ncol(d)]],pred)),
       points=data.frame(original_id=d$id,time=d$time,actual=d[[ncol(d)]],predicted=pred,
                         error=pred-d[[ncol(d)]]))
}
main <- function() {
  dir.create('results',showWarnings=FALSE)
  all_cv <- all_best <- metrics <- direct <- features <- list()
  for(label in c('UP','DOWN','STABLE')) {
    d <- read_data(label); ids <- sort(unique(d$id)); n <- length(ids)
    span <- function(a,b) ids[seq.int(floor(a*n)+1L,floor(b*n))]
    windows <- list(c(0,.5,.5,.6),c(.1,.6,.6,.7),c(.2,.7,.7,.8))
    records <- list()
    for(k in 1:3) {
      w <- windows[[k]]; train_ids <- span(w[1],w[2]); val_ids <- span(w[3],w[4])
      stopifnot(!length(intersect(train_ids,val_ids)))
      tr <- d[d$id %in% train_ids,]; va <- d[d$id %in% val_ids,]; p <- prepare(tr)
      for(psi in 10^seq(-2,10,length.out=15)) {
        cat(label,'fold',k,'psi',psi,'\n'); flush.console()
        result <- tryCatch({
          design <- fhuber_design(p,psi)
          fit <- fhuber_pgd_batch(design$operator,tr[[ncol(tr)]],p$group,
            rep(1e-4,4),c(.01,.1,.5,1),max_iter=500000L,
            initialization='training_unpenalized_huber')
          xval <- fhuber_design(p,psi,va[,-ncol(va),drop=FALSE],design$M)$X
          pred <- xval %*% fit$gamma
          data.frame(dataset=label,fold=k,psi=psi,alpha=1e-4,lambda=c(.01,.1,.5,1),
            RMSE=vapply(1:4,function(j) fhuber_metrics(va[[ncol(va)]],pred[,j])['RMSE'],numeric(1)),
            status=fit$status,iterations=fit$iterations)
        },error=function(e) data.frame(dataset=label,fold=k,psi=psi,alpha=1e-4,
          lambda=c(.01,.1,.5,1),RMSE=NA_real_,status=conditionMessage(e),iterations=NA_integer_))
        records[[length(records)+1L]] <- result
        write_csv(do.call(rbind,records),file.path('results',paste0(label,'_fold_results.csv')))
      }
    }
    records <- do.call(rbind,records)
    cv <- do.call(rbind,lapply(split(records,interaction(records$psi,records$lambda)),function(r) {
      r <- r[order(r$fold),]; ok <- nrow(r)==3L && all(is.finite(r$RMSE)) &&
        all(r$status %in% c('converged','max_iter'))
      data.frame(dataset=label,psi=r$psi[1],alpha=1e-4,lambda=r$lambda[1],
        fold1_RMSE=r$RMSE[1],fold2_RMSE=r$RMSE[2],fold3_RMSE=r$RMSE[3],
        mean_validation_RMSE=if(ok) mean(r$RMSE) else Inf,eligible=ok)
    }))
    eligible <- cv[cv$eligible,]; stopifnot(nrow(eligible)>0)
    best <- eligible[order(eligible$mean_validation_RMSE,eligible$psi,eligible$lambda),][1,]
    all_cv[[label]] <- cv; all_best[[label]] <- best
    tr <- d[d$id %in% span(0,.8),]; te <- d[d$id %in% span(.8,1),]
    groups <- list(); groups[[label]] <- tr
    if(label=='STABLE') groups[['STABLE_full']] <- d else {
      membership <- read.csv(paste0('coaldata_',tolower(label),'_unit1.csv'))
      stopifnot(all(vapply(split(membership$cluster,membership$id),function(z)
        length(unique(z))==1L,logical(1))))
      for(g in sort(unique(membership$cluster)))
        groups[[paste0(label,g)]] <- d[d$id %in% membership$id[membership$cluster==g],]
    }
    for(name in names(groups)) {
      data <- groups[[name]]; model <- fit_final(data,best)
      saveRDS(model,file.path('results',paste0(name,'_model.rds')))
      features[[name]] <- data.frame(dataset=name,feature=model$selected_features)
      for(period in if(name==label) c('train','test') else 'all_sample_fit') {
        current <- if(period=='test') te else data; value <- score(model,current)
        row <- data.frame(dataset=if(name=='STABLE_full') 'STABLE' else name,
          period=period,trajectories=length(unique(current$id)),
          as.data.frame(value$metrics),status=model$status,iterations=model$iterations)
        write_csv(value$points,file.path('results',paste0(name,'_',period,'_predictions.csv')))
        if(name==label) metrics[[length(metrics)+1L]] <- row else direct[[name]] <- row
      }
    }
    write_csv(do.call(rbind,all_cv),'results/cv_results.csv')
    write_csv(do.call(rbind,all_best),'results/best_parameters.csv')
    write_csv(do.call(rbind,metrics),'results/metrics.csv')
    write_csv(do.call(rbind,direct),'results/direct_fit_metrics.csv')
    write_csv(do.call(rbind,features),'results/selected_features.csv')
  }
}
if(sys.nframe()==0L) main()
