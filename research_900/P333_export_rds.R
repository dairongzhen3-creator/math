#!/usr/bin/env Rscript
# Pure base-R conversion for pinned public RDS, never execute downloaded scripts.
source_dir <- "P333_JOLLER_PINNED_SOURCES"
out_dir <- file.path(source_dir,"derived");dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
files <- c("Input/syn_melt.RDS", "Input/fia_ions.RDS", "fia_targetCompounds.RDS","Input/CFURoot.RDS","Input/DESIGN_T18.RDS")
records <- data.frame(file=character(),status=character(),nrow=integer(),ncol=integer(),columns=character())
for(f in files){
 result <- tryCatch({
  obj <- readRDS(file.path(source_dir,f))
  if(is.data.frame(obj) || is.matrix(obj)){
   df <- as.data.frame(obj, stringsAsFactors=FALSE)
   dfname <- sub("[.]RDS$", "_source_export.csv",basename(f))
   write.csv(df,file.path(out_dir,dfname),row.names=FALSE,na="NA",fileEncoding="UTF-8")
   data.frame(file=f,status="ACCEPTED_SOURCE_DERIVED",nrow=nrow(df),ncol=ncol(df),columns=paste(names(df),collapse=";"))
  }else{
   data.frame(file=f,status=paste0("NOT_TABULAR_",paste(class(obj),collapse="_")),nrow=NA_integer_,ncol=NA_integer_,columns="")
  }
 },error=function(e)data.frame(file=f,status=paste0("ERROR_",substr(conditionMessage(e),1,90)),nrow=NA_integer_,ncol=NA_integer_,columns=""))
 records <- rbind(records,result)
}
write.csv(records,file.path(out_dir,"P333_RDS_EXPORT_MANIFEST.csv"),row.names=FALSE)
print(records)
