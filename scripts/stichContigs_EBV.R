library(data.table)
library(seqinr)


getArgs <- function() {
    myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
      myargs <- lapply(myargs.list,function(x) x[2] )
      names(myargs) <- lapply(myargs.list,function(x) x[1])
      return (myargs)
  }

myArgs <- getArgs()

print(myArgs)

if ('arg1' %in% names(myArgs)) arg1 <- as.character( myArgs[[ 'arg1' ]])
if ('arg2' %in% names(myArgs)) arg2 <- as.character( myArgs[[ 'arg2' ]])
if ('arg3' %in% names(myArgs)) arg3 <- as.character( myArgs[[ 'arg3' ]])
if ('arg4' %in% names(myArgs)) arg4 <- as.character( myArgs[[ 'arg4' ]])

print(arg1)
print(arg2)
print(arg3)
print(arg4)


genome.recov<-matrix(0, nrow=1, ncol=2)
colnames(genome.recov)<-c("type1", "type2")

type.length<-list()
type.length[["type1"]] <- 171823
type.length[["type2"]] <- 172764

pseudogenome.L<-list()
pseudogenome.L[["type1"]] <- NULL
pseudogenome.L[["type2"]] <- NULL

final.contigs.L<-list()
final.contigs.L[["type1"]] <- NULL
final.contigs.L[["type2"]] <- NULL

totalLength<-list()
totalLength[["type1"]] <- NULL
totalLength[["type2"]] <- NULL

clist<-c(arg1, arg2)


for (file in clist){
  print(file)
  if (grepl(pattern='type1', x=file)){
    type<-'type1'}
  else {
    type<-'type2'
  }
    
  message("working with type", type)
    
  ordered.contigs<-fread(file, sep="\t", select=c(1,2,3,4))
  setnames(ordered.contigs, old=c("id","start", "end", "seq" ))
  ordered.contigs.df<-as.data.frame(ordered.contigs)

  test<-ordered.contigs.df[,1:3]

  for (i in 1:(nrow(ordered.contigs.df)-1)){
    if ( (ordered.contigs.df[i+1,"start"]>=ordered.contigs.df[i,"start"])  &  (ordered.contigs.df[i+1,"end"] <=ordered.contigs.df[i, "end"]) ) {
      ordered.contigs.df[i+1,]<-ordered.contigs.df[i,]
    } else {
      ordered.contigs.df[i+1,]<-ordered.contigs.df[i+1,]
    }
  }

  ordered.contigs.df<-unique(ordered.contigs.df)


  overlap<-matrix(0, nrow=nrow(ordered.contigs.df)-1, ncol=1)
  updated.contigs<-ordered.contigs.df

  for (i in 1:(nrow(updated.contigs)-1)) {
    print(i)
    overlap[i,]<-updated.contigs[i, "end"]-updated.contigs[i+1,"start"]
  
    if (overlap[i,]>=0) {
      message("overlap")
      print(overlap[i,])
      updated.contigs[i,"id"]<-paste(updated.contigs[i,"id"], updated.contigs[i+1,"id"], sep="")
      updated.contigs[i,"start"]<-updated.contigs[i, "start"]
      updated.contigs[i,"end"]<-updated.contigs[i+1, "end"]
      updated.contigs[i, "seq"]<-paste(substr(updated.contigs[i,"seq"], start=1, stop=(nchar(updated.contigs[i,"seq"])-overlap[i,] -1) ), updated.contigs[i+1,"seq"], sep="")

      updated.contigs[i+1,]<-updated.contigs[i,]
    } else {
      message("gap")
      print(overlap[i,])
      updated.contigs[i,"id"]<-updated.contigs[i,"id"]
      updated.contigs[i,"start"]<-updated.contigs[i,"start"]
      updated.contigs[i,"end"]<-updated.contigs[i,"end"]
      updated.contigs[i,"seq"]<-updated.contigs[i,"seq"]

    }
  
  }

  updated.contigs<-unique(updated.contigs)

  keep<-vector()

  for (i in 1:(nrow(updated.contigs))) {
    print(i)
    keep[i]<-all(updated.contigs[i, "start"]<updated.contigs[i+1,"start"] & updated.contigs[i,"end"]<updated.contigs[i+1, "end"])
    keep[which(is.na(keep))]<-TRUE
  }


  final.contigs<- updated.contigs[keep, ]
  final.contigs$gap<-0

  if (nrow(final.contigs)!=1){  
    for (i in 1:(nrow(final.contigs)-1)) {
      final.contigs[i,"gap"]<-(final.contigs[i+1,"start"] - final.contigs[i,"end"])-1
    
      if (final.contigs[i,"gap"]!=0){
        final.contigs[i+1,"seq"]<- paste(final.contigs[i,"seq"], final.contigs[i+1,"seq"], sep=paste((rep("N",final.contigs[i,"gap"])), collapse=""))
      } else {
        final.contigs[i+1,"seq"]<- paste(final.contigs[i,"seq"], final.contigs[i+1,"seq"], sep="")      
      }
    
    }

  } else {
    final.contigs<-final.contigs
  }

  message("1stDone")



  pseudogenome.tmp<-final.contigs[nrow(final.contigs),"seq"]
  message("2ndDone")

  
  if (final.contigs[1,"start"]==1){
    pseudogenome<-pseudogenome.tmp
  } else {
    Ns<-paste(rep("N", (final.contigs[1,"start"] - 1)), collapse="")
    pseudogenome<-paste( Ns,pseudogenome.tmp, sep="")
  }

  pseudogenome.L[[type]]<-pseudogenome
    message("3rdDone")

  final.contigs.L[[type]]<-final.contigs
    message("4thDone")

  totalLength[[type]]<-sum(final.contigs.L[[type]][,3]-final.contigs.L[[type]][,2])
  print(type)
  print(totalLength)
  message("5thDone")

  print(genome.recov[1,type])
  print(totalLength[[type]])
  print(type.length[[type]])
  genome.recov[1,type]<-totalLength[[type]]/type.length[[type]]

  message("6thDone")

}

  save.image( paste(arg3,"/pseudogenome.RData", sep=""))
  message("7Done")


if (genome.recov[,"type1"] > genome.recov[,"type2"]) {
  message("type1 genotype")

  print(genome.recov)
  write.table(x="type1", file=paste(arg3, "/genotype.tab", sep=""), quote=F, row.names=F, col.names=F)

  write.fasta(sequences=pseudogenome.L[["type1"]], names=arg4, file.out=paste(arg3,"/pseudogenome.fasta", sep=""),   as.string = TRUE)
  write.table(file=paste(arg3,"/totalLength.tab", sep=""), x=totalLength[["type1"]], quote=F)

  write.table(file=paste(arg3, "/numberofContigs.tab", sep=""), x=nrow(final.contigs.L[["type1"]]), quote=F)

  write.table(file=paste(arg3,"/", arg4,"_contigs.tab", sep=""), x=final.contigs.L[["type1"]][,1:3], sep="\t", quote=F)

  write.table(file=paste(arg3,"/maxLength.tab", sep=""), x=max(final.contigs.L[["type1"]][,3]-final.contigs.L[["type1"]][,2]), quote=F)



}  else {
  message("type2 genotype")
  print(genome.recov)
  write.table(x="type2", file=paste(arg3, "/genotype.tab", sep=""), quote=F, row.names=F, col.names=F)

  write.fasta(sequences=pseudogenome.L[["type2"]], names=arg4, file.out=paste(arg3,"/pseudogenome.fasta", sep=""),   as.string = TRUE)
  write.table(file=paste(arg3,"/totalLength.tab", sep=""), x=totalLength[["type2"]], quote=F)

  write.table(file=paste(arg3, "/numberofContigs.tab", sep=""), x=nrow(final.contigs.L[["type2"]]), quote=F)

  write.table(file=paste(arg3,"/", arg4,"_contigs.tab", sep=""), x=final.contigs.L[["type2"]][,1:3], sep="\t", quote=F)

  write.table(file=paste(arg3,"/maxLength.tab", sep=""), x=max(final.contigs.L[["type2"]][,3]-final.contigs.L[["type2"]][,2]), quote=F)

      
}    


