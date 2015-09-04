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

print(arg1)
print(arg2)
print(arg3)


ordered.contigs<-fread(arg1, sep="\t", select=c(1,2,3,4))
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




pseudogenome.tmp<-final.contigs[nrow(final.contigs),"seq"]

if (final.contigs[1,"start"]==1){
   pseudogenome<-pseudogenome.tmp
 } else {
     Ns<-paste(rep("N", (final.contigs[1,"start"] - 1)), collapse="")
       pseudogenome<-paste( Ns,pseudogenome.tmp, sep="")
   }



write.fasta(sequences=pseudogenome, names=arg3, file.out=paste(arg2,"/pseudogenome.fasta", sep=""),   as.string = TRUE)



write.table(file=paste(arg2,"/totalLength.tab", sep=""), x=sum(final.contigs[,3]-final.contigs[,2]), quote=F)

write.table(file=paste(arg2, "/numberofContigs.tab", sep=""), x=nrow(final.contigs), quote=F)

write.table(file=paste(arg2,"/", arg3,"_contigs.tab", sep=""), x=final.contigs[,1:3], sep="\t", quote=F)

write.table(file=paste(arg2,"/maxLength.tab", sep=""), x=max(final.contigs[,3]-final.contigs[,2]), quote=F)


save.image( paste(arg2,"/pseudogenome.RData", sep=""))
