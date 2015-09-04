
                                        #library(Biostrings)
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

contig.lengths<-matrix(0, nrow=8, ncol=2)
colnames(contig.lengths)<-c("H1N1", "H3N2")
genome.recov<-matrix(0, nrow=1, ncol=2)
colnames(genome.recov)<-c("H1N1", "H3N2")
H1N1.length <- 13200
H3N2.length <- 13500


arg4<-"H1N1"

if (arg4=="H1N1"){

  ordered.contigs<-fread(arg1, sep="\t", select=c(1,2,3,4, 5))
  setnames(ordered.contigs, old=c("id","start", "end", "seq", "seg" ))
  ordered.contigs.df<-as.data.frame(ordered.contigs)
  updated.contigs<-ordered.contigs.df



###H1N1
  H1N1.seg1<-"gi|758967842|ref|NC_026438.1|"
  H1N1.seg2<-"gi|758899361|ref|NC_026435.1|"
  H1N1.seg3<-"gi|758967835|ref|NC_026437.1|"
  H1N1.seg4<-"gi|758899355|ref|NC_026433.1|"
  H1N1.seg5<-"gi|758899363|ref|NC_026436.1|"
  H1N1.seg6<-"gi|758899359|ref|NC_026434.1|"
  H1N1.seg7<-"gi|758899349|ref|NC_026431.1|"
  H1N1.seg8<-"gi|758899352|ref|NC_026432.1|"

  H1N1.segs<-list(H1N1.seg1, H1N1.seg2, H1N1.seg3, H1N1.seg4, H1N1.seg5, H1N1.seg6, H1N1.seg7, H1N1.seg8)
  
  seg1<-updated.contigs[which(updated.contigs[,5]==H1N1.seg1),]
  seg2<-updated.contigs[which(updated.contigs[,5]==H1N1.seg2),]
  seg3<-updated.contigs[which(updated.contigs[,5]==H1N1.seg3),]
  seg4<-updated.contigs[which(updated.contigs[,5]==H1N1.seg4),]
  seg5<-updated.contigs[which(updated.contigs[,5]==H1N1.seg5),]
  seg6<-updated.contigs[which(updated.contigs[,5]==H1N1.seg6),]
  seg7<-updated.contigs[which(updated.contigs[,5]==H1N1.seg7),]
  seg8<-updated.contigs[which(updated.contigs[,5]==H1N1.seg8),]



  seg1.ord<-seg1[order(seg1[,"start"], seg1[,"end"]),]
  seg2.ord<-seg2[order(seg2[,"start"], seg2[,"end"]),]
  seg3.ord<-seg3[order(seg3[,"start"], seg3[,"end"]),]
  seg4.ord<-seg4[order(seg4[,"start"], seg4[,"end"]),]
  seg5.ord<-seg5[order(seg5[,"start"], seg5[,"end"]),]
  seg6.ord<-seg6[order(seg6[,"start"], seg6[,"end"]),]
  seg7.ord<-seg7[order(seg7[,"start"], seg7[,"end"]),]
  seg8.ord<-seg8[order(seg8[,"start"], seg8[,"end"]),]


  segments<-list(seg1.ord, seg2.ord, seg3.ord, seg4.ord, seg5.ord, seg6.ord, seg7.ord, seg8.ord )

for (j in 1:8 ){

    overlap<-rep(0, n=nrow(segments[[j]])-1)

  
  if (nrow(segments[[j]])==1) {
    
    segments[[j]]<-segments[[j]]
    
  } else if (nrow(segments[[j]])==0) {

        nocontig<-as.data.frame(matrix(0, ncol=5, nrow=1))
        colnames(nocontig)<-c("id", "start", "end", "seq", "seg")
        segments[[j]]<-nocontig
        segments[[j]]["seg"]<-H1N1.segs[[j]]



    
  } else {
    
    for (i in 1:(nrow(segments[[j]])-1)) {
      
      print(i)
      print(j)
      print(segments[[j]])

      overlap[i]<-segments[[j]][i, "end"]-segments[[j]][i+1,"start"]
      if ( overlap[i]>0 & (segments[[j]][i+1,"end"] > segments[[j]][i,"end"]) ) {
        message("overlap")
        print(overlap[i])
        segments[[j]][i,"id"]<-paste(segments[[j]][i,"id"], segments[[j]][i+1,"id"], sep="")
        segments[[j]][i,"start"]<-segments[[j]][i, "start"]
        segments[[j]][i,"end"]<-segments[[j]][i+1, "end"]
        segments[[j]][i, "seq"]<-paste(substr(segments[[j]][i,"seq"], start=1, stop=(nchar(segments[[j]][i,"seq"])-overlap[i] -1) ), segments[[j]][i+1,"seq"], sep="")
        segments[[j]][i+1,]<-segments[[j]][i,]

      
##       overlap[i,]<-segments[[j]][i, "end"]-segments[[j]][i+1,"start"]
##       if ( overlap[j,]>0 & (segments[[j]][i+1,"end"] > segments[[j]][i,"end"]) ) {
##         message("overlap")
##         print(overlap[j,])
##         segments[[j]][i,"id"]<-paste(segments[[j]][i,"id"], segments[[j]][i+1,"id"], sep="")
##         segments[[j]][i,"start"]<-segments[[j]][i, "start"]
##         segments[[j]][i,"end"]<-segments[[j]][i+1, "end"]
##         segments[[j]][i, "seq"]<-paste(substr(segments[[j]][i,"seq"], start=1, stop=(nchar(segments[[j]][i,"seq"])-overlap[j,] -1) ), segments[[j]][i+1,"seq"], sep="")
##         segments[[j]][i+1,]<-segments[[j]][i,]

     } else {

        segments[[j]][i+1,]<-segments[[j]][i,]
      }     
    
    }
    
   
  }

    
  segments[[j]]<-unique(segments[[j]])
  H1N1.segments<-segments


}




for (j in 1:8 ){
  final.row<-nrow(segments[[j]])
  
  H1N1.segments[[j]]<-H1N1.segments[[j]][final.row,] 
  contig.lengths[j,"H1N1"]<- as.numeric(segments[[j]][final.row, "end"]-segments[[j]][final.row, "start"]+1)
}

print(arg4)
print(contig.lengths)
print(sum(contig.lengths[,"H1N1"])/H1N1.length)

genome.recov[1,"H1N1"]<-sum(contig.lengths[,"H1N1"])/H1N1.length
  

}


arg4<-"H3N2"

if (arg4=="H3N2"){

  ordered.contigs<-fread(arg2, sep="\t", select=c(1,2,3,4, 5))
  setnames(ordered.contigs, old=c("id","start", "end", "seq", "seg" ))
  ordered.contigs.df<-as.data.frame(ordered.contigs)
  overlap<-matrix(0, nrow=nrow(ordered.contigs.df)-1, ncol=1)
  updated.contigs<-ordered.contigs.df

  H3N2.seg1<-"gi|73919059|ref|NC_007373.1|"
  H3N2.seg2<-"gi|73919148|ref|NC_007372.1|"
  H3N2.seg3<-"gi|73919133|ref|NC_007371.1|"
  H3N2.seg4<-"gi|73919206|ref|NC_007366.1|"
  H3N2.seg5<-"gi|73919146|ref|NC_007369.1|"
  H3N2.seg6<-"gi|73919135|ref|NC_007368.1|"
  H3N2.seg7<-"gi|73919151|ref|NC_007367.1|"
  H3N2.seg8<-"gi|73919211|ref|NC_007370.1|"


  H3N2.segs<-list(H3N2.seg1, H3N2.seg2, H3N2.seg3, H3N2.seg4, H3N2.seg5, H3N2.seg6, H3N2.seg7, H3N2.seg8)

  
  seg1<-updated.contigs[which(updated.contigs[,5]==H3N2.seg1),]
  seg2<-updated.contigs[which(updated.contigs[,5]==H3N2.seg2),]
  seg3<-updated.contigs[which(updated.contigs[,5]==H3N2.seg3),]
  seg4<-updated.contigs[which(updated.contigs[,5]==H3N2.seg4),]
  seg5<-updated.contigs[which(updated.contigs[,5]==H3N2.seg5),]
  seg6<-updated.contigs[which(updated.contigs[,5]==H3N2.seg6),]
  seg7<-updated.contigs[which(updated.contigs[,5]==H3N2.seg7),]
  seg8<-updated.contigs[which(updated.contigs[,5]==H3N2.seg8),]



  seg1.ord<-seg1[order(seg1[,"start"], seg1[,"end"]),]
  seg2.ord<-seg2[order(seg2[,"start"], seg2[,"end"]),]
  seg3.ord<-seg3[order(seg3[,"start"], seg3[,"end"]),]
  seg4.ord<-seg4[order(seg4[,"start"], seg4[,"end"]),]
  seg5.ord<-seg5[order(seg5[,"start"], seg5[,"end"]),]
  seg6.ord<-seg6[order(seg6[,"start"], seg6[,"end"]),]
  seg7.ord<-seg7[order(seg7[,"start"], seg7[,"end"]),]
  seg8.ord<-seg8[order(seg8[,"start"], seg8[,"end"]),]


  segments<-list(seg1.ord, seg2.ord, seg3.ord, seg4.ord, seg5.ord, seg6.ord, seg7.ord, seg8.ord )

  for (j in 1:8 ){
    overlap<-rep(0, n=nrow(segments[[j]])-1)

    if (nrow(segments[[j]])==1) {
    
      segments[[j]]<-segments[[j]]
    
    } else if (nrow(segments[[j]])==0) {

      nocontig<-as.data.frame(matrix(0, ncol=5, nrow=1))
      colnames(nocontig)<-c("id", "start", "end", "seq", "seg")
      segments[[j]]<-nocontig
      segments[[j]]["seg"]<-H3N2.segs[[j]]



    
    } else {
    
      for (i in 1:(nrow(segments[[j]])-1)) {
        print(i)

        overlap[i]<-segments[[j]][i, "end"]-segments[[j]][i+1,"start"]
        if ( overlap[i]>0 & (segments[[j]][i+1,"end"] > segments[[j]][i,"end"]) ) {
          message("overlap")
          print(overlap[i])
          segments[[j]][i,"id"]<-paste(segments[[j]][i,"id"], segments[[j]][i+1,"id"], sep="")
          segments[[j]][i,"start"]<-segments[[j]][i, "start"]
          segments[[j]][i,"end"]<-segments[[j]][i+1, "end"]
          segments[[j]][i, "seq"]<-paste(substr(segments[[j]][i,"seq"], start=1, stop=(nchar(segments[[j]][i,"seq"])-overlap[i] -1) ), segments[[j]][i+1,"seq"], sep="")
          segments[[j]][i+1,]<-segments[[j]][i,]


##         overlap[j,]<-segments[[j]][i, "end"]-segments[[j]][i+1,"start"]
##         if ( overlap[j,]>0 & (segments[[j]][i+1,"end"] > segments[[j]][i,"end"]) ) {
##           message("overlap")
##           print(overlap[j,])
##           segments[[j]][i,"id"]<-paste(segments[[j]][i,"id"], segments[[j]][i+1,"id"], sep="")
##           segments[[j]][i,"start"]<-segments[[j]][i, "start"]
##           segments[[j]][i,"end"]<-segments[[j]][i+1, "end"]
##           segments[[j]][i, "seq"]<-paste(substr(segments[[j]][i,"seq"], start=1, stop=(nchar(segments[[j]][i,"seq"])-overlap[j,] -1) ), segments[[j]][i+1,"seq"], sep="")
##           segments[[j]][i+1,]<-segments[[j]][i,]
        } else {

          segments[[j]][i+1,]<-segments[[j]][i,]
        }     
    
      }
    
   
    }

    
    segments[[j]]<-unique(segments[[j]])
    H3N2.segments<-segments


  }


  for (j in 1:8 ){
    final.row<-nrow(segments[[j]])
    H3N2.segments[[j]]<-H3N2.segments[[j]][final.row,] 

    contig.lengths[j,"H3N2"]<- as.numeric(segments[[j]][final.row, "end"]-segments[[j]][final.row, "start"]+1)
}

  print(arg4)
  print(contig.lengths)
  print(sum(contig.lengths[,"H3N2"])/H3N2.length)


  genome.recov[1,"H3N2"]<-sum(contig.lengths[,"H3N2"])/H3N2.length

}



if (genome.recov[,"H1N1"] > genome.recov[,"H3N2"]) {
  message("H1N1 genotype")
  print(genome.recov)
  print(contig.lengths)
  write.table(x="H1N1", file=paste(arg3, "/genotype.tab", sep=""), quote=F, row.names=F, col.names=F)
  for (j in 1:8){

    write.fasta(sequences=H1N1.segments[[j]]["seq"], names=H1N1.segments[[j]]["seg"], file.out=paste(arg3,"/H1N1segment", j, ".fasta", sep=""),   as.string = TRUE)
 }

}  else {
  message("H3N2 genotype")
  print(genome.recov)
  print(contig.lengths)
  write.table(x="H3N2", file=paste(arg3, "/genotype.tab", sep=""), quote=F, row.names=F, col.names=F)
  for (j in 1:8){
    
    write.fasta(sequences=H3N2.segments[[j]]["seq"], names=H3N2.segments[[j]]["seg"], file.out=paste(arg3,"/H3N2segment", j, ".fasta", sep=""),   as.string = TRUE)

  }
 

}


save.image("ordered_contigs.RData")
