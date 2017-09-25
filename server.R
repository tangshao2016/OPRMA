options(shiny.maxRequestSize=200*1024^2) 
library(DT)

## input the data
load("D:\\research\\Collaborations\\Ruiqi\\distribution\\distribution\\hg19.gtf.input.wattr.RData")

load("D:\\research\\Collaborations\\Ruiqi\\distribution\\distribution\\mm10.gtf.input.wattr.RData")

load("D:\\research\\Collaborations\\Ruiqi\\distribution\\distribution\\hg19.exons.keep.first.isoform.proteincoding.RData")

load("D:\\research\\Collaborations\\Ruiqi\\distribution\\distribution\\mm10.exons.keep.first.isoform.proteincoding.RData")

find.peak.to.stopcodon.distance.by.strand<-function(data.matrix, stop.codon, pivot, strand){
  distance.to.pivot=0
  start.count=FALSE
  if(pivot<stop.codon){  ## pivot on the left, upstream
    for(i in 1:dim(data.matrix)[1]){
      sstart=data.matrix[i, 1]
      send=data.matrix[i, 2]
      if((sstart<pivot & pivot<send) & (sstart<stop.codon & stop.codon<send)){
        distance.to.pivot=stop.codon-pivot+1
        start.count=FALSE
        break
      }			
      if(sstart<=stop.codon & stop.codon<=send & start.count==TRUE){
        distance.to.pivot=distance.to.pivot+stop.codon-sstart+1
        start.count=FALSE
      }
      if(start.count==TRUE){
        distance.to.pivot=distance.to.pivot+(send-sstart+1)
      }
      
      if((sstart<=pivot & pivot<=send) & !(sstart<=stop.codon & stop.codon<=send)){
        distance.to.pivot=distance.to.pivot+send-pivot+1
        start.count=TRUE
      }
    }
    if(strand=="+"){
      distance.to.pivot=distance.to.pivot*(-1)
    }	
  }else{  ## pivot on the right
    for(i in 1:dim(data.matrix)[1]){
      sstart=data.matrix[i, 1]
      send=data.matrix[i, 2]
      if((sstart<pivot & pivot<send) & (sstart<stop.codon & stop.codon<send)){
        distance.to.pivot=pivot-stop.codon+1
        start.count=FALSE
        break
      }			
      if((sstart<pivot & pivot<send) & !(sstart<stop.codon & stop.codon<send)){ # reading the end, stop
        distance.to.pivot=distance.to.pivot+pivot-sstart+1
        start.count=FALSE
      }			
      if(start.count==TRUE){
        distance.to.pivot=distance.to.pivot+(send-sstart+1)
      }
      if(sstart<stop.codon & stop.codon<send){ # stop and pivot not on the same frame
        distance.to.pivot=distance.to.pivot+send-stop.codon+1
        start.count=TRUE
      }
      
    }
    if(strand=="-"){
      distance.to.pivot=distance.to.pivot*(-1)
    }
  }
  return(distance.to.pivot)
}



##----------------------
## new 01-17-2017 to adjust for strand and orientation  
find.peak.to.lastexon.distance.by.strand<-function(data.matrix, exon.pos, pivot, strand){
  distance.to.pivot=0
  start.count=FALSE
  if(pivot<exon.pos){
    for(i in 1:dim(data.matrix)[1]){
      sstart=data.matrix[i, 1]
      send=data.matrix[i, 2]
      if(sstart==exon.pos | send==exon.pos){  # first, determine when to stop
        start.count=FALSE
      }
      if(start.count==TRUE){
        distance.to.pivot=distance.to.pivot+(send-sstart+1)
      }
      if(sstart<=pivot & pivot<=send){ 
        distance.to.pivot=send-pivot+1
        if(send!=exon.pos){ # same exon, no need to set true
          start.count=TRUE
        }
      }
    }
    if(strand=="+"){
      distance.to.pivot=distance.to.pivot*(-1)
    }	
  }else{
    for(i in 1:dim(data.matrix)[1]){
      sstart=data.matrix[i, 1]
      send=data.matrix[i, 2]
      if(send==exon.pos | sstart==exon.pos){  # first, determine where to stat
        start.count=TRUE
      }			
      if(sstart<=pivot & pivot<=send){
        distance.to.pivot=distance.to.pivot+pivot-sstart+1
        start.count=FALSE
      }			
      if(start.count==TRUE & sstart==exon.pos){
        distance.to.pivot=distance.to.pivot+(send-sstart+1)
      }
    }
    if(strand=="-"){
      distance.to.pivot=distance.to.pivot*(-1)
    }
  }
  return(distance.to.pivot)
}





######### GET UTR
get.utr.side<-function(start.end.all, row.vector, types, isof){
  match.row.index=0
  for(i in 1:nrow(start.end.all)){
    if(all.equal(as.numeric(as.character(as.matrix(start.end.all[i,]))), as.numeric(as.character(as.matrix(row.vector))))==TRUE){
      match.row.index=i
      break
    }	
  }
  cds.index=which(types=='CDS')
  utr.index=which(types=='UTR')
  #if(length(cds.index)==0){print(isof)}
  utr.5.index=utr.index[which(utr.index<min(cds.index))]
  utr.3.index=utr.index[which(utr.index>max(cds.index))]
  utr.index=utr.5.index
  side=1
  if(length(which(utr.3.index==match.row.index))==1){
    side=0
    utr.index=utr.3.index
  }	
  return(list(side, utr.index))
}


#### binary search algorithm
binary_search_test<-function(target,data){
  min=0
  max=length(data)
  if(target>=min(data)&target<=max(data)){
    while(max-min>1){
      guess=floor((min+max)/2)
      if(data[guess]==target){
        return(guess)
      }
      else if(data[guess]<=target){
        min=guess
      }
      else if (data[guess]>=target){
        max=guess
      }
    }
    if(data[max]>=target& data[min]<=target){
      return(max)
    }
    else if(data[min]>=target){
      return(min)
    }
  }
  else{
    return(NA)
  }
}


#### find normalized position
find.normalized.position<-function(data.matrix, row.vector, pivot){
  if(dim(data.matrix)[1]>1){data.matrix=data.matrix[order(data.matrix[,1]),]}
  cds.length=0
  distance.to.pivot=0
  for(i in 1:dim(data.matrix)[1]){
    #if(length(all.equal(data.matrix[i,], row.vector)==TRUE)==1){
    if(all.equal(as.numeric(as.character(as.matrix(data.matrix[i,]))), as.numeric(as.character(as.matrix(row.vector))))==TRUE){
      distance.to.pivot=cds.length+abs(pivot-as.numeric(row.vector[1])+1)
    }
    cds.length=cds.length+(data.matrix[i,2]-data.matrix[i,1]+1)
  }
  distance=100*distance.to.pivot/cds.length
  if(distance>100){distance=100}
  return(distance)
}

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  dataInput <- reactive({
    infile<-input$file$datapath
    data.date="01212017"
    table.2.01212017=read.table(infile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(table.2.01212017)=c("Chr",     "Start",       "End",     "Sumit")
    table.2=table.2.01212017
    if(input$org=="hg19"){
      hg19.CDS.UTR.keep.first.isoform<-hg19.gtf.input.wattr[which((hg19.gtf.input.wattr$X16 %in% as.vector(unique(hg19.exons.keep.first.isoform.proteincoding$X16)))
                                                                  & (hg19.gtf.input.wattr$V3 == 'UTR' | hg19.gtf.input.wattr$V3 == 'CDS')), ]					
      table.1<-hg19.CDS.UTR.keep.first.isoform[,1:49]  ##  wrong, must use CDS+UTR, table.1=mm.exons.keep.first.isoform.proteincoding[,1:49]
      org.gtf.input.wattr<-hg19.gtf.input.wattr
      mean.utr.3.size.allgenes=1363.202
      mean.utr.5.size.allgenes=254.5384
      mean.cds.size.allgenes=1638.583
    }
    if(input$org=="mm10"){
      mm10.CDS.UTR.keep.first.isoform<-mm10.gtf.input.wattr[which((mm10.gtf.input.wattr$X16 %in% as.vector(unique(mm.exons.keep.first.isoform.proteincoding$X16)))
                                                                  & (mm10.gtf.input.wattr$V3 == 'UTR' | mm10.gtf.input.wattr$V3 == 'CDS')), ]					
      table.1<-mm10.CDS.UTR.keep.first.isoform[,1:49]  ##  wrong, must use CDS+UTR, table.1=mm.exons.keep.first.isoform.proteincoding[,1:49]
      org.gtf.input.wattr<-mm10.gtf.input.wattr
      mean.utr.3.size.allgenes=1134.835
      mean.utr.5.size.allgenes=218.4317
      mean.cds.size.allgenes=1572.418
    }
    
    runs=c("f1")
    for(run in runs){
      table.1.2.merged=matrix("", dim(table.1)[1], dim(table.1)[2]+dim(table.2)[2])
      m=1
      for(j in 1:length(table(table.2$Chr))){
        name=names(table(table.2$Chr))[j]
        table2<-table.2[which(table.2$Chr==name),]
        table1<-table.1[which(table.1$V1==name),]
        table1<-table1[order(table1$V4),]
        for(i in 1:dim(table2)[1]){
          target=table2[i,4]
          if(target>=min(table1$V4)& target<=max(table1$V4)){
            index1=binary_search_test(target = target,table1$V4)
            table3=table1[1:(index1-1),]
            table3=table3[order(table3$V5),]
            index2=binary_search_test(target = target,table3$V5)
            if(is.integer(index2)){
              matched.ids=row.names(table3[index2:dim(table3)[1],])
              for(k in matched.ids){
                table.1.2.merged[m,]=c(as.character(as.matrix(table1[which(row.names(table1)==k), ])), as.character(as.matrix(table2[i,])))
                m=m+1
              }	
            }
          }
        }
      }
      table.1.2.merged=table.1.2.merged[1:(m-1),]
      table.1.2.merged=data.frame(table.1.2.merged, stringsAsFactors=FALSE)
      assign(paste("table.1.2.merged.peak.", run, ".", data.date, sep=""), table.1.2.merged)
    }
    
    
    
    if(input$analysis=="STOPCODON"){
      distance.to.stop.codon=c()
      for(i in 1:dim(table.1.2.merged.peak.f1.01212017)[1]){
        peak.summit=as.numeric(table.1.2.merged.peak.f1.01212017[i, 'X53'])
        transcript.name=table.1.2.merged.peak.f1.01212017[i, 'X25']
        row.index=which(org.gtf.input.wattr$V3=='stop_codon' & org.gtf.input.wattr$X16==transcript.name)[1]
        if(!is.na(row.index)){
          stop.codon = org.gtf.input.wattr[row.index,'V4']
          strand=org.gtf.input.wattr[row.index, 'V7']
          
          start.end.all.exons=org.gtf.input.wattr[which(org.gtf.input.wattr$V3=='exon'  & org.gtf.input.wattr$X16==transcript.name), 4:5] #table.1[which(exons$X16==isof), 4:5]
          start.end.order.index=order(start.end.all.exons[,1])
          start.end.all.exons=start.end.all.exons[start.end.order.index,]
          
          distance.to.stop.codon=c(distance.to.stop.codon, 
                                   find.peak.to.stopcodon.distance.by.strand(start.end.all.exons, stop.codon, peak.summit, strand))
        }
      }
      human.01212017.distance.to.stop.codon=distance.to.stop.codon
      human.01212017.distance.to.stop.codon.chosen=human.01212017.distance.to.stop.codon[which(human.01212017.distance.to.stop.codon>(-1000)&human.01212017.distance.to.stop.codon<1000)]
    }
    
    
    if(input$analysis=="LASTEXON"){
      
      ## distance from peak summit to last exon junction
      isof.norm.dist.matrix=matrix("", nrow(table.1.2.merged.peak.f1.01212017), 4)
      distance.to.last.exon=c()
      for(i in 1:dim(table.1.2.merged.peak.f1.01212017)[1]){
        
        summit.peak=as.numeric(table.1.2.merged.peak.f1.01212017[i, 'X53'])  # c(4,5,54)  for cirseq (because the peak list is diff))
        
        isof=table.1.2.merged.peak.f1.01212017[i,'X25']
        strand=table.1.2.merged.peak.f1.01212017[i,'X7']
        
        start.end.all.exons=org.gtf.input.wattr[which(org.gtf.input.wattr$V3=='exon'  & org.gtf.input.wattr$X16==isof), 4:5] #table.1[which(exons$X16==isof), 4:5]
        start.end.order.index=order(start.end.all.exons[,1])
        start.end.all.exons=start.end.all.exons[start.end.order.index,]
        
        last.exon=c()
        if(strand=='-'){
          last.exon=start.end.all.exons[1,]
        }else{
          last.exon=start.end.all.exons[dim(start.end.all.exons)[1],]			
        }
        
        last.exon.pos=if(summit.peak>last.exon[2]){last.exon[2]}else if(summit.peak<last.exon[1]){last.exon[1]}else{0}
        last.exon.pos=as.numeric(last.exon.pos)
        
        normalized.distance.to.pivot=find.peak.to.lastexon.distance.by.strand(start.end.all.exons, last.exon.pos, summit.peak, strand)
        distance.to.last.exon=c(distance.to.last.exon, normalized.distance.to.pivot)
        
      }
      human.01212017.distance.to.last.exon=distance.to.last.exon
      human.01212017.distance.to.last.exon.chosen=human.01212017.distance.to.last.exon[which(human.01212017.distance.to.last.exon>(-1000)&human.01212017.distance.to.last.exon<1000)]
      #human.01212017.distance.to.last.exon.chosen=human.01212017.distance.to.last.exon.chosen[which(human.01212017.distance.to.last.exon.chosen!=0)]
    }
    
    
    if(input$analysis=="53UTRCDS"){
      
      ## need to make sure column 53 (or 54)  
      for(run in runs){  
        isof.norm.dist.matrix=c()	
        table.1.2.merged=get(paste("table.1.2.merged.peak.", run, ".", data.date, sep=""))
        unique.peaks=c()
        #for(isof in gene.peak.isof){
        #	isof.row.indices=which(table.1.2.merged$X25==isof)
        #	for(index in isof.row.indices){
        for(i in 1:dim(table.1.2.merged)[1]){
          #isof.info=table.1.2.merged[index, c(4,5,54)]  # c(4,5,54)  for cirseq (because the peak list is diff)
          index=i
          #isof.info=table.1.2.merged[index, c(4,5,ncol(table.1.2.merged))]  # c(4,5,54)  for cirseq (because the peak list is diff)
          isof.info=table.1.2.merged[index, c(4,5,53)]  # c(4,5,54)  for cirseq (because the peak list is diff)
          
          isof=table.1.2.merged[index,'X25']
          pivot=as.numeric(isof.info[3])
          this.type=table.1.2.merged[index, 3]
          transcript.type=table.1.2.merged[index, 'X21']
          unique.peak=paste(c(table.1.2.merged[index, c('X1' , 'X19')], pivot), collapse="-")  # c('X1' , 'X4', 'X5')
          if(length(which(unique.peaks==unique.peak))==0){		
            start.end.all=table.1[which(table.1$X16==isof), 4:5]
            start.end.order.index=order(start.end.all[,1])
            start.end.all=start.end.all[start.end.order.index,]
            types=table.1[which(table.1$X16==isof), 3]
            types=types[start.end.order.index]
            
            #start.end=table.1[which(table.1$X16==isof & table.1$V3==type), 4:5]
            #start.end=start.end[order(start.end[,1]),]
            strand=unique(table.1[which(table.1$X16==isof), 'V7'])
            row.vector=data.matrix(isof.info[1:2])
            UTR.5=1
            start.end=c()
            if(this.type=="UTR"){
              side.and.utr.index=get.utr.side(start.end.all, row.vector, types, isof)
              UTR.5=side.and.utr.index[[1]]
              utr.index=side.and.utr.index[[2]]
              start.end=start.end.all[utr.index,]
            }else{
              start.end=start.end.all[which(types=='CDS'),]
            }
            data.matrix=start.end
            if(nrow(data.matrix)>0){
              normalized.distance.to.pivot=find.normalized.position(data.matrix, row.vector, pivot)	
              if(strand=="-"){
                UTR.5=1-UTR.5
                normalized.distance.to.pivot=100-normalized.distance.to.pivot
              }
              isof.norm.dist.matrix=rbind(isof.norm.dist.matrix, c(isof, this.type, UTR.5, normalized.distance.to.pivot, transcript.type))
            }
          }
          unique.peaks=c(unique.peaks, unique.peak)
        }
        name=paste("isof.norm.dist.matrix.unique", run, ".", data.date, sep=".")
        assign(name, isof.norm.dist.matrix)
      }
      
      #name=paste("isof.norm.dist.matrix.", run, sep=".")
      #name=paste("isof.norm.dist.matrix.unique.", run, sep=".")  #OLD VERSION. TO BE DELETED.
      name=paste("isof.norm.dist.matrix.unique", run, ".", data.date, sep=".")
      isof.norm.dist.matrix=get(name)
      
      UTR.5.list=as.numeric(as.character(as.matrix(isof.norm.dist.matrix[which(isof.norm.dist.matrix[,5]=="protein_coding;" & isof.norm.dist.matrix[,2]=="UTR" & isof.norm.dist.matrix[,3]=="1"), 4])))
      UTR.3.list=as.numeric(as.character(as.matrix(isof.norm.dist.matrix[which(isof.norm.dist.matrix[,5]=="protein_coding;" & isof.norm.dist.matrix[,2]=="UTR" & isof.norm.dist.matrix[,3]=="0"), 4])))
      CDS.list  =as.numeric(as.character(as.matrix(isof.norm.dist.matrix[which(isof.norm.dist.matrix[,5]=="protein_coding;" & isof.norm.dist.matrix[,2]=="CDS"), 4])))
      
      peak.normalized.list=c(UTR.5.list*mean.utr.5.size.allgenes/100, 
                             mean.utr.5.size.allgenes+CDS.list*mean.cds.size.allgenes/100, 
                             mean.utr.5.size.allgenes+mean.cds.size.allgenes+UTR.3.list*mean.utr.3.size.allgenes/100)
      #hist(peak.normalized.list, breaks=20, prob=TRUE, col="grey")
      #lines(density(peak.normalized.list, adjust=2), type="l", lty="dotted", col="darkgreen", lwd=2) 
      
      
      dvals=density(peak.normalized.list, adjust=0.5)
      x=dvals$x[which(dvals$x>=0 & dvals$x<=3300)]
      y=dvals$y[which(dvals$x>=0 & dvals$x<=3300)]
      x=c(0,x)
      y=c(0,y)
      
    }
    
    
    
    if(input$analysis=="STOPCODON"){
      plotdata<-human.01212017.distance.to.stop.codon.chosen
    }
    else if(input$analysis=="LASTEXON"){
      plotdata<-human.01212017.distance.to.last.exon.chosen
    }
    else if(input$analysis=="53UTRCDS"){
      plotdata<-peak.normalized.list
    }
    list(plotdata=plotdata,mydata=table.1.2.merged)
    
  })
  
  output$rect1 <- renderPlot({
  dat<-dataInput()$mydata 
  sub<-dat[which(as.character(strsplit(dat$X19,";"))==input$gene),c(19,51:52)]
  sub$X51<-as.numeric(sub$X51)
  sub$X52<-as.numeric(sub$X52)
  plot(c(min(sub$X51)-50,max(sub$X52)+50),c(0,2),type = "n",xlab = "location",ylab = "",yaxt='n')
  y<-seq(1,10,by=1)
  for(i in 1:dim(sub)[1]){
    rect(sub$X51[i],0,sub$X52[i],0.5,border=rainbow(10)[i],lty = y[i],lwd = y[i])
  }
  
  })
  
  
  index<-0
  
 
  dataInput2 <- reactive({
    if(!is.null(input$file2)){
    infile<-input$file2$datapath
    data.date="01212017"
    
    table.2.01212017=read.table(infile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    index<-1
    colnames(table.2.01212017)=c("Chr",     "Start",       "End",     "Sumit")
    table.2=table.2.01212017
    if(input$org=="hg19"){
      hg19.CDS.UTR.keep.first.isoform<-hg19.gtf.input.wattr[which((hg19.gtf.input.wattr$X16 %in% as.vector(unique(hg19.exons.keep.first.isoform.proteincoding$X16)))
                                                                  & (hg19.gtf.input.wattr$V3 == 'UTR' | hg19.gtf.input.wattr$V3 == 'CDS')), ]					
      table.1<-hg19.CDS.UTR.keep.first.isoform[,1:49]  ##  wrong, must use CDS+UTR, table.1=mm.exons.keep.first.isoform.proteincoding[,1:49]
      org.gtf.input.wattr<-hg19.gtf.input.wattr
      mean.utr.3.size.allgenes=1363.202
      mean.utr.5.size.allgenes=254.5384
      mean.cds.size.allgenes=1638.583
    }
    if(input$org=="mm10"){
      mm10.CDS.UTR.keep.first.isoform<-mm10.gtf.input.wattr[which((mm10.gtf.input.wattr$X16 %in% as.vector(unique(mm.exons.keep.first.isoform.proteincoding$X16)))
                                                                  & (mm10.gtf.input.wattr$V3 == 'UTR' | mm10.gtf.input.wattr$V3 == 'CDS')), ]					
      table.1<-mm10.CDS.UTR.keep.first.isoform[,1:49]  ##  wrong, must use CDS+UTR, table.1=mm.exons.keep.first.isoform.proteincoding[,1:49]
      org.gtf.input.wattr<-mm10.gtf.input.wattr
      mean.utr.3.size.allgenes=1134.835
      mean.utr.5.size.allgenes=218.4317
      mean.cds.size.allgenes=1572.418
    }
    
    runs=c("f1")
    for(run in runs){
      table.1.2.merged=matrix("", dim(table.1)[1], dim(table.1)[2]+dim(table.2)[2])
      m=1
      for(j in 1:1:length(table(table.2$Chr))){
        name=names(table(table.2$Chr))[j]
        table2<-table.2[which(table.2$Chr==name),]
        table1<-table.1[which(table.1$V1==name),]
        table1<-table1[order(table1$V4),]
        for(i in 1:dim(table2)[1]){
          target=table2[i,4]
          if(target>=min(table1$V4)& target<=max(table1$V4)){
            index1=binary_search_test(target = target,table1$V4)
            table3=table1[1:(index1-1),]
            table3=table3[order(table3$V5),]
            index2=binary_search_test(target = target,table3$V5)
            if(is.integer(index2)){
              matched.ids=row.names(table3[index2:dim(table3)[1],])
              for(k in matched.ids){
                table.1.2.merged[m,]=c(as.character(as.matrix(table1[which(row.names(table1)==k), ])), as.character(as.matrix(table2[i,])))
                m=m+1
              }	
            }
          }
        }
      }
      table.1.2.merged=table.1.2.merged[1:(m-1),]
      table.1.2.merged=data.frame(table.1.2.merged, stringsAsFactors=FALSE)
      assign(paste("table.1.2.merged.peak.", run, ".", data.date, sep=""), table.1.2.merged)
    }
    
    
    
    if(input$analysis=="STOPCODON"){
      distance.to.stop.codon=c()
      for(i in 1:dim(table.1.2.merged.peak.f1.01212017)[1]){
        peak.summit=as.numeric(table.1.2.merged.peak.f1.01212017[i, 'X53'])
        transcript.name=table.1.2.merged.peak.f1.01212017[i, 'X25']
        row.index=which(org.gtf.input.wattr$V3=='stop_codon' & org.gtf.input.wattr$X16==transcript.name)[1]
        if(!is.na(row.index)){
          stop.codon = org.gtf.input.wattr[row.index,'V4']
          strand=org.gtf.input.wattr[row.index, 'V7']
          
          start.end.all.exons=org.gtf.input.wattr[which(org.gtf.input.wattr$V3=='exon'  & org.gtf.input.wattr$X16==transcript.name), 4:5] #table.1[which(exons$X16==isof), 4:5]
          start.end.order.index=order(start.end.all.exons[,1])
          start.end.all.exons=start.end.all.exons[start.end.order.index,]
          
          distance.to.stop.codon=c(distance.to.stop.codon, 
                                   find.peak.to.stopcodon.distance.by.strand(start.end.all.exons, stop.codon, peak.summit, strand))
        }
      }
      human.01212017.distance.to.stop.codon=distance.to.stop.codon
      human.01212017.distance.to.stop.codon.chosen=human.01212017.distance.to.stop.codon[which(human.01212017.distance.to.stop.codon>(-1000)&human.01212017.distance.to.stop.codon<1000)]
    }
    
    
    if(input$analysis=="LASTEXON"){
      
      ## distance from peak summit to last exon junction
      isof.norm.dist.matrix=matrix("", nrow(table.1.2.merged.peak.f1.01212017), 4)
      distance.to.last.exon=c()
      for(i in 1:dim(table.1.2.merged.peak.f1.01212017)[1]){
        
        summit.peak=as.numeric(table.1.2.merged.peak.f1.01212017[i, 'X53'])  # c(4,5,54)  for cirseq (because the peak list is diff))
        
        isof=table.1.2.merged.peak.f1.01212017[i,'X25']
        strand=table.1.2.merged.peak.f1.01212017[i,'X7']
        
        start.end.all.exons=org.gtf.input.wattr[which(org.gtf.input.wattr$V3=='exon'  & org.gtf.input.wattr$X16==isof), 4:5] #table.1[which(exons$X16==isof), 4:5]
        start.end.order.index=order(start.end.all.exons[,1])
        start.end.all.exons=start.end.all.exons[start.end.order.index,]
        
        last.exon=c()
        if(strand=='-'){
          last.exon=start.end.all.exons[1,]
        }else{
          last.exon=start.end.all.exons[dim(start.end.all.exons)[1],]			
        }
        
        last.exon.pos=if(summit.peak>last.exon[2]){last.exon[2]}else if(summit.peak<last.exon[1]){last.exon[1]}else{0}
        last.exon.pos=as.numeric(last.exon.pos)
        
        normalized.distance.to.pivot=find.peak.to.lastexon.distance.by.strand(start.end.all.exons, last.exon.pos, summit.peak, strand)
        distance.to.last.exon=c(distance.to.last.exon, normalized.distance.to.pivot)
        
      }
      human.01212017.distance.to.last.exon=distance.to.last.exon
      human.01212017.distance.to.last.exon.chosen=human.01212017.distance.to.last.exon[which(human.01212017.distance.to.last.exon>(-1000)&human.01212017.distance.to.last.exon<1000)]
      #human.01212017.distance.to.last.exon.chosen=human.01212017.distance.to.last.exon.chosen[which(human.01212017.distance.to.last.exon.chosen!=0)]
    }
    
    
    if(input$analysis=="53UTRCDS"){
      
      ## need to make sure column 53 (or 54)  
      for(run in runs){  
        isof.norm.dist.matrix=c()	
        table.1.2.merged=get(paste("table.1.2.merged.peak.", run, ".", data.date, sep=""))
        unique.peaks=c()
        #for(isof in gene.peak.isof){
        #	isof.row.indices=which(table.1.2.merged$X25==isof)
        #	for(index in isof.row.indices){
        for(i in 1:dim(table.1.2.merged)[1]){
          #isof.info=table.1.2.merged[index, c(4,5,54)]  # c(4,5,54)  for cirseq (because the peak list is diff)
          index=i
          #isof.info=table.1.2.merged[index, c(4,5,ncol(table.1.2.merged))]  # c(4,5,54)  for cirseq (because the peak list is diff)
          isof.info=table.1.2.merged[index, c(4,5,53)]  # c(4,5,54)  for cirseq (because the peak list is diff)
          
          isof=table.1.2.merged[index,'X25']
          pivot=as.numeric(isof.info[3])
          this.type=table.1.2.merged[index, 3]
          transcript.type=table.1.2.merged[index, 'X21']
          unique.peak=paste(c(table.1.2.merged[index, c('X1' , 'X19')], pivot), collapse="-")  # c('X1' , 'X4', 'X5')
          if(length(which(unique.peaks==unique.peak))==0){		
            start.end.all=table.1[which(table.1$X16==isof), 4:5]
            start.end.order.index=order(start.end.all[,1])
            start.end.all=start.end.all[start.end.order.index,]
            types=table.1[which(table.1$X16==isof), 3]
            types=types[start.end.order.index]
            
            #start.end=table.1[which(table.1$X16==isof & table.1$V3==type), 4:5]
            #start.end=start.end[order(start.end[,1]),]
            strand=unique(table.1[which(table.1$X16==isof), 'V7'])
            row.vector=data.matrix(isof.info[1:2])
            UTR.5=1
            start.end=c()
            if(this.type=="UTR"){
              side.and.utr.index=get.utr.side(start.end.all, row.vector, types, isof)
              UTR.5=side.and.utr.index[[1]]
              utr.index=side.and.utr.index[[2]]
              start.end=start.end.all[utr.index,]
            }else{
              start.end=start.end.all[which(types=='CDS'),]
            }
            data.matrix=start.end
            if(nrow(data.matrix)>0){
              normalized.distance.to.pivot=find.normalized.position(data.matrix, row.vector, pivot)	
              if(strand=="-"){
                UTR.5=1-UTR.5
                normalized.distance.to.pivot=100-normalized.distance.to.pivot
              }
              isof.norm.dist.matrix=rbind(isof.norm.dist.matrix, c(isof, this.type, UTR.5, normalized.distance.to.pivot, transcript.type))
            }
          }
          unique.peaks=c(unique.peaks, unique.peak)
        }
        name=paste("isof.norm.dist.matrix.unique", run, ".", data.date, sep=".")
        assign(name, isof.norm.dist.matrix)
      }
      
      #name=paste("isof.norm.dist.matrix.", run, sep=".")
      #name=paste("isof.norm.dist.matrix.unique.", run, sep=".")  #OLD VERSION. TO BE DELETED.
      name=paste("isof.norm.dist.matrix.unique", run, ".", data.date, sep=".")
      isof.norm.dist.matrix=get(name)
      
      UTR.5.list=as.numeric(as.character(as.matrix(isof.norm.dist.matrix[which(isof.norm.dist.matrix[,5]=="protein_coding;" & isof.norm.dist.matrix[,2]=="UTR" & isof.norm.dist.matrix[,3]=="1"), 4])))
      UTR.3.list=as.numeric(as.character(as.matrix(isof.norm.dist.matrix[which(isof.norm.dist.matrix[,5]=="protein_coding;" & isof.norm.dist.matrix[,2]=="UTR" & isof.norm.dist.matrix[,3]=="0"), 4])))
      CDS.list  =as.numeric(as.character(as.matrix(isof.norm.dist.matrix[which(isof.norm.dist.matrix[,5]=="protein_coding;" & isof.norm.dist.matrix[,2]=="CDS"), 4])))
      
      peak.normalized.list=c(UTR.5.list*mean.utr.5.size.allgenes/100, 
                             mean.utr.5.size.allgenes+CDS.list*mean.cds.size.allgenes/100, 
                             mean.utr.5.size.allgenes+mean.cds.size.allgenes+UTR.3.list*mean.utr.3.size.allgenes/100)
      #hist(peak.normalized.list, breaks=20, prob=TRUE, col="grey")
      #lines(density(peak.normalized.list, adjust=2), type="l", lty="dotted", col="darkgreen", lwd=2) 
      
      
      dvals=density(peak.normalized.list, adjust=0.5)
      x=dvals$x[which(dvals$x>=0 & dvals$x<=3300)]
      y=dvals$y[which(dvals$x>=0 & dvals$x<=3300)]
      x=c(0,x)
      y=c(0,y)
      
    }
    
    
    
    if(input$analysis=="STOPCODON"){
      plotdata<-human.01212017.distance.to.stop.codon.chosen
    }
    else if(input$analysis=="LASTEXON"){
      plotdata<-human.01212017.distance.to.last.exon.chosen
    }
    else if(input$analysis=="53UTRCDS"){
      plotdata<-peak.normalized.list
    }
    list(plotdata=plotdata,mydata=table.1.2.merged)
    }
  })
  
  
  
  output$rect2 <- renderPlot({
    dat<-dataInput2()$mydata 
    sub<-dat[which(as.character(strsplit(dat$X19,";"))==input$gene),c(19,51:52)]
    sub$X51<-as.numeric(sub$X51)
    sub$X52<-as.numeric(sub$X52)
    
    plot(c(min(sub$X51)-50,max(sub$X52)+50),c(0,2),type = "n",xlab = "location",ylab = "",yaxt='n')
    y<-seq(1,10,by=1)
    for(i in 1:dim(sub)[1]){
      rect(sub$X51[i],0,sub$X52[i],0.5,border=rainbow(10)[i],lty = y[i],lwd = y[i])
    }
    
  })
  options(DT.options = list(pageLength = 5))
  output$x1 = DT::renderDataTable(datasetInput(), server = FALSE)
  
  datasetInput <- reactive({
    switch(input$dataset,
           "merged_case" = as.matrix(dataInput()$mydata)[,-9],
           "peak_case" =as.matrix(dataInput()$plotdata),
             "merged_con"=as.matrix(dataInput2()$mydata)[,-9],            
           "peak_con"=as.matrix(dataInput2()$plotdata))
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$dataset, '.csv', sep='')},
    content = function(file) {
      write.csv(as.data.frame(datasetInput()), sep="\t", file)
    }
  )
  
  

  
  
  output$distPlot <- renderPlot({
    if(input$type=="density"){
      if(input$analysis!="53UTRCDS"){
        plot(density(dataInput()$plotdata, bw=50), col=input$color,xlim = c(input$range) ,ylim=c(input$rangey), main = paste(input$org,input$analysis,"density plot",sep = " "),xlab =paste("the distance in base pairs between the peak and the\n",input$analysis,
                                                                                                                                                                                   "(negative for upstream peak and positive for downstream peak)",sep = " "))
        if(!is.null(dataInput2()$plotdata)){
        lines(density(dataInput2()$plotdata,bw=50),col=input$color2,xlim=c(input$range),ylim=c(input$rangey)) 
        }
      }
      else{
        if(input$org=="hg19"){
          mean.utr.3.size.allgenes=1363.202
          mean.utr.5.size.allgenes=254.5384
          mean.cds.size.allgenes=1638.583
          plot(density(dataInput()$plotdata, bw=50), xaxt="n",col=input$color,xlim = c(0,3200) ,ylim=c(input$rangey), main = paste(input$org,input$analysis,"density plot",sep = " "),xlab = "The normalized distribtion of the peaks in any given gene \n (5'UTR, CDS and 3'UTR in each gene are normalized to be the averages of their values using all genes)")
          if(!is.null(dataInput2()$plotdata)){
          lines(density(dataInput2()$plotdata,bw=50),col=input$color2,xlim=c(0,3200))
          }
            axis(1, c(mean.utr.5.size.allgenes, mean.utr.5.size.allgenes+mean.cds.size.allgenes, mean.utr.5.size.allgenes+mean.cds.size.allgenes+mean.utr.3.size.allgenes),
               labels=c("5-UTR", "CDS", "3-UTR")) 
        }
        else if (input$org=="mm10"){
          mean.utr.3.size.allgenes=1134.835
          mean.utr.5.size.allgenes=218.4317
          mean.cds.size.allgenes=1572.418
          plot(density(dataInput()$plotdata, bw=50), xaxt="n",col=input$color,xlim = c(0,3200) , ylim=c(input$rangey),main = paste(input$org,input$analysis,"density plot",sep = " "),xlab = "The normalized distribtion of the peaks in any given gene \n (5'UTR, CDS and 3'UTR in each gene are normalized to be the averages of their values using all genes)")
          if(!is.null(dataInput2()$plotdata)){
          lines(density(dataInput2()$plotdata,bw=50),col=input$color2,xlim=c(0,3200))
          }
            axis(1, c(mean.utr.5.size.allgenes, mean.utr.5.size.allgenes+mean.cds.size.allgenes, mean.utr.5.size.allgenes+mean.cds.size.allgenes+mean.utr.3.size.allgenes),
               labels=c("5-UTR", "CDS", "3-UTR"))
        }
      }
    }
    else if (input$type=="histogram(only for case group)"){
      if(input$analysis!="53UTRCDS"){
        hist(dataInput()$plotdata,col=input$color,xlim = c(input$range) , main = paste(input$org,input$analysis,"density plot",sep = " "),xlab =paste("the distance in base pairs between the peak and the\n",input$analysis,
                                                                                                                                                      "(negative for upstream peak and positive for downstream peak)",sep = " "))
      }
      else{
        hist(dataInput()$plotdata, breaks=20, prob=TRUE, col=input$color,xlab = "The normalized distribtion of the peaks in any given gene \n (5'UTR, CDS and 3'UTR in each gene are normalized to be the averages of their values using all genes)")
      }
    }
  })
  
  
})