memory.size(max = TRUE)
memory.limit(size = NA)

library(parallel)
library(doParallel)
library(MASS)

numCores<-detectCores()
numCores
options(cores=8)
registerDoParallel(8)
getDoParWorkers()

library(latticeExtra)
library(TCGAbiolinks)
library(TCGAbiolinksGUI.data)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library("edgeR")
library(sjstats)
library(dplyr)


lgg.gbm.subtype<- TCGAquery_subtype(tumor = "lgg")
query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts",
                  barcode = lgg.gbm.subtype$patient)
GDCdownload(query)
data <- GDCprepare(query)
GenewiseCounts <- data.frame(data@assays@data@listData[["HTSeq - Counts"]],
                             row.names = data@rowRanges@ranges@NAMES)

colnames(GenewiseCounts)<-data@colData@rownames
colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,12)
GenewiseCounts<-cbind(data@rowRanges@ranges@width,GenewiseCounts)
colnames(GenewiseCounts)[1]<-paste("Length")
y <- DGEList(GenewiseCounts[,-1], group=NULL, genes=GenewiseCounts[,1,drop=FALSE])
options(digits=3)
y <- y[!is.na(y$genes$Length), ]
dim(y)
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y), keytype="ENTREZID", column="SYMBOL")
keep <- rowSums(cpm(y) > 0.5) >= 2
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
options(digits=3)
y$samples

num<-as.matrix(1:56457)
num<- num^2/20

row.names(num)<-rownames (GenewiseCounts)
y2=y$counts
logCMP<- cpm(y2, prior.count=2, log=F)
logCMP<-as.data.frame(logCMP)
num<-as.data.frame(num)
NewW<-merge(logCMP,num, by="row.names")
row.names(NewW)<-NewW$Row.names

NumDat<- cbind(NewW$V1)
row.names(NumDat)<-NewW$Row.names 
NewW$Row.names<-NULL
NewW$V1<-NULL
Datmulnumnum<- NumDat*NewW

#require(openxlsx)
#list_of_datasets <- list("Datmulnumnum"=Datmul)
#write.xlsx(list_of_datasets, file = "Datmulbuen.xlsx")

#clear ram 
rm(Datmul)
rm(GenewiseCounts)


#tablas 1 y 2 
tabla1<- data.frame(ENSG=row.names(tabla1),tabla1)
tabla2<- data.frame(ENSG=row.names(tabla2),tabla2)


#crear la nueva tabla de contingencia de las muestras 
targets<-data.frame(data@colData@listData[["paper_Grade"]],
                    data@colData@listData[["paper_Histology"]],
                    stringsAsFactors = F)


row.names(targets)<- data@colData@rownames
colnames(targets)<-c("Grade","Histology")

#Crear la tabla de grupos 
group <- paste(targets$Grade,targets$Histology,sep="/") 
group<-as.data.frame(group)

#diplr tablas de conteo por count 
# tabla de conteos y fracciones 
labels<-read.csv("labels5.csv")
graphic<- cbind(labels, group)
graphic$X<-NULL
colnames(graphic)<- c("labels","group")

tabla<-table(graphic$labels)
tabla<-as.data.frame(tabla) 

graphic2 <- paste(graphic$labels,graphic$group,sep="-") 
graphic2<-as.data.frame(graphic2)


tgraph1<-table(group$group)
tgraph1<-as.data.frame(tgraph1) 

tgraph2<-table(graphic2)
tgraph2<-as.data.frame(tgraph2)


################## aqui empiezan las grafias############3
library(ggplot2)

ggplot(graphic, aes(labels))+
  geom_histogram()

ggplot(tabla, aes(Var1,Freq, label= Freq, fill =Var1))+
  geom_bar(stat = "identity")+
  geom_text(size=4,vjust=2,color="white")


ggplot(tgraph1, aes(Var1,Freq, label= Freq, fill =Var1))+
  geom_bar(stat = "identity")+
  geom_text(size=4,vjust=2,color="white")


ggplot(tgraph2, aes(graphic2,Freq, label= Freq, fill =graphic2))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))+
  geom_text(size=4,vjust=2,color="white")

Q <- paste(graphic$labels,graphic$group,sep="-") 
#Q<-cbind(Q,graphic)
Q<-table(Q)
Q<-as.data.frame(Q)
E<-levels(tgraph1$Var1)
as.data.frame(E)
E<-c(E,E)
#View(E)
#E<-E[-c(6,16,32)]
#labels2E<-E[-c(22,30)]
W<-cbind(E,Q)

theme_Publication <- function(base_size=14, base_family="Times New Roman") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(text = element_text(),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "black", size = rel(2)),
            axis.text = element_text(size = rel(2)),
            axis.ticks = element_line(size = rel(2)),
            axis.ticks.length = unit(-0.25, "cm"),
            axis.text.x = element_text(margin = margin(t = .8, unit = "cm", l = 50),
                                       angle = -55, size = 12, hjust = 0),
            axis.text.y = element_text(margin = margin(r = .8, unit = "cm")),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour="#CCCCCC", size = rel(1)),
            panel.grid.minor.x = element_line(colour="#CCCCCC", linetype = "dashed", size = rel(1))
    ))
  
}


ggplot(data = W , aes(x = Q, y =Freq , fill = E)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(breaks=seq(0,1900,100))+
  scale_fill_discrete(name = "",
                      labels = c("G2/astrocytoma",
                                 "G2/oligoastrocytoma","G2/oligodendroglioma"
                                 ,"G3/astrocytoma	","G3/oligoastrocytoma",
                                 "G3/oligodendroglioma","NA/NA")) +
  theme_Publication() 


#----------pocentaje

X0<-W[c(1:7),]
X1<-W[c(8:14),]

X1$percentage<-(round(X0$Freq/sum(X0$Freq),3)*100)
X0$percentage<-(round(X1$Freq/sum(X1$Freq),3)*100)
row.names(X0)<-X0$E



pcent <- as.data.frame(cbind(row.names(X0),
                             X0$percentage, X1$percentage))

colnames(pcent) <- c("Tipo","0","1")

pcent2 <- tidyr::gather(pcent, key = type, value = value, -Tipo)
pcent2$value <- as.numeric(pcent2$value)

pcent2 <- pcent2[order(pcent2$value, decreasing = T),]
#pcent2<-pcent2[-c(6,16,32),]
str(pcent2)

#grafica por tipo
ggplot(data = pcent2 , aes(x =type, y =value,fill = Tipo)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=sprintf("%0.2f",round(value,digits = 2))),
            position=position_dodge(width=1),
            size=3, angle=90,
            hjust=1,vjust=.4,color="white")+
  scale_y_continuous(breaks=seq(0,100,10))+
  scale_fill_discrete(name = "", labels = c("G2/astrocytoma",
                                            "G2/oligoastrocytoma","G2/oligodendroglioma"
                                            ,"G3/astrocytoma	","G3/oligoastrocytoma",
                                            "G3/oligodendroglioma","NA/NA")) +
  theme_Publication()


#grafica por numero 

#grafia de solo oligo y astro 

grouphist<- paste(targets$Histology) 
grouphist<-as.data.frame(grouphist)
graphichist<- cbind(labels, grouphist)
graphichist$X<-NULL
colnames(graphichist)<- c("labels","grouphist")

graphichist<-table(graphichist)
graphichist<-as.data.frame(graphichist)

#rescue each label that has 1 and 0 and create an especific DF
#hacer un combine c()
hist0<-graphichist[graphichist[,1]=="0",]
hist0<-as.data.frame(hist0)
hist1<-graphichist[graphichist[,1]=="1",]
hist1<-as.data.frame(hist1)


hist0$percentage<-(round(hist0$Freq/sum(hist0$Freq),2)*100)
hist1$percentage<-(round(hist1$Freq/sum(hist1$Freq),2)*100)
row.names(hist0)<-hist0$grouphist
pcenthist <- as.data.frame(cbind(row.names(hist0),
                                 hist0$percentage, hist1$percentage))

colnames(pcenthist) <- c("histology","0","1")
pcenthist2 <- tidyr::gather(pcenthist, key = type,
                            value = value, -histology)

pcenthist2$value <- as.numeric(pcenthist2$value)
pcenthist2 <- pcenthist2[order(pcenthist2$value, decreasing = T),]
str(pcent2)

ggplot(data = pcenthist2 , aes(x =type, y =value,fill = histology)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=sprintf("%0.2f",round(value,digits = 2))),
            position=position_dodge(width=1),
            size=3, angle=90,
            hjust=1,vjust=.4,color="white")+
  scale_y_continuous(breaks=seq(0,100,10))+
  scale_fill_discrete(name = "", labels = c("astrocytoma",
                                            "NA/NA",
                                            "oligoastrocytoma",
                                            "oligodendroglioma")) +
  theme_Publication()




hist0notna<-hist0[-c(2),]
hist1notna<-hist1[-c(2),]
hist0notna$percentage<-(round(hist0notna$Freq/sum(hist0notna$Freq),3)*100)
hist1notna$percentage<-(round(hist1notna$Freq/sum(hist1notna$Freq),3)*100)

pcenthistnotna <- as.data.frame(cbind(row.names(hist0notna),
                                      hist0notna$percentage,
                                      hist1notna$percentage))


colnames(pcenthistnotna) <- c("histology","0","1")
pcenthist2notna <- tidyr::gather(pcenthistnotna, key = type,
                                 value = value, -histology)

pcenthist2notna$value <- as.numeric(pcenthist2notna$value)
pcenthist2notna <- pcenthist2notna[order(pcenthist2notna$value,
                                         decreasing = T),]
str(pcent2notna)

ggplot(data = pcenthist2notna , aes(x =type, y =value,fill = histology)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=sprintf("%0.2f",round(value,digits = 0))),
            position=position_dodge(width=1),
            size=3.5, angle=90,
            hjust=1,vjust=.4,color="white")+
  scale_y_continuous(breaks=seq(0,100,10))+
  scale_fill_discrete(name = "", labels = c("astrocytoma",
                                            "oligoastrocytoma",
                                            "oligodendroglioma")) +
  labs(x="Grupo",title = "% grupos post analisis")+
  theme_Publication()

tgraph1hist<-table(grouphist$grouphist)
tgraph1hist<-as.data.frame(tgraph1hist) 

ggplot(tgraph1hist, aes(Var1,Freq, label= Freq, fill =Var1))+
  geom_bar(stat = "identity")+
  geom_text(size=6,vjust=2,color="white")+
  labs(x="Grupo",title = "Frecuencias grupos sin analisis")+
  theme_Publication() 

tgraph1hist$percentage<-(round(tgraph1hist$Freq/sum(tgraph1hist$Freq),2)*100)

ggplot(tgraph1hist, aes(Var1,percentage, label= percentage,
                        fill =Var1))+
  geom_bar(stat = "identity")+
  geom_text(size=6,vjust=2,color="white")+
  labs(x="Grupo",title = "% grupos sin analisis")+
  theme_Publication() 

tabla$percentage<-(round(tabla$Freq/sum(tabla$Freq),2)*100)

ggplot(tabla, aes(Var1,percentage, label= percentage,
                  fill= Var1))+
  geom_bar(stat = "identity")+
  geom_text(size=8,vjust=2,color="white")+
  labs(x="Grupo",title = "% grupos con analisis")+
  scale_fill_discrete(name = "", labels = c("0-oligodendroglioma",
                                            "1-astrocytoma")) +
  theme_Publication() 

#tablas 
tabladeconteos<- y$counts

grouphist<- paste(targets$Histology) 















hist0notna<-hist0[-c(2),]
hist1notna<-hist1[-c(2),]
hist0notna$percentage<-(round(hist0notna$Freq/sum(hist0notna$Freq),3)*100)
hist1notna$percentage<-(round(hist1notna$Freq/sum(hist1notna$Freq),3)*100)

pcenthistnotna <- as.data.frame(cbind(row.names(hist0notna),
                                      hist0notna$percentage,
                                      hist1notna$percentage))


colnames(pcenthistnotna) <- c("histology","0","1")
pcenthist2notna <- tidyr::gather(pcenthistnotna, key = type,
                                 value = value, -histology)

pcenthist2notna$value <- as.numeric(pcenthist2notna$value)
pcenthist2notna <- pcenthist2notna[order(pcenthist2notna$value,
                                         decreasing = T),]
str(pcent2notna)

ggplot(data = pcenthist2notna , aes(x =type, y =value,fill = histology)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=sprintf("%0.2f",round(value,digits = 0))),
            position=position_dodge(width=1),
            size=3.5, angle=90,
            hjust=1,vjust=.4,color="white")+
  scale_y_continuous(breaks=seq(0,100,10))+
  scale_fill_discrete(name = "", labels = c("astrocytoma",
                                            "oligodendroglioma")) +
  labs(x="Grupo",title = "% grupos post analisis")+
  theme_Publication()




hist0<-graphichist[graphichist[,1]=="0",]
hist0<-as.data.frame(hist0)
hist1<-graphichist[graphichist[,1]=="1",]
hist1<-as.data.frame(hist1)


hist0$percentage<-(round(hist0$Freq/sum(hist0$Freq),2)*100)
hist1$percentage<-(round(hist1$Freq/sum(hist1$Freq),2)*100)
row.names(hist0)<-hist0$grouphist
pcenthist <- as.data.frame(cbind(row.names(hist0),
                                 hist0$percentage, hist1$percentage))

colnames(pcenthist) <- c("histology","0","1")
pcenthist2 <- tidyr::gather(pcenthist, key = type,
                            value = value, -histology)

pcenthist2$value <- as.numeric(pcenthist2$value)
pcenthist2 <- pcenthist2[order(pcenthist2$value, decreasing = T),]
str(pcent2)

ggplot(data = pcenthist2 , aes(x =type, y =value,fill = histology)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=sprintf("%0.2f",round(value,digits = 2))),
            position=position_dodge(width=1),
            size=3, angle=90,
            hjust=1,vjust=.4,color="white")+
  scale_y_continuous(breaks=seq(0,100,10))+
  scale_fill_discrete(name = "", labels = c("astrocytoma",
                                            "NA/NA",
                                            "oligoastrocytoma",
                                            "oligodendroglioma")) +
  theme_Publication()
