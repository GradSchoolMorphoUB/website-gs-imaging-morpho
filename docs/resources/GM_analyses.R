#load libraries
library(scatterplot3d)
library(shapes)
library(ade4)
library(Morpho)
library(rgl)
library(MASS)
library(pairwiseAdonis)
library(vegan)

#load data set (semilandmarks.csv)
data<-read.table("Semilandmarks.csv",header=FALSE,sep=";")

#define some parameters
#number of dimensions
dimensions=3

#a group with genera names
name=as.factor(data[,1])

#number of individuals
nbIndividu=nrow(data)

#number of landmarks
nbLandmark=(ncol(data)-1)/3

#selecting only semilandmark data (i.e., removing the first column that include genera names)
data2=as.matrix(data[,-1])

#reorganizing the data as an array
data2<-aperm(array(t(data2),dim=c(dimensions,nbLandmark,nbIndividu)),c(2,1,3))

####running the Procrustes superimposition
gp1<-procGPA(data2,scale=TRUE)

####PCA plots
#isolating percentage of variance
x=round(gp1$percent, digits=2)
x1=as.factor(paste("PC1 (",x[1],"%)")) 
x2=as.factor(paste("PC2 (",x[2],"%)"))

#PC1 and PC2 plot
plot(gp1$scores[,c(1,2)],xlab=x1,ylab=x2,asp=1,pch=as.numeric(name),col=as.numeric(name))

#vertical and horizontal dotted lines passing by zero
abline(v=0,h=0,lwd=.5,lty=3)

#adding labels
text(gp1$scores[,c(1,2)],labels=name,pos=1,offset=0.1,cex=0.5)

#creating convex hulls for groups
s.chull (gp1$scores[,c(1,2)],fac=name,xax=1,optchull=1,clabel=0,add.plot=TRUE)

#now, let's look at PC2 vs. PC3
x3=as.factor(paste("PC3 (",x[3],"%)"))
plot(gp1$scores[,c(2,3)],xlab=x2,ylab=x3,asp=1,pch=as.numeric(name),col=as.numeric(name))
abline(v=0,h=0,lwd=.5,lty=3)
text(gp1$scores[,c(2,3)],labels=name,pos=1,offset=0.1,cex=0.5)
s.chull (gp1$scores[,c(2,3)],fac=name,xax=1,optchull=1,clabel=0,add.plot=TRUE)

#you can also explore other PC combinations and further PCs

####extreme shapes PC1
source(file.choose()) # open functions.R

#creating an object with rotated data for visualization
procrustes<-array2mat(gp1$rotated)

#computing extreme shapes along PC1
#minimum shape
minPC1<-extshapes(procrustes,"pca",gp1,1,"min")
array_minPC1<- aperm(array(minPC1,c(3,length(minPC1)/3)),perm=c(2,1))

#visualize in 3D the minimum shape
points3d(array_minPC1,size=5,col="blue")

#maximal shape
maxPC1<-extshapes(procrustes,"pca",gp1,1,"max")
array_maxPC1<- aperm(array(maxPC1,c(3,length(maxPC1)/3)),perm=c(2,1))

#visualize in 3D the maximum shape
points3d(array_maxPC1,size=5,col="red")

#change to ortogonal view (by default it is in perspective view)
par3d(FOV=0)

#save pdf with extreme shape points (as visible on screen)
#rgl.postscript("PC1_min-max.pdf","pdf")

#close the rgl window
rgl.close()


####between group PCA
#computing again Procrustes superimposition (with Morpho function procSym)
gp2=procSym(data2)
group=name

#computing the between-group PCA
bgPCA=groupPCA(dataarray=gp2$orpdata,groups=group,rounds=1000,cv=T,weighting=T)

#calculating the percentage of variance along bgPC1 and bgPC2
X1=as.factor(paste("bgPC1 (",round(bgPCA$Variance$exVar[1]*100,digits=2),"%)",sep=""))
X2=as.factor(paste("bgPC2 (",round(bgPCA$Variance$exVar[2]*100,digits=2),"%)",sep=""))

#normal bgPC1 and bgPC2 plot
plot(bgPCA$Scores[,1:2],xlab=X1,ylab=X2,pch=as.numeric(group),col=as.numeric(group),asp=1)

#vertical and horizontal dotted lines passing by zero
abline(v=0,h=0,lwd=.5,lty=3)

#adding labels
text(bgPCA$Scores[,1:2],labels=name,pos=1,offset=0.1,cex=0.5)

#creating convex hulls for groups
s.chull (bgPCA$Scores[,1:2],fac=group,xax=1,optchull=1,clabel=0,add.plot=TRUE)

#now try to plot bgPC2 and bgPC3

#cross-validated bgPCA plot
plot(bgPCA$CV[,1:2],xlab="cv-bgPC1",ylab="cv-bgPC2",pch=as.numeric(group),col=as.numeric(group),asp=1)
abline(v=0,h=0,lwd=.5,lty=3)
text(bgPCA$CV[,1:2],labels=name,pos=1,offset=0.1,cex=0.5)
s.chull (bgPCA$CV[,1:2],fac=group,xax=1,optchull=1,clabel=0,add.plot=TRUE)

#checking the degree of correct classification
bgPCA


####extreme shapes bgPC1
#min bgPC1
minbgPC1<-extshapes(procrustes,"bgpca",bgPCA,1,"min")
array_minbgPC1<- aperm(array(minbgPC1,c(3,length(minbgPC1)/3)),perm=c(2,1))

#visualize in 3D the minimum shape
points3d(array_minbgPC1,size=5,col="blue")

#max bgPC1
maxbgPC1<-extshapes(procrustes,"bgpca",bgPCA,1,"max")
array_maxbgPC1<- aperm(array(maxbgPC1,c(3,length(maxbgPC1)/3)),perm=c(2,1))

#visualize in 3D the maximum shape
points3d(array_maxbgPC1,size=5,col="red")

#change to ortogonal view (by default it is in perspective view)
par3d(FOV=0)

#save pdf with extreme shape points (as visible on screen)
rgl.postscript("bgPC1_min-max.pdf","pdf")

#close the rgl window
rgl.close()


###CVA with Morpho
# input data
ncomp=40
dat <- as.data.frame(cbind(Group = group, gp1$scores, 1:ncomp))

# plot the highest cross-classification accuracy 
accuracies <- rep(NA, ncomp - 1)
names(accuracies) <- paste(2:ncomp, "axes")
for (k in 3:ncol(dat)) {
    mod <- lda(dat$Group ~ ., data = dat[, 1:k],
               prior = rep(1/4, 4),
               CV = TRUE)
    accuracies[k - 2] <- sum(mod$class == dat$Group) / nrow(dat)
}
plot(x = 2:ncomp, y = accuracies,
     type = "b", pch = 16,
     xlab = "Number of PCA axes in LDA model",
     ylab = "X-val classification accuracy")
grid()

#display the cumulative percentage of variation along PC scores 
cumsum(gp1$percent[1:40])

#create data objects based on the number of PC scores relevant for CVA
resDonnees2=gp1$scores[,1:8]

#run CVA
cva=CVA(resDonnees2, group, cv=TRUE, prior=rep(1/4,4), weighting=TRUE)

#isolating percentage of variance for CV1 and CV2
percx=round(cva$Var[1,2],digits=2)
percy=round(cva$Var[2,2],digits=2)

#plot CV1 and CV2
plot(cva$CVscores,xlab=paste("CV1 (",percx,"%)"),ylab=paste("CV2 (",percy,"%)"),asp=1,cex=1,col=as.numeric(group),pch=as.numeric(group))

#vertical and horizontal dotted lines passing by zero
abline(v=0,h=0,lwd=1,lty=3)

#adding labels
text(cva$CVscores[,c(1,2)],labels=name,pos=1,offset=0.1,cex=0.5)

#creating convex hulls for groups
s.chull(cva$CVscores,fac=group,xax=1,optchull=1,clabel=0,add.plot=TRUE)

#check classification accuracy
classify(cva)

####extreme shapes CVA
#min CV1
minCV1<-extshapes(procrustes,"lda",cva,1,"min")
array_minCV1<- aperm(array(minCV1,c(3,length(minCV1)/3)),perm=c(2,1))

#visualize in 3D the minimum shape
points3d(array_minCV1,size=5,col="blue")

#max CV1
maxCV1<-extshapes(procrustes,"lda",cva,1,"max")
array_maxCV1<- aperm(array(maxCV1,c(3,length(maxCV1)/3)),perm=c(2,1))

#visualize in 3D the maximum shape
points3d(array_maxCV1,size=5,col="red")

#change to ortogonal view (by default it is in perspective view)
par3d(FOV=0)

#save pdf with extreme shape points (as visible on screen)
rgl.postscript("CV1_min-max.pdf","pdf")

###allometry and size effect

#define object for centroid size
centroid.size=gp1$size

#plot PC1 and centroid size
plot(gp1$scores[,1],centroid.size,col=as.numeric(group))
abline(coef(lm(centroid.size~gp1$scores[,1])),lwd=2,col="red")
text(gp1$scores[,1],centroid.size,labels=name,pos=1,offset=0.1,cex=0.5)

#try to investigate other PC scores, as well as bgPCs and CVs

#statistically testing for allometry
perm1=adonis(formula=gp1$scores~centroid.size*group,permutations=9999,method="euc")
perm1