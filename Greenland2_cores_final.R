#3-2-2022
library(rioja)
library(analogue)
library(vegan)

setwd("C:\\Users\\Dal User\\Desktop\\northwestern\\greenland\\final")
species3<-read.csv(file="Greenland2021_conservative2.csv", row.names=1) #read the csv file

#load file that shows different regions
Kterms<-read.csv(file="Kterms.csv", row.names=1)
Cat <- Kterms$Cat
species3$Cat <- paste(Kterms$Cat)

#calculate the row sums because data is counts - check for low counts
species3$sums <-paste(rowSums(species3[1:124]))
species3$sums <- as.numeric(species3$sums)

species7<-na.omit(species3)

speciesX <- cbind(name = rownames(species7), species7)
rownames(speciesX) <- NULL

library(dplyr)
df<-speciesX %>% filter(sums >= 40) #use a threshold of 40, higharctic sites

library(tibble)
df<-column_to_rownames(df, var = "name")

Kterms<-df[, -cbind(1:124)]

#here I calculate relative abundance before removing all of the weird taxa
species8<-df[, -cbind(124:128) ] #removes Kterms

species9 <- species8 / rowSums(species8) * 100

species1<-species9[, -cbind(118:128) ] #this removes the taxa that are not good for modelling, undiff tanytarsini and such

write.csv(Kterms, "Kterms2.csv") #double check things

env<-read.csv(file="env2021b.csv", row.names=1) #read the envcsv file
rows_to_keep <- intersect(rownames(species1),rownames(env)) #ensure sites deleted from species list are deleted from env

env1 <- env[rows_to_keep, ,drop=FALSE]


species<-na.omit(species1)

tspecies <- decostand(species, method="hellinger") #transform data

#dca (b-diversity)
vare.dca <- decorana(tspecies)
vare.dca
summary(vare.dca)
plot(vare.dca)
plot(vare.dca, type = "n", xlim=c(-3,3), ylim=c(3,-3))
points(vare.dca, display = "sites", col = "black", cex = 1, pch = 21, bg = "red")
#text(vare.dca, display="species", pos=2)

Kterms2<-read.csv(file="Kterms2.csv", row.names=1, stringsAsFactors=T)
Cat <- Kterms2$Cat
with(Kterms2, levels(Cat))
library(ggplot2)
library(RColorBrewer)
myColors <- brewer.pal(7,"Set1")
myColors3 <- c("black", "red", "green", "blue", "#66FFFF")

names(myColors) <- levels(Kterms2$Cat)

dune.dist2 <- vegdist(tspecies)

mod4 <- betadisper(dune.dist2, Kterms2$Cat, type = "centroid")
mod4
permutest(mod4, permutations = 999)
anova(mod4)
#tiff("figS1_b.tiff", width = 6, height = 6, units = 'in', res = 150)
plot(mod4, ellipse = TRUE, hull = FALSE, conf = 0.95)
#dev.off()
boxplot(mod4)
plot(TukeyHSD(mod4))

#figure2b
#tiff("fig2_b.tiff", width = 6, height = 6, units = 'in', res = 150)
plot(vare.dca, type = "n", xlim=c(-3,3), ylim=c(-3,3))
points(vare.dca, display = "sites", col = "black",
       scaling = scl, cex = 1.5, pch = 21, bg = myColors3[Cat])
with(Kterms2, legend("bottomleft", legend = levels(Cat), bty = "n",
                     col = "black", pch = 21, pt.bg = myColors3))

#text(vare.dca, display="sites", pos=2)

mod.st <- scores(vare.dca, scaling=2, display="sites", choices=1:2)
#clusterplot
clus <- kmeans(mod.st,centers= 3, iter.max=1000, nstart=10000)
pl <- ordiellipse(vare.dca, scaling=2, clus$cluster, kind="sd", conf=0.95, lwd=2, col="gray", lty="dotted")


which(rownames(species) == "BLACK") 
#write.csv(species, "species.csv")

rows_to_keep <- intersect(rownames(species),rownames(Kterms2))

Kterms3 <- Kterms2[rows_to_keep, ,drop=FALSE]

Cat <- Kterms3$Cat
with(Kterms3, levels(Cat))
library(ggplot2)
library(RColorBrewer)
myColors2 <- brewer.pal(5,"Set1")
names(myColors) <- levels(Kterms3$Cat)

library(rioja)
library(ggpalaeo)
require("analogue")
require("ggplot2")
library(tibble)
 

spp1b<- species[,-(which(colSums(species)==0))]
spp1<-spp1b

rows_to_keep <- intersect(rownames(species),rownames(env1))

env5 <- env1[rows_to_keep, ,drop=FALSE]

spp <- decostand(spp1, method="hellinger")

mod <- MAT(spp, env5$July, lean=FALSE)
fit.md<-rioja:::crossval(mod, cv.method="boot", nboot=100)
autoplot(mod)
fit.md

#
fit <- WAPLS(spp, env5$WC_SUMMER, nlps=2)
fit
# cross-validate model
#rioja:::crossval(mod)
fit.cv <- rioja:::crossval(fit, cv.method="boot", nboot=1000)
# How many components to use?
rand.t.test(fit.cv)
screeplot(fit.cv)
fit.cv

autoplot(fit.cv, npls=2)

fit <- WA(spp, env5$WC_SUMMER, tolDW=TRUE)
# plot predicted vs. observed
plot(fit)
plot(fit, resid=TRUE)

fit.xv <- rioja:::crossval(fit, cv.method="boot", nboot=1000)
rand.t.test(fit.xv)
autoplot(fit.xv)

nn2<-n2(spp, "species")

output5<-cbind(fit.xv$coefficients, nn2)

#write.csv(output5, file="Supplementary_Table1b.csv")
#half of table s1

plot(fit.cv)
plot(fit, resid=TRUE)

#tiff("Fig3a.tiff", width = 6, height = 6, units = 'in', res = 300)

plot(fit.cv, xval=TRUE,pch = 21, bg = myColors[Cat])

with(Kterms3, legend("bottomright", legend = levels(Cat), bty = "n",
                     col = "black", pch = 21, pt.bg = myColors))

#dev.off()
#tiff("Fig3b.tiff", width = 6, height = 6, units = 'in', res = 300)
plot(fit.cv, xval=TRUE, resid=TRUE,pch = 21, bg = myColors[Cat])
#dev.off()

fit.xv
rand.t.test(fit.xv)

#code below needs to be edited for cores
#setwd("C:\\Users\\Dal User\\Desktop\\northwestern\\greenland")
#coreinput2<-read.csv(file="delta4.csv", row.names=1)
#coreinput2<-read.csv(file="north4.csv", row.names=1)
coreinput2<-read.csv(file="lc4.csv", row.names=1)
#coreinput2<-read.csv(file="fish4.csv", row.names=1)
#coreinput2<-read.csv(file="N14_2022.csv", row.names=1)

coresp<-coreinput2[, -cbind(1:3) ] #use this for cores
#coresp<-coreinput2[, -cbind(1:1) ] #use this for surface

coresp1 <- coresp / rowSums(coresp) * 100
coresp2<-coresp1[, -cbind(118:124) ]
cols_to_keep <- intersect(colnames(spp1),colnames(coresp2))
core <- coresp2[,cols_to_keep, drop=FALSE]

AD <- analogue_distances(spp1, core)
autoplot(AD, df = data.frame(age = as.numeric(coreinput2$age)), x_axis = "age") +
  labs(y = "Squared-chord distance")

df<-env1[(as.character(env$name) %in% as.character(species$name)), ]


library(palaeoSig)
library(analogue)
library(ggpalaeo)
library(rioja)
library(tidyverse)


env7<-env5$WC_SUMMER

fit <- WAPLS(sqrt(spp1), env7, nlps=2, nboot=9999)
fit
# cross-validate model
#rioja:::crossval(mod)
fit.cv <- rioja:::crossval(fit, cv.method="loo")
# How many components to use?
rand.t.test(fit.cv)
screeplot(fit.cv)
autoplot(fit.cv, npls=2)

table1b<-cbind(fit.cv$fitted.values, env7)

write.csv(table1b, file="table1b.csv")

#predict the core
pred <- predict(fit, core, npls=2)

#plot predictions - depths are in rownames
depth <- as.numeric(coreinput2$age)
plot(depth, pred$fit[, 1], type="b", ylab="Predicted SumSST", las=1)

# predictions with sample specific errors
# }
# NOT RUN {
pred <- predict(fit, core, npls=2, sse=TRUE, nboot=1000)
pred
plot(depth, pred$fit[, 2], type="b", ylab="Predicted Mean Summer Temperature", las=1)
arrows(depth, pred$fit[,2]-pred$v1.boot[,2], depth, pred$fit[,2]+pred$v1.boot[,2], length=0.05, angle=90, code=3)

# }

plot(depth, pred$fit[, 1], type="b", ylab="Predicted Mean Summer Temperature", las=1)
arrows(depth, pred$fit[,1]-pred$v1.boot[,1], depth, pred$fit[,1]+pred$v1.boot[,1], length=0.05, angle=90, code=3)



output<-cbind(pred$fit, pred$v1.boot)

write.csv(output, file="reconstruction2.csv")



#tiff("Fig3a.tiff", width = 6, height = 6, units = 'in', res = 300)

plot(fit.cv, xval=TRUE,pch = 21, bg = myColors[Cat])

with(Kterms3, legend("bottomright", legend = levels(Cat), bty = "n",
                     col = "black", pch = 21, pt.bg = myColors))

#dev.off()
#tiff("Fig3b.tiff", width = 6, height = 6, units = 'in', res = 300)
plot(fit.cv, xval=TRUE, resid=TRUE,pch = 21, bg = myColors[Cat])


#dev.off()


#temp<-pred$fit.boot

#WA model if wanting to compare
#fit <- WA(spp1, env6, tolDW=TRUE)
# plot predicted vs. observed
#plot(fit)
#plot(fit, resid=TRUE)
#plot the reconstructio
#plot(depth, pred$fit[, 1], type="b")
#arrows(depth, pred$fit[,1]-pred$v1.boot[,1], depth, pred$fit[,1]+pred$v1.boot[,1], length=0.05, angle=90, code=3)


###########################################################################
### this script is to run goodness-of-fit and analogue analyses 

spp<-spp1
core_analog <- core

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(spp >= M) >= N
spp_red <- spp[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

#here I identify the columns that we keep based on our deletion criteria 
#this allows us to filter the core such that we only include taxa in the model
#the rest are deleted.

cols_to_keep <- intersect(colnames(spp_red),colnames(core_analog))

core_gof <- core_analog[,cols_to_keep, drop=FALSE]
ncol(core_gof)


#env7 is the recon variable

# chronological axis for plotting analogue results
chron <- cbind(coreinput2[,1] , coreinput2[,2])
colnames(chron) <- c("Depth","Age")
chron=as.data.frame(chron)

require("ggpalaeo")
require("analogue")
require("ggplot2")

spp_T <- decostand(spp_red, method="hellinger")
core_T <- decostand(core_gof, method="hellinger")

#goodness-of-fit residuals
rlens <- residLen(spp_T, env7, core_T)
autoplot(rlens, df = data.frame(age = as.numeric(chron$Age)), x_axis = "age", fill = c("azure4", "grey","white")) +
  labs(x = "age", y = "Squared residual distance", fill = "Goodness of fit", categories = c("Good", "Fair", "Poor"))


#tiff("GOF_core2.tiff", width = 6, height = 6, units = 'in', res = 150)
#rlens <- residLen(spp_T, env7, core_T)
#autoplot(rlens, df = data.frame(age = as.numeric(chron$Age)), x_axis = "age", fill = c("azure4", "grey","white")) +
#  labs(x = "age", y = "Squared residual distance", fill = "Goodness of fit", categories = c("Good", "Fair", "Poor"))
#dev.off()

#analogue distance
## squared chord lengths for Core 
AD <- analogue_distances(spp, core_analog)
autoplot(AD, df = data.frame(age = as.numeric(chron$Age)), x_axis = "age", fill = c("azure4", "grey","white")) +
  labs(x = "Year", y = "Squared-chord distance")

#tiff("ANALOGUE_core2.tiff", width = 6, height = 6, units = 'in', res = 150)
#AD <- analogue_distances(spp, core_analog)
#autoplot(AD, df = data.frame(age = as.numeric(chron$Age)), x_axis = "age", fill = c("azure4", "grey","white")) +
#  labs(x = "Year", y = "Squared-chord distance")
#dev.off()

mod <- timetrack(spp, core_analog, transform = "hellinger",
                 method = "rda", scaling = "symmetric", correlation=TRUE)
mod
plot(mod, pch=21, bg=myColors3[Cat])

ordisurf(mod, env5$July, add=T, col="black", lwd=1)

## illustrating use of the formula

mod2 <- timetrack(spp, core_analog, env = data.frame(WC_Temp = env5$WC_SUMMER),
                  transform = "hellinger", method = "cca",
                  formula = ~ WC_Temp)
mod2


plot(mod2, ptype="b",  lwd = 2, pch=21, bg=myColors3[Cat])
ordisurf(mod2, env5$WC_SUMMER, add=T, col="black", lwd=1)

spe.cca.pars <-cca(spp ~ WC_SUMMER + Elev, data=env5)

fit <- envfit(spe.cca.pars ~ Elev+WC_SUMMER, data = env5, perm=999, display="lc")
plot(fit, lwd=2, col = "blue", arrow=1.5, lwd=2, lty="dashed", add=TRUE)


legend("bottomleft", inset=.02, title="Lakes",
       c("Canada", "Iceland", "Greenland", "Svalbard"), fill=c("black", "green", "red", "lightblue"), horiz=FALSE, cex=0.8)


library(survival)
library(grid)
library(gridGraphics)

plot(mod2, type = "n", ptype = "b", xlim=c(-3,3), ylim=c(-3,3))

# capture the plotted output as a grob
grid.echo()
grid.grab() -> k

# pull out the data from the grob..
k$children$`graphics-plot-1-points-1`$x -> x
k$children$`graphics-plot-1-points-1`$y -> y

#plot(mod, type = "n", ptype = "b")
#points(mod, which = "ordination", col = "grey", pch = 19, cex = 0.7)
#lines(x,y, col=1)
#points(x,y)
plot(mod2, xlim=c(-3,3), ylim=c(-3,3), type="n", ptype="n", lwd=3)
points(mod2, which = "ordination",  col=myColors3[Cat], pch = 21, bg=myColors3[Cat], cex = 0.7)
lines(x,y, col="magenta4", lty=1, lwd=2)


ordisurf(mod2, env5$WC_SUMMER, add=T, col="black", lwd=1)

spe.cca.pars <-cca(spp ~ WC_SUMMER + Elev, data=env5)

fit <- envfit(spe.cca.pars ~ Elev+WC_SUMMER, data = env5, perm=999, display="lc")
plot(fit, lwd=2, col = "blue", arrow=1.5, lwd=2, lty="dashed", add=TRUE)


legend("bottomleft", inset=.02, title="Lakes",
       c("Canada", "Iceland", "Greenland", "Svalbard"), fill=c("black", "green", "red", "lightblue"), horiz=FALSE, cex=0.8)




#tiff("CORETRAJECTORY.tiff", width = 6, height = 6, units = 'in', res = 150)
#plot(mod2, xlim=c(-3,3), ylim=c(-3,3), type="n", ptype="n", lwd=3)
#points(mod2, which = "ordination",  col=myColors3[Cat], pch = 21, bg="white", cex = 0.7)
#lines(x,y, col="magenta4", lty=1, lwd=2)
#ordisurf(mod2, env5$WC_SUMMER, add=T, col="black", lwd=1)
#spe.cca.pars <-cca(spp ~ WC_SUMMER + Elev, data=env5)
#fit <- envfit(spe.cca.pars ~ Elev+WC_SUMMER, data = env5, perm=999, display="lc")
#plot(fit, lwd=2, col = "blue", arrow=1.5, lwd=2, lty="dashed", add=TRUE)
#legend("bottomleft", inset=.02, title="Lakes",
#       c("Canada", "Iceland", "Greenland", "Svalbard"), fill=c("black", "green", "red", "lightblue"), horiz=FALSE, cex=0.8)
#dev.off()



#Run Telford analysis
#Telford, R. J., & Birks, H. J. B. (2011). 
#A novel method for assessing the statistical significance of quantitative reconstructions 
#inferred from biotic assemblages. Quaternary Science Reviews, 30(9-10), 1272-1278.
library(palaeoSig)
require(rioja)

rlghr <- randomTF(spp = spp_T, env =  env7,
                  fos = core_T, n = 999, fun = WAPLS,  col = 2)
rlghr$sig
plot(rlghr, "Temp")

require("ggplot2")
autoplot(rlghr, "Temp")
#
#
######################################################
# analogue, passively plot core
#build core model using MAT
mod <- MAT(spp_T, env7)
autoplot(mod)

#timetrack supplemental

mod <- timetrack(spp, core_gof, transform = "hellinger",
                 method = "rda")

mod
plot(mod)
ordisurf(mod, env7, add=T, col="black", lwd=2)

#build the model
fit <- WAPLS(spp_T, env7, nlps=2)

fit
# cross-validate model
#rioja:::crossval(mod)
fit.cv <- rioja:::crossval(fit, cv.method="loo")
# How many components to use?
rand.t.test(fit.cv)
screeplot(fit.cv)


tiff("model_c1.tiff", width = 6, height = 6, units = 'in', res = 150)

plot(fit.cv, xval=TRUE)

dev.off()


tiff("model_c2.tiff", width = 6, height = 6, units = 'in', res = 150)

plot(fit.cv, xval=TRUE, resid=TRUE)


dev.off()



#predict the core
pred <- predict(fit, core_T, npls=2, sse=TRUE, nboot=1000)
Age<-chron$Age
plot(Age, pred$fit[,1], type="b", ylim=range(c(pred$fit[,1]-pred$v1.boot[,1], pred$fit[,1]+pred$v1.boot[,1])))
arrows(Age, pred$fit[,1]-pred$v1.boot[,1], Age, pred$fit[,1]+pred$v1.boot[,1], length=0.05, angle=90, code=3)

output<-cbind(pred$fit, pred$v1.boot)
write.csv(output, file="reconstruction2.csv")
#
#Pete asked for
#plot pete asked for
vare.dca3 <- decorana(core)
vare.dca.sc <- scores(vare.dca3, display=c("sites"), choices=c(1,2))
dat3<-cbind(pred$fit, vare.dca.sc)

tiff("delta_supp2.tiff", width = 10, height = 8, units = 'in', res = 150)

plot(dat3[,1],dat3[,6],pch = 20,col = 1, type="b", ylab="DCA1", xlab="Predicted Mean Summer Temperature")

dev.off()
#
