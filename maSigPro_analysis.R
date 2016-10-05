
# /////////////////////////////////////////////////////////////////////////////
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#                             Erycina - maSigPro

# /////////////////////////////////////////////////////////////////////////////
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# Michelle Hwang
# September 2016

library(maSigPro)
library(Mfuzz)
library(ggplot2)
library(reshape2)

options(stringsAsFactors=F)
setwd("~/Google Drive/UGA/Research/Erycina\ 2/")

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# The see.genes function in maSigPro has an error that does not allow it
# to perform fuzzy clustering. This is a fix provided by Karolina Heyduck.

see.genes.kh<-function (data, edesign = data$edesign, time.col = 1, repl.col = 2, 
                        group.cols = c(3:ncol(edesign)), 
                        names.groups = colnames(edesign)[3:ncol(edesign)], 
                        cluster.data = 1, groups.vector = data$groups.vector, 
                        k = 9, m = 1.45, cluster.method = "hclust", 
                        distance = "cor", agglo.method = "ward.D", 
                        show.fit = FALSE, dis = NULL, step.method = "backward", 
                        min.obs = 3, alfa = 0.05, nvar.correction = FALSE, 
                        show.lines = TRUE, iter.max = 500, 
                        summary.mode = "median", color.mode = "rainbow", 
                        cexlab = 1, legend = TRUE, newX11 = TRUE, ylim = NULL, 
                        main = NULL, 
                        ...) 
{
  time = edesign[, time.col]
  repvect = edesign[, repl.col]
  groups = edesign[, group.cols]
  narrays <- length(time)
  if (!is.null(dim(data))) {
    dat <- as.data.frame(data)
    clusterdata <- data
  }
  else {
    clusterdata <- data[[cluster.data]]
    dat <- as.data.frame(data$sig.profiles)
  }
  clusterdata <- clusterdata
  if (nrow(dat) > 1) {
    dat <- as.data.frame(dat[, (ncol(dat) - length(time) + 
                                  1):ncol(dat)])
    count.na <- function(x) length(x[is.na(x)])
    NAs <- apply(as.matrix(dat), 1, count.na)
    count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= 
                       2), ]
  }
  else {
    NAs <- 1
  }
  kdata <- NULL
  dcorrel<-NULL
  clus2<-NULL
  out <- TRUE
  if (nrow(dat) > 1) {
    if (cluster.data != 1 || cluster.data != "sig.profiles") {
      if (any(is.na(clusterdata))) 
        clusterdata[is.na(clusterdata)] <- 0
    }
    else if (is.na(all(dist(clusterdata) > 0)) || (cluster.method == 
                                                   "kmeans" & any(is.na(clusterdata))) || (distance == 
                                                                                           "cor" & any(sd(t(clusterdata), na.rm = TRUE) == 0))) {
      if (!is.null(kdata)) {
        clusterdata <- kdata
      }
      else {
        clusterdata <- NULL
      }
    }
    clusterdata <- clusterdata
    if (!is.null(clusterdata)) {
      k <- min(k, nrow(dat), na.rm = TRUE)
      if (cluster.method == "hclust") {
        if (distance == "cor") {
          dcorrel <- matrix(rep(1, nrow(clusterdata)^2), nrow(clusterdata), 
                            nrow(clusterdata)) - cor(t(clusterdata), 
                            use = "pairwise.complete.obs")
          clust <- hclust(as.dist(dcorrel), method = agglo.method)
          c.algo.used = paste(cluster.method, "cor", 
                              agglo.method, sep = "_")
        }
        else {
          clust <- hclust(dist(clusterdata, method = distance), 
                          method = agglo.method)
          c.algo.used = paste(cluster.method, distance, 
                              agglo.method, sep = "_")
        }
        cut <- cutree(clust, k = k)
      }
      else if (cluster.method == "kmeans") {
        cut <- kmeans(clusterdata, k, iter.max)$cluster
        c.algo.used = paste("kmeans", k, iter.max, sep = "_")
      }
      else if (cluster.method == "mfuzz") {
        n <- dim(clusterdata)[2]
        clusterdata[is.na(clusterdata)] <- 0
        temp <- tempfile()
        write.table(clusterdata, temp, quote = FALSE, 
                    sep = "\t", row.names = TRUE, col.names = TRUE)
        signif <- readExpressionSet(temp)
        signif.r<-filter.NA(signif, thres=0.25)
        signif.f<-fill.NA(signif.r, mode="mean")
        signif.s<-standardise(signif.f)
        cl <- mfuzz(signif.s, c = k, m = m)
        clus <- acore(signif.s, cl = cl, min.acore = (1/k))
        for (i in 1:k) {
          clus[[i]] <- transform(clus[[i]], cluster = i)
        }
        cut0 <- clus[[1]][, c(1, 3)]
        for (i in 2:k) {
          cut0 <- rbind(cut0, clus[[i]][, c(1, 3)])
        }
        cut <- transform(clusterdata, name = "")
        cut <- transform(cut, cluster = 0)
        cut <- cut[, c(n + 1, n + 2)]
        cut[, 1] <- rownames(cut)
        for (i in 1:dim(clusterdata)[1]) {
          cut[i, 2] <- cut0[cut[i, 1], 2]
        }
        cut <- cut[, 2]
        c.algo.used = paste("mfuzz", k, m, sep = "_")
      }
      else stop("Invalid cluster algorithm")
      if (newX11) 
        X11()
      groups <- as.matrix(groups)
      colnames(groups) <- names.groups
      if (k <= 4) 
        par(mfrow = c(2, 2))
      else if (k <= 6) 
        par(mfrow = c(3, 2))
      else if (k > 6) 
        rows<-round(k/4)
        par(mfrow = c(2+1, 4))
      for (i in 1:(k)) {
        PlotProfiles(data = dat[cut == i, ], repvect = repvect, 
                     main = i, ylim = ylim, color.mode = color.mode, 
                     cond = rownames(edesign), ...)
      }
      if (newX11) 
        X11()
      if (k <= 4) {
        par(mfrow = c(2, 2))
        cexlab = 0.6
      }
      else if (k <= 6) {
        par(mfrow = c(3, 2))
        cexlab = 0.6
      }
      else if (k > 6) {
        rows<-round(k/4)
        par(mfrow = c(rows+1, 4))
        cexlab = 0.35
      }
      for (j in 1:(k)) {
        PlotGroups(data = dat[cut == j, ], show.fit = show.fit, 
                   dis = dis, step.method = step.method, min.obs = min.obs, 
                   alfa = alfa, nvar.correction = nvar.correction, 
                   show.lines = show.lines, time = time, groups = groups, 
                   repvect = repvect, summary.mode = summary.mode, 
                   xlab = "time", main = paste("Cluster", j, sep = " "), 
                   ylim = ylim, cexlab = cexlab, legend = legend, 
                   groups.vector = groups.vector, ...)
      }
    }
    else {
      print("warning: impossible to compute hierarchical clustering")
      c.algo.used <- NULL
      cut <- 1
    }
  }
  else if (nrow(dat) == 1) {
    if (newX11) 
      X11()
    PlotProfiles(data = dat, repvect = repvect, main = NULL, 
                 ylim = ylim, color.mode = color.mode, cond = rownames(edesign), 
                 ...)
    if (newX11) 
      X11()
    PlotGroups(data = dat, show.fit = show.fit, dis = dis, 
               step.method = step.method, min.obs = min.obs, alfa = alfa, 
               nvar.correction = nvar.correction, show.lines = show.lines, 
               time = time, groups = groups, repvect = repvect, 
               summary.mode = summary.mode, xlab = "time", main = main, 
               ylim = ylim, cexlab = cexlab, legend = legend, groups.vector = groups.vector, 
               ...)
    c.algo.used <- NULL
    cut <- 1
  }
  else {
    print("warning: NULL data. No visualization possible")
    c.algo.used <- NULL
    cut <- NULL
  }
  OUTPUT <- list(cut, c.algo.used, groups, dcorrel)
  names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups", "correl")
  OUTPUT}

# ==============================================================================
# ANALYSIS
# ==============================================================================

# -----------
# Erycina C3 |
# ------------------------------------------------------------------------------

ec3 <- read.delim("../Erycina-C3/C3-fpkm4.isoforms.TMM.fpkm.matrix", row.names=1)
ec3[is.na(ec3)] <- 0
ec3  <- ec3[rowSums(ec3)!=0, ]

ec3_names    <- c('E10',	'E11',	'E12',	'E13',	'E14',	'E15',	'E16',	
                  'E17',	'E18',	'E19',	'E1',	'E20',	'E21',	'E22',	
                  'E23',	'E24',	'E25',	'E28',	'E3',	'E5',	'E6',	'E7',	
                  'E8',	'E9')
ec3_times    <- c('6AM', '10PM', '6PM', '2AM', '6PM', '2PM', '10AM', '6AM', 
                  '6PM', '10PM', '6AM', '10PM', '2AM', '2AM', '2AM',  '10AM', 
                  '2PM', '2PM', '2AM', '10AM', '10PM', '6AM', '10AM', '6PM')
ec3_times_v2 <- c(6, 22, 18, 2, 18, 14, 10, 6, 18, 22, 6, 22, 2, 2, 2, 10, 
                  14, 14, 2, 10, 22, 6, 10, 18)

colnames(ec3) <- ec3_names

info       <- as.data.frame(cbind(colnames(ec3), ec3_times_v2))
time_order <- info[order(ec3_times_v2),]


## Make design matrix
time             <- as.numeric(time_order[,2])
reps             <- c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,5,5,6,6,6,6)
group            <- rep(1,24)
design           <- NULL
design           <- cbind(time, reps, group)
rownames(design) <- time_order[,1]
designmat        <- make.design.matrix(design, degree=4)

## Calculate polynomials
fitc3 <- p.vector(ec3, designmat, family=negative.binomial(10))
fitc3$i # num of sig genes
fitc3$SELEC # exp matrix

## Remove influential genes
fitted         <- T.fit(fitc3)
influential    <- fitted$influ.info
inf.genenames  <- colnames(influential)
ec3 <- ec3[!rownames(ec3) %in% inf.genenames, ]

## Pick k clusters
newfit <- p.vector(ec3, designmat, family=negative.binomial(1))
wss    <- (nrow(newfit$SELEC)-1)*sum(apply(newfit$SELEC,2,var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(newfit$SELEC, centers=i, iter.max=20)$withinss)
}
quartz()
plot(1:15, wss, type="b")

## Estimate m for fuzzy clustering
NBgenes   <- ExpressionSet(assayData=newfit$SELEC)
NBgenes.f <- filter.std(NBgenes, min.std=0)
NBgenes.s <- standardise(NBgenes.f)
m1        <- mestimate(NBgenes.s)
m1
finalfit  <- T.fit(newfit)

## Plot
out               <- see.genes.kh(newfit$SELEC, edesign=design, 
                                  cluster.method="mfuzz", cluster.data=1, 
                                  k=6, m=m1)
out$cut           <- as.data.frame(out$cut)
rownames(out$cut) <- rownames(newfit$SELEC)
finaloutC3 <- out
finalfitC3 <- newfit$SELEC
out$cut["TR38229|c1_g2_i1",] #2, best pepc
out$cut["TR57996|c1_g1_i1",] #2
out$cut["TR74256|c0_g1_i1",] #1
out$cut["TR78785|c0_g1_i1",] #4
out$cut["TR50499|c0_g1_i1",] #NA, pepck
out$cut["TR67386|c0_g1_i3",] #NA, pepck
out$cut["TR51044|c0_g1_i1",] #1, ppdk
out$cut["TR26187|c7_g1_i1",] #5, ppdk

write.table(finaloutC3$cut, file="maSigPro-C3.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

# ------------
# Set1 + ND-M |
# ------------------------------------------------------------------------------
d_set1_ndm_all  <- read.delim("data/feb9/isoforms_v4.TMM.fpkm.matrix", row.names=1) 
d_set1_ndm_all[is.na(d_set1_ndm_all)] <- 0
d_set1_ndm_all  <- d_set1_ndm_all[ rowSums(d_set1_ndm_all)!=0, ]

sample_names_v3 <- c('F1','F3','F4','F5','F6','F7','F8','F9','F11','F13',
                     'F15','F18','F25','F26','F27','F28','F29','F30','F31',
                     'F32','F33','F34','F35','F36', 'E3','E4','E5','E6',
                     'E7','E8','E9','E10','E12','E14','E15','E16','E17',
                     'E18','E19','E20','E21','E43','E44','E45','E59','E60',
                     'label')
sample_times    <- c(c('10PM', '12AM', '4AM', '6AM', '10AM', '2AM', '4PM', 
                       '10PM','8PM', '6PM', '6AM', '12PM', '4AM', '4PM', 
                       '12PM', '2AM', '8AM', '2PM', '10AM', '8AM', '2PM', 
                       '6PM', '8PM', '12AM'), 
                     c('2AM', '2AM','6PM','10PM','2PM','2AM','2PM','2PM',
                       '2PM','6PM','6PM','6PM','6PM','10PM','10PM','10PM',
                       '10PM','10AM','10AM','10AM','6AM','6AM'))
sample_times_v2 <- c(c(22, 0, 4, 6, 10, 2, 16, 22, 20, 18, 6, 12, 4, 16, 
                       12, 2, 8, 14, 10, 8, 14, 18, 20, 0), 
                     c(2, 2, 18, 22, 14, 2, 14, 14, 14, 18, 18, 18, 18, 
                       22, 22, 22, 22, 10, 10, 10, 6, 6))

info       <- as.data.frame(cbind(colnames(d_set1_ndm_all), sample_times_v2))
time_order <- info[order(sample_times_v2),]

## Make design matrix
time             <- as.numeric(time_order[,2])
reps             <- c(1,1,2,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,6,7,7,8,8,8,
                      8,8,8,9,9,10,10,10,10,10,10,10,11,11,12,12,12,12,
                      12,12,12)
group            <- rep(1,46)
group2           <- rep(0,46)
design           <- NULL
design           <- cbind(time, reps, group)
rownames(design) <- time_order[,1]
designmat        <- make.design.matrix(design, degree=4)

## Calculate polynomials
fit <- p.vector(d_set1_ndm_all, designmat, family=negative.binomial(10))
fit$i # num of sig genes
fit$SELEC # exp matrix

## Remove influential genes
fitted         <- T.fit(fit)
influential    <- fitted$influ.info
inf.genenames  <- colnames(influential)
d_set1_ndm_all <- d_set1_ndm_all[!rownames(d_set1_ndm_all) %in% inf.genenames,]

## Pick k clusters
newfit <- p.vector(d_set1_ndm_all, designmat, family=negative.binomial(1))
wss    <- (nrow(newfit$SELEC)-1)*sum(apply(newfit$SELEC,2,var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(newfit$SELEC, centers=i, iter.max=20)$withinss)
}
quartz()
plot(1:15, wss, type="b")

## Estimate m for fuzzy clustering
NBgenes   <- ExpressionSet(assayData=newfit$SELEC)
NBgenes.f <- filter.std(NBgenes, min.std=0)
NBgenes.s <- standardise(NBgenes.f)
m1        <- mestimate(NBgenes.s)
m1
finalfit  <- T.fit(newfit)

## Plot
out               <- see.genes.kh(newfit$SELEC, edesign=design, 
                                  cluster.method="mfuzz", cluster.data=1, 
                                  k=10, m=1.05)
out$cut           <- as.data.frame(out$cut)
rownames(out$cut) <- rownames(newfit$SELEC)
out$cut["TR50389|c5_g3_i1",]
out$cut["TR87206|c1_g1_i1",]
out$cut["TR50976|c0_g1_i1",] #10
out$cut["TR99121|c1_g1_i1",]
out$cut["TR109755|c6_g1_i1",] #3


# -----
# Med | 
# ------------------------------------------------------------------------------

med <- read.delim("data/feb9/isoforms_med.TMM.fpkm.matrix", 
                  row.names=1, sep="\t")
med[is.na(med)] <- 0

sample_names_med         <- c('F1','F3','F4','F5','F6','F7','F8','F9','F11',
                              'F13','F15','F18','F25','F26','F27','F28','F29',
                              'F30','F31','F32','F33','F34','F35','F36','label') 
sample_times             <- c(c('10PM', '12AM', '4AM', '6AM', '10AM', '2AM', 
                                '4PM', '10PM','8PM', '6PM', '6AM', '12PM', 
                                '4AM', '4PM', '12PM', '2AM', '8AM', '2PM', 
                                '10AM', '8AM', '2PM', '6PM', '8PM', '12AM'))
sample_times_ordered_med <- c(22, 0, 4, 6, 10, 2, 16, 22, 20, 18, 6, 12, 4, 16, 
                              12, 2, 8, 14, 10, 8, 14, 18, 20, 0)

info       <- as.data.frame(cbind(colnames(med), sample_times_ordered_med))
time_order <- info[order(sample_times_ordered_med),]
data2      <- NULL
data2      <- med[,match(time_order[,1], colnames(med))]

## Make design matrix
time             <- as.numeric(time_order[,2])
reps             <- rep(1:12, each=2) 
group            <- rep(1,24)
group2           <- rep(0,24)
design           <- NULL
design           <- cbind(time, reps, group, group2)
rownames(design) <- time_order[,1]
designmat        <- make.design.matrix(design, degree=4)


## Calculate polynomials
fit2 <- p.vector(data2, designmat, family=negative.binomial(1))
fit2$i # num of sig genes
fit2$SELEC # exp matrix

## Remove influential genes
fitted2       <- T.fit(fit2)
influential   <- fitted2$influ.info
inf.genenames <- colnames(influential)
data2         <- data2[!rownames(data2) %in% inf.genenames, ]

## Pick k clusters
newfit2 <- p.vector(data2, designmat, family=negative.binomial(1))
wss2    <- (nrow(newfit2$SELEC)-1)*sum(apply(newfit2$SELEC,2,var))
for (i in 2:15) {
  wss2[i]<- sum(kmeans(newfit2$SELEC, centers=i, iter.max=20)$withinss)
}
quartz()
plot(1:15, wss, type="b")

## Estimate m for fuzzy clustering
NBgenes2   <- ExpressionSet(assayData=newfit2$SELEC)
NBgenes2.f <- filter.std(NBgenes2, min.std=0)
NBgenes2.s <- standardise(NBgenes2.f)
m2         <- mestimate(NBgenes2.s)
finalfit2  <- T.fit(newfit2)

## Plot
out2 <- see.genes.kh(newfit2$SELEC, edesign=design, cluster.method="mfuzz", 
                     cluster.data=1, k=10, m=1.1)
out2$cut           <- as.data.frame(out2$cut)
rownames(out2$cut) <- rownames(newfit2$SELEC)
out2$cut["TR50389|c5_g3_i1",]
out2$cut["TR87206|c1_g1_i1",]
out2$cut["TR50976|c0_g1_i1",] #3
out2$cut["TR99121|c1_g1_i1",]
out2$cut["TR109755|c6_g1_i1",] #10


# ------------------
# Med remove bad TA |
# ------------------------------------------------------------------------------

med <- read.delim("data/feb9/isoforms_med.TMM.fpkm.matrix", 
                  row.names=1, sep="\t")
med[is.na(med)] <- 0

sample_names_med         <- c('F1','F3','F4','F5','F6','F7','F8','F9',
                              'F11','F13','F15','F18','F25','F26','F27',
                              'F28','F29','F30','F31','F32','F33','F34',
                              'F35','F36','label') 
sample_times             <- c(c('10PM', '12AM', '4AM', '6AM', '10AM', '2AM', 
                                '4PM', '10PM','8PM', '6PM', '6AM', '12PM', 
                                '4AM', '4PM', '12PM', '2AM', '8AM', '2PM', 
                                '10AM', '8AM', '2PM', '6PM', '8PM', '12AM'))
sample_times_ordered_med <- c(22, 0, 4, 6, 10, 2, 16, 22, 20, 18, 6, 12, 4, 
                              16, 12, 2, 8, 14, 10, 8, 14, 18, 20, 0)

info  <- as.data.frame(cbind(colnames(med), sample_times_ordered_med))
data3 <- NULL
data3 <- med[,match(time_order[,1], colnames(med))]
data3['F3']  <- NULL
data3['F28'] <- NULL
data3['F6']  <- NULL
data3['F26'] <- NULL
time_order <- info[order(sample_times_ordered_med),]
time_order <- time_order[-grep("F3$", time_order$V1),]
time_order <- time_order[-grep("F28$", time_order$V1),]
time_order <- time_order[-grep("F6$", time_order$V1),]
time_order <- time_order[-grep("F26$", time_order$V1),]
data3[is.na(data3)] <- 0
data3 <- data3[rowSums(data3)!=0,]

## Make design matrix
time             <- as.numeric(time_order[,2])
reps             <- c(1,2,3,3,4,4,5,5,6,7,7,8,8,9,10,10,11,11,12,12)
group            <- rep(1,20)
design           <- NULL 
design           <- cbind(time, reps, group)
rownames(design) <- time_order[,1]
designmat        <- make.design.matrix(as.data.frame(design), degree=4)

## Calculate polynomials
# Doesn't work with negative binomial...stops @5100
fit3 <- p.vector(data3, designmat)
fit3$i # num of sig genes
fit3$SELEC # exp matrix

## Remove influential genes
fitted3       <- T.fit(fit3)
influential   <- fitted3$influ.info
inf.genenames <- colnames(influential)
data3         <- data3[!rownames(data3) %in% inf.genenames, ]

## Pick k clusters
newfit3 <- p.vector(data3, designmat)
wss3 <- (nrow(newfit3$SELEC)-1) * sum(apply(newfit3$SELEC,2,var))
for (i in 2:15) {
  wss3[i]<- sum(kmeans(newfit3$SELEC, centers=i, iter.max=20)$withinss)
}
quartz()
plot(1:15, wss3, type="b")

## Estimate m for fuzzy clustering
NBgenes3   <- ExpressionSet(assayData=newfit3$SELEC)
NBgenes3.f <- filter.std(NBgenes3, min.std=0)
NBgenes3.s <- standardise(NBgenes3.f)
m3         <- mestimate(NBgenes3.s)
finalfit3  <- T.fit(newfit3)
m3 #1.127358

## Plot
out3 <- see.genes.kh(newfit3$SELEC, edesign=design, cluster.method="mfuzz", 
                     cluster.data=1, k=6, m=m3)
out3$cut <- as.data.frame(out3$cut)
rownames(out3$cut) <- rownames(newfit3$SELEC)
out3$cut["TR50389|c5_g3_i1",]
out3$cut["TR87206|c1_g1_i1",]
out3$cut["TR50976|c0_g1_i1",] #6
out3$cut["TR99121|c1_g1_i1",]
out3$cut["TR109755|c6_g1_i1",] 

finaloutCAM <- out3
finaloutCAM$cut["TR50976|c0_g1_i1",] #6
finalfitCAM <- newfit3$SELEC

write.table(finaloutCAM$cut, file="maSigPro-NDM.txt", quote=FALSE, sep="\t", 
            row.names=TRUE, col.names=TRUE)


# ==============================================================================
# MY OWN PLOTS
# ==============================================================================

## Expression plots

# PEPC
time = c("12AM", "2AM", "4AM", "6AM", "8AM", "10AM", "12PM", "2PM", "4PM", 
         "6PM", "8PM", "10PM")
pepc.med.m <- rbind(
  cbind(rep("TR50976|c0_g1_i1",20), 
        t(data3["TR50976|c0_g1_i1",]), 
        time_order[,2]),
  cbind(rep("TR99121|c1_g1_i1",20), 
        t(data3["TR99121|c1_g1_i1",]), 
        time_order[,2]),
  cbind(rep("TR109755|c6_g1_i1",20), 
        t(data3["TR109755|c6_g1_i1",]), 
        time_order[,2])
)
rownames(pepc.med.m) <- 1:60
colnames(pepc.med.m) <- c('transcript', 'value', 'variable')
pepc.med.m <- as.data.frame(pepc.med.m)
pepc.med.m$value <- as.numeric(pepc.med.m$value)
pepc.med.m$variable <- as.numeric(pepc.med.m$variable)


(ggplot(pepc.med.m, aes(x=factor(variable), y=value, group=transcript)) +
  aes(colour=transcript) +
  stat_summary(fun.y="mean", geom="line") +
  geom_point(size=0.2, shape=c(0)) +
  labs(x="Time", y="Expression") +
  ggtitle("Expression of PEPC copies in E. pusilla (CAM)") +
  scale_x_discrete(labels=time) +
  scale_y_continuous(labels=seq(0, 6000, 1000), breaks=seq(0, 6000, 1000)) +
  theme(plot.margin=unit(c(1,1,1,1), "cm"),
      panel.background=element_rect(fill="white"),
      panel.grid.major=element_line(colour="grey", size=0.1),
      panel.border=element_rect(fill=NA, colour="black"),
      legend.background=element_rect(fill="white"),
      axis.title.x=element_text(vjust=0.5, hjust=0.5),
      legend.title=element_blank()
  ) 
)


# PEPC C3
time = c("2AM",  "6AM",  "10AM", "2PM",  "6PM", "10PM")
ec3 = ec3[,order(ec3_times_v2)]
pepc.c3.m <- rbind(
  cbind(rep("TR38229|c1_g2_i1",24), 
        t(ec3["TR38229|c1_g2_i1",]), 
        time_order[,2]),
  cbind(rep("TR57996|c1_g1_i1",24), 
        t(ec3["TR57996|c1_g1_i1",]), 
        time_order[,2]),
  cbind(rep("TR78785|c0_g1_i1",24), 
        t(ec3["TR78785|c0_g1_i1",]), 
        time_order[,2])
)
rownames(pepc.c3.m) <- 1:72
colnames(pepc.c3.m) <- c('transcript', 'value', 'variable')
pepc.c3.m <- as.data.frame(pepc.c3.m)
pepc.c3.m$value <- as.numeric(pepc.c3.m$value)
pepc.c3.m$variable <- as.numeric(pepc.c3.m$variable)


(ggplot(pepc.c3.m, aes(x=factor(variable), y=value, group=transcript)) +
  aes(colour=transcript) +
  stat_summary(fun.y="mean", geom="line") +
  geom_point(size=0.2, shape=c(0)) +
  labs(x="Time", y="Expression") +
  ggtitle("Expression of PEPC copies\n in E. crista-galli (C3)") +
  scale_x_discrete(labels=time) +
  scale_y_continuous(labels=seq(0, 1000, 100), breaks=seq(0, 1000, 100)) +
  theme(plot.margin=unit(c(1,1,1,1), "cm"),
        panel.background=element_rect(fill="white"),
        panel.grid.major=element_line(colour="grey", size=0.1),
        panel.border=element_rect(fill=NA, colour="black"),
        legend.background=element_rect(fill="white"),
        axis.title.x=element_text(vjust=0.5, hjust=0.5),
        legend.title=element_blank()
  ) 
)


## Cluster plot
#6
ndm.6.mat <- med[rownames(finaloutCAM$cut)
                 [sapply(finaloutCAM$cut, function(x) x == 6)],]
c3.6.mat  <- ec3[rownames(finaloutC3$cut)
                 [sapply(finaloutC3$cut, function(x) x == 6)],]

write.table(ndm.6.mat, "maSigPro-NDM-6.matrix", 
            row.names=TRUE, col.names=TRUE, quote=FALSE) 
write.table(c3.6.mat, "maSigPro-C3-6.matrix", 
            row.names=TRUE, col.names=TRUE, quote=FALSE) 


# ==============================================================================
# SAVED DATA
# ==============================================================================

m <- read.table("maSigPro-NDM.txt", header=FALSE, row.names=1)
n <- read.table("maSigPro-C3.txt", header=FALSE, row.names=1)


m["TR50976|c0_g1_i1",] #6, pepc
m["TR48647|c0_g1_i1",] #6, mdh
m["TR50499|c0_g1_i3",] #7, pepck

n["TR38229|c1_g2_i1",] #2, best pepc
n["TR57996|c1_g1_i1",] #2
n["TR78785|c0_g1_i1",] #4
n["TR51044|c0_g1_i1",] #1, ppdk
n["TR26187|c7_g1_i1",] #5, ppdk

canonical <- read.delim("canonical.txt", row.names=1) 
canonical <- canonical[order(row.names(canonical)),,drop=FALSE]
m         <- m[order(row.names(m)),,drop=FALSE]
goi.CAM   <- merge(m, canonical, by.x=0, by.y=0)
goi.CAM
n         <- n[order(row.names(n)),,drop=FALSE]
goi.C3    <- merge(n, canonical, by.x=0, by.y=0)
goi.C3


