mydata <- read.csv("putative_enhancers.csv", header = TRUE)
fix(mydata)

#function from little book of r for multivariate analysis
calcWithinGroupsVariance <- function(variable,groupvariable)
  {
     # find out how many values the group variable can take
     groupvariable2 <- as.factor(groupvariable[[1]])
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # get the mean and standard deviation for each group:
     numtotal <- 0
     denomtotal <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidata <- variable[groupvariable==leveli,]
        levelilength <- length(levelidata)
        # get the standard deviation for group i:
        sdi <- sd(levelidata)
        numi <- (levelilength - 1)*(sdi * sdi)
        denomi <- levelilength
        numtotal <- numtotal + numi
        denomtotal <- denomtotal + denomi
     }
     # calculate the within-groups variance
     Vw <- numtotal / (denomtotal - numlevels)
     return(Vw)
  }
  
  #function from little book of r for multivariate analysis
  calcSeparations <- function(variables,groupvariable)
  {
     # find out how many variables we have
     variables <- as.data.frame(variables)
     numvariables <- length(variables)
     # find the variable names
     variablenames <- colnames(variables)
     # calculate the separation for each variable
     for (i in 1:numvariables)
     {
        variablei <- variables[i]
        variablename <- variablenames[i]
        Vw <- calcWithinGroupsVariance(variablei, groupvariable)
        Vb <- calcBetweenGroupsVariance(variablei, groupvariable)
        sep <- Vb/Vw
        print(paste("variable",variablename,"Vw=",Vw,"Vb=",Vb,"separation=",sep))
     }
  }
  
  #function from little book of r for multivariate analysis
  calcWithinGroupsCovariance <- function(variable1,variable2,groupvariable)
  {
     # find out how many values the group variable can take
     groupvariable2 <- as.factor(groupvariable[[1]])
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # get the covariance of variable 1 and variable 2 for each group:
     Covw <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidata1 <- variable1[groupvariable==leveli,]
        levelidata2 <- variable2[groupvariable==leveli,]
        mean1 <- mean(levelidata1)
        mean2 <- mean(levelidata2)
        levelilength <- length(levelidata1)
        # get the covariance for this group:
        term1 <- 0
        for (j in 1:levelilength)
        {
           term1 <- term1 + ((levelidata1[j] - mean1)*(levelidata2[j] - mean2))
        }
        Cov_groupi <- term1 # covariance for this group
        Covw <- Covw + Cov_groupi
     }
     totallength <- nrow(variable1)
     Covw <- Covw / (totallength - numlevels)
     return(Covw)
  }
  calcBetweenGroupsVariance <- function(variable,groupvariable)
  {
     # find out how many values the group variable can take
     groupvariable2 <- as.factor(groupvariable[[1]])
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # calculate the overall grand mean:
     grandmean <- mean(variable)
     # get the mean and standard deviation for each group:
     numtotal <- 0
     denomtotal <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidata <- variable[groupvariable==leveli,]
        levelilength <- length(levelidata)
        # get the mean and standard deviation for group i:
        meani <- mean(levelidata)
        sdi <- sd(levelidata)
        numi <- levelilength * ((meani - grandmean)^2)
        denomi <- levelilength
        numtotal <- numtotal + numi
        denomtotal <- denomtotal + denomi
     }
     # calculate the between-groups variance
     Vb <- numtotal / (numlevels - 1)
     Vb <- Vb[[1]]
     return(Vb)
  }
  
  #function from little book of r for multivariate analysis
  mosthighlycorrelated <- function(mydataframe,numtoreport)
  {
     # find the correlations
     cormatrix <- cor(mydataframe)
     # set the correlations on the diagonal or lower triangle to zero,
     # so they will not be reported as the highest ones:
     diag(cormatrix) <- 0
     cormatrix[lower.tri(cormatrix)] <- 0
     # flatten the matrix into a dataframe for easy sorting
     fm <- as.data.frame(as.table(cormatrix))
     # assign human-friendly names
     names(fm) <- c("First.Variable", "Second.Variable","Correlation")
     # sort and print the top n correlations
     head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
  }
library(corrgram)
corrgram(mydata[6:13], order=TRUE, lower.panel=panel.shade, upper.panel=panel.conf, diag.panel = NULL, label.pos = 0.5, font.labels = 3, cex.labels = 1, col.regions = colorRampPalette(c("green", "lightgreen", "orange", "darkred")), main = "Correlations between Datasets")

standardisedconcentrations <- as.data.frame(scale(mydata[5:13]))
mydata.pca <- prcomp(standardisedconcentrations) 
summary(mydata.pca)
screeplot(mydata.pca, type="lines")
mydata.pca$rotation[,1]
sum((mydata.pca$rotation[,1])^2)

#function from little book of r for multivariate analysis
calcpc <- function(variables,loadings)
  {
     # find the number of samples in the data set
     as.data.frame(variables)
     numsamples <- nrow(variables)
     # make a vector to store the component
     pc <- numeric(numsamples)
     # find the number of variables
     numvariables <- length(variables)
     # calculate the value of the component for each sample
     for (i in 1:numsamples)
     {
        valuei <- 0
        for (j in 1:numvariables)
        {
           valueij <- variables[i,j]
           loadingj <- loadings[j]
           valuei <- valuei + (valueij * loadingj)
        }
        pc[i] <- valuei
     }
     return(pc)
  }
  
  calcpc(standardisedconcentrations, mydata.pca$rotation[,1])
  
  plot(mydata.pca$x[,1],mydata.pca$x[,2], col = c("black", "darkred",  "forestgreen")[mydata$status], pch = 16, xlab = "PC1", ylab = "PC2", main = "PC1 vs PC2")
  #text(mydata.pca$x[,1],mydata.pca$x[,2], mydata$gene, cex=0.7, pos=4, col="blue") 
  
  
  library("MASS") 

knowndata <- subset(mydata, mydata$status!="putative")

standardisedconcentrations <- as.data.frame(scale(knowndata[5:13]))
knowndata.pca <- prcomp(standardisedconcentrations)
screeplot(knowndata.pca, type="lines")
plot(knowndata.pca$x[,1],knowndata.pca$x[,2], pch = 16, col = c("green", "blue", "purple", "red")[knowndata$expr], xlab = "PCA 1", ylab = "PCA 2", main = "PCA, functional vs. nonfunctional enhancers")
text(knowndata.pca$x[,1],knowndata.pca$x[,2], knowndata$expr, cex=0.5, pos=4, col="black")


knowndata.lda <- lda(knowndata$status ~ knowndata$distance_from_TSS + knowndata$Eisen_120min_Zelda_chip.seq + 
                                   knowndata$Furlong_Snail_chip.seq + knowndata$Furlong_Twist_chip.seq + 
								   knowndata$MacArthur_Dorsal_chip.chip + knowndata$MacArthur_Snail_chip.chip +
                                   knowndata$MacArthur_Twist_chip.chip + knowndata$Zeitlinger1_Twist_chip.seq + 
                                   knowndata$Zeitlinger2_Twist_chip.seq)

knowndata.lda
knowndata.lda$scaling[,1]


#function from little book of r for multivariate analysis
 calclda <- function(variables,loadings)
  {
     # find the number of samples in the data set
     as.data.frame(variables)
     numsamples <- nrow(variables)
     # make a vector to store the discriminant function
     ld <- numeric(numsamples)
     # find the number of variables
     numvariables <- length(variables)
     # calculate the value of the discriminant function for each sample
     for (i in 1:numsamples)
     {
        valuei <- 0
        for (j in 1:numvariables)
        {
           valueij <- variables[i,j]
           loadingj <- loadings[j]
           valuei <- valuei + (valueij * loadingj)
        }
        ld[i] <- valuei
     }
     # standardise the discriminant function so that its mean value is 0:
     ld <- as.data.frame(scale(ld, center=TRUE, scale=FALSE))
     ld <- ld[[1]]
     return(ld)
  }
  
 knowndata.lda.values$x[,1] 
 calclda(knowndata[5:13], knowndata.lda$scaling[,1])
 knowndata.lda.values <- predict(knowndata.lda, knowndata[5:13])
 
 
 
 
 #function from little book of r for multivariate analysis
 groupStandardise <- function(variables, groupvariable)
  {
     # find out how many variables we have
     variables <- as.data.frame(variables)
     numvariables <- length(variables)
     # find the variable names
     variablenames <- colnames(variables)
     # calculate the group-standardised version of each variable
     for (i in 1:numvariables)
     {
        variablei <- variables[i]
        variablei_name <- variablenames[i]
        variablei_Vw <- calcWithinGroupsVariance(variablei, groupvariable)
        variablei_mean <- mean(variablei)
        variablei_new <- (variablei - variablei_mean)/(sqrt(variablei_Vw))
        data_length <- nrow(variablei)
        if (i == 1) { variables_new <- data.frame(row.names=seq(1,data_length)) }
        variables_new[`variablei_name`] <- variablei_new
     }
     return(variables_new)
  }
  
  #groupstandardisedconcentrations <- groupStandardise(knowndata[5:13], knowndata$status)
  
  #knowndata.lda2 <- lda(knowndata$status ~ groupstandardisedconcentrations$distance_from_TSS + groupstandardisedconcentrations$Eisen_120min_Zelda_chip.seq  +
  #                           groupstandardisedconcentrations$Furlong_Snail_chip.seq + groupstandardisedconcentrations$Furlong_Twist_chip.seq  +
  #                           groupstandardisedconcentrations$MacArthur_Dorsal_chip.chip + groupstandardisedconcentrations$MacArthur_Snail_chip.chip +
  #                           groupstandardisedconcentrations$MacArthur_Twist_chip.chip + groupstandardisedconcentrations$knowndata$Zeitlinger1_Twist_chip.seq  +
  #                           groupstandardisedconcentrations$Zeitlinger2_Twist_chip.seq)

  calcSeparations(knowndata.lda.values$x,knowndata$status)

plot(knowndata.lda.values$x[,1], pch = 16, col = c("black", "red")[knowndata$status], ylab = "lda 1", main = "Linear discriminate analysis of functional vs. non-functional enhancers", cex.main = 0.85)
 
