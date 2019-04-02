library("bazar")
library("clue")
library("graphics")
library("stats")
library("rpart")
library("rattle")
library("mvc")

Heaviside <- 
    function(x, a = 0) 
{ 
	result = (sign(x-a) + 1)/2
    
    # Return Value:	
    result
}


###################################################################################################################################

#Bayesian Inference##

#Current Directory###
current_dir <- "C:/Users/Surya/Desktop/Landslides/Bayesian/" 

print(getwd())
setwd(current_dir)
print(getwd())

#Defining constants#######

 HOUR_TO_SECONDS = 3600
 lockBinding("HOUR_TO_SECONDS" ,globalenv())




inputFile <- paste(current_dir,"INPUT/ED.csv",sep="")
data <- read.csv(inputFile)
data_yes <- subset(data, Landslide == "YES")
data_no <- subset(data, Landslide == "NO")

duration_vector <- data$Duration
rainfall_intensity <- data$Intensity

duration_vector_yes <- data_yes$Duration
rainfall_intensity_yes <- data_yes$Intensity

duration_vector_no<- data_no$Duration
rainfall_intensity_no <- data_no$Intensity

#print(class(duration_vector))
#print(duration_vector)

#print(class(rainfall_intensity))
#print(rainfall_intensity_yes)
#print(rainfall_intensity_no)

#max_d <- max(data$Duration_h)
#max_i <- max(data$Intensity)

 for(i in 1:length(duration_vector))
 {
 	duration_vector[i] <- HOUR_TO_SECONDS*duration_vector[i]
 }

 for(i in 1:length(duration_vector_yes))
 {
 	duration_vector_yes[i] <- HOUR_TO_SECONDS*duration_vector_yes[i]
 }

 for(i in 1:length(duration_vector_no))
 {
 	duration_vector_no[i] <- HOUR_TO_SECONDS*duration_vector_no[i]
 }



 plot(duration_vector_yes, rainfall_intensity_yes,  main = "Duration vs Rainfall Intensity", xlab = "Duration(hours)",
  	ylab = "Rainfall Intensity(mm)", type = "p", pch = 21,col = "violet", bg = "blue", xlim = c(10000, 10000000), ylim = c(0.1,1000), log="xy")

 points(duration_vector_no, rainfall_intensity_no, pch = 21 ,col="red", bg = "red")

	c = 9.0909091*(duration_vector^(-0.2))
	print(c)

 points(duration_vector,c, pch =21, col = "black", bg = "black")
#legend(10000,1000, legend=c("Landslide", "No Landslide"), col=c("red", "blue"), lty=1:2, cex=0.8)

################################################################################################################################
N = length(duration_vector)
a = 1
b = 1
x = seq(from = 0, to = 1, by = 1/N)
Alpha = seq( from = 0.01 , to = 10 , by = 1/N ) 
Beta = seq(from = 0.1, to = 2, by = 1/N )
pAlpha = dunif(Alpha, min = 0.1, max=100)
pBeta = dunif(Beta	, min = 0.1, max = 2)
maxY1 = max(pAlpha)
maxY2 = max(pBeta)

#plot( Alph , pAlpha , type="l" , lwd=3 ,
#	 xlim=c(0,200) , ylim=c(0,maxY1) , cex.axis=1.2 ,
#	 xlab=bquote(theta) , ylab=bquote(p(theta)) , cex.lab=1.5 ,
#	 main="Prior" , cex.main=1.5 )

#plot( Beta , pBeta , type="l" , lwd=3 ,
#	 xlim=c(0,3) , ylim=c(0,maxY2) , cex.axis=1.2 ,
#	 xlab=bquote(theta) , ylab=bquote(p(theta)) , cex.lab=1.5 ,
#	 main="Prior" , cex.main=1.5 )
#################################################################################################################################



zID <- vector(mode="numeric", length = N)
muID <- vector(mode="numeric", length = N)
pIDGivenAB <- vector(mode = "integer", length = N)
pID_1 <- vector(mode = "integer",length = N)
pGivenData <- vector(mode = "integer",length = N)

for(j in 1:N)
{
	zID[j] = (1-(Alpha[j]*rainfall_intensity[j]*(duration_vector[j]^(-1*Beta[j]))))
	 print(zID[j])
	# print(Alpha[j])
	# print(Beta[j])
}

for(j in 1:N)
{
	muID[j] = (0.9*Heaviside(zID[j], a = 0) + 0.1)*exp(-0.5*abs(zID[j]))
	print(muID[j])
}



for(j in 1:N)
{
	pIDGivenAB[j] = dbern(x[j],muID[j])
}
print(pIDGivenAB)
for(j in 1:N)
{
	pID_1[j] = pIDGivenAB[j]*pAlpha[j]*pBeta[j]
}
pID = sum(pID_1)

print(pID)

for(j in 1:N)
{
	pGivenData[j] = (pAlpha[j]*pBeta[j]*pIDGivenAB[j])/pID
}

ind = which.max(pGivenData)

print(pGivenData)

print(Alpha[ind])
print(Beta[ind])

persp( x = Alpha , y = Beta,  z = outer(x,y,pGivenData)  ,phi = 45, theta = 45 )




