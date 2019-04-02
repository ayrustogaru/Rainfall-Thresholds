#Current Directory###
current_dir <- "C:/Users/Surya/Desktop/Landslides/" 

print(getwd())
setwd(current_dir)
print(getwd())

inputFile <- paste(current_dir,"OUTPUT/Reconstructed rainfall conditions/MPRC.csv",sep="")
sheet <- read.csv(inputFile, header=TRUE, sep=";")

duration_total <- sheet$D_L
# duration_c <- subset(sheet, Place == "Chukha")$Duration_h
# duration_g <- subset(sheet, Place == "Gedu")$Duration_h
# duration_m <- subset(sheet, Place == "Malbase")$Duration_h

event_total <- sheet$E_L
# event_c <- subset(sheet, Place == "Chukha")$Intensity
# event_g <- subset(sheet, Place == "Gedu")$Intensity
# event_m <- subset(sheet, Place == "Malbase")$Intensity

print(duration_total)
print(event_total)

#Defining constants#######

 HOUR_TO_SECONDS = 3600
 lockBinding("HOUR_TO_SECONDS" ,globalenv())

#Conversion of the vectors#################################
 

##Plotting the events on an E-D Graph####################
plot(duration_total,event_total, xaxt = "n", yaxt = "n", log = "xy", xlab = "Duration, D(h)", ylab = "Cumulated Event Rainfall, E(mm)", pch = 21,
	col = "darkorange", bg = "darkorange", xlim = c(1,1000), ylim = c(1,1000))

axis(side=1, at=c(1,10,100,1000), tck=-0.04,las=0)
axis(side=1, at=c(seq(2,9,1),seq(20,90,10),seq(200,900,100)),  tck=-0.02, labels=FALSE)

axis(side=2, at=c(1,10,100,1000), tck=-0.04,las=0)
axis(side=2, at=c(seq(2,9,1),seq(20,90,10),seq(200,900,100)),  tck=-0.02, labels=FALSE)


#points(duration_g, event_g, pch = 21, col = "red", bg = "red")
#points(duration_m, event_m, pch = 21, col = "gray", bg = "gray")

#Finding the Linear fit using Least Squares##############
fit <- lm(log10(event_total) ~ log10(duration_total))
alpha50 = fit$coefficients[[1]]
beta = fit$coefficients[[2]]
#abline(fit)
coefficients(fit)

# #Calculation of the Intensity along the linear fit (T50)######
N = length(event_total)
lgIf <- vector(mode = "numeric", length = N)
del <- vector(mode = "numeric", length = N)

for(i in 1:N)
{
	lgIf[i] <- fit$coefficients[[2]]*log10(duration_total[i]) + fit$coefficients[[1]]
}
#Calculating the delta#########################################
for(i in 1:N)
{
	del[i] <- log10(event_total[i]) - lgIf[i]
}

#Probability Density Function of del##################################
pdf_del <- density(del, bw = "nrd0", adjust = 1, kernel = "gaussian")
a = min(del)
b = max(del)
print(pdf_del)
#plot(pdf_del) #Uncomment to plot the pdf

# #Delta fitted with a Gaussian function##########################
yl<-pdf_del$y
t<-pdf_del$x

# #z = min(pdf_del$x)
# #y = k*exp(-0.5*((z-m)/s)^2)
# #print(y)

nlm <- nls(yl ~ k*exp(-0.5*((t-m)/s)^2), start = list(k = 1.1, m = 0.0, s = 0.4))
#plot(t,predict(nlm))
summary(nlm)

k <- coef(nlm)[1]
m <- coef(nlm)[2]
s <- coef(nlm)[3]

func<- function(x) 
{
	k*exp(-0.5*((x-m)/s)^2)
}
#curve(k*exp(-0.5*((x-m)/s)^2), from = -1.05, to = 1.05, n = 100)


total_area <- integrate(function(x) k*exp(-0.5*((x-m)/s)^2), lower = -1, upper = 1)$value
lower <- -1
dx <- 0.01
upper_start <- lower+dx
area <- 0
while(area/total_area <=0.05)
{
	area <- integrate(function(x) k*exp(-0.5*((x-m)/s)^2), 
				lower = -1, upper = upper_start)$value
	upper_start <- upper_start + dx
}
area_5percent <- upper_start - dx
print(area_5percent)
del_star <- m - area_5percent
print(del_star)
alpha5 <- alpha50 - del_star
print(alpha5)
print(beta)
xco<-c(3,3,400,400)
yco<-c(10.5,5.4,75,190)
polygon(xco,yco,col=rgb(1.0,0.6,0.1,0.15),border = NA)
clip(3,400,1,1000)
abline(0.6127838, 0.5585236,col = "darkorange",lwd= "2")
sa=4.1
si=1.5
sg=0.55
sib=0.09
mylabel = bquote(italic(E)== ~ "(" ~ .(format(sa, digits = 2,nsmall = 1)) ~ "±"~ .(format(si, digits = 1,nsmall = 1)) ~ ")" ~ 
                           italic(D)^( .(format(sg, digits = 3,nsmall = 1)) ~ "±" ~ .(format(sib, digits = 1,nsmall = 1)))            )
        
        
        
        text(md,1.41, mylabel , lty=1, bty="n", col="darkred",adj=1,cex=1.2)

#eq<-expression(paste(bold("E = 3.72.D")^0.67))
#mtext(eq,3,line=-16,cex=1.5)


# ecdf_D = ecdf(duration_total)
#print(ecdf_D)
#plot(ecdf_D, xlab = "Duration(s)", ylab = "ECDF")
#dev.off()