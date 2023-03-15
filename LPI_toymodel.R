# First install the package 'devtools' needed to install the 'RLPI' package from GitHub
install.packages("devtools")
# Load 'devtools' package
library(devtools)

# Install from main ZSL repository online
install_github("Zoological-Society-of-London/rlpi", dependencies=TRUE)
# Load the 'rlpi' package
library(rlpi)

###########################################################################
# Set the random seed so that the stochastic process is identical to the main manuscript
set.seed(41)

# Years of the simulation
years <- 1970:2020

# Number of species in the simulation
S <- 5000

# Create a vector for multiplicative and additive error
mult.e <- rep(c(0,0.01,0.03,0.05), each=4)
add.e <- rep(c(0,0.01,0.03,0.05), 4)

# Placeholders for LPI outputs (mean, upper and lower confidence bands)
lpi.means <- matrix(NA,nrow=16,ncol=length(years))
lpi.upper <- matrix(NA,nrow=16,ncol=length(years))
lpi.lower <- matrix(NA,nrow=16,ncol=length(years))

# Simulate a loop for each combination of additive and multiplicative error
for (m in 1:16){

  # Create a matrix for the 5000 populations, and define starting populations as 100 individuals
  N.mat <- matrix(NA,ncol=length(years), nrow=S); N.mat[,1] <- 100

  # Simulate population fluctuations with multiplicative and additive error for each species
  for (n in 1:S){
    m.err <- rnorm(length(years),0,mult.e[m])
    a.err <- rnorm(length(years),0,add.e[m])

    for (k in 2:length(years)) {
      N.mat[n,k] <- (N.mat[n,k-1]*exp(m.err[k]))+(N.mat[n,k-1]*a.err[k])
    }
  }

  # This code is to calculate the Living Planet Index for the simulated data.

  # !!! Important, you must have a folder named "LPI_files" in your working directory for this to work. The 'rlpi' code save the population infiles in this folder !!!

  # Give each time-series a unique integer ID
  ID <- 1:S
  # Give each species a unique integer ID
  Species <- as.factor(1:S)

  # Add the time-series ID to the simulated population data
  Pop_N <- cbind(ID,Species,N.mat)

  # Add column names to the dataset. The column names are the defaults needed by the 'rlpi' package.
  colnames(Pop_N) <- (c("ID","Binomial",paste0("X",as.factor(years))))

  # This is just an index vector for wich time-serie hould be included in calculations
  # Here, we include all time series (Note: this command in needed when you want to use a subset of the data)
  index_vector <- rep(TRUE, S)

  #This creates an 'infile' for calculating the LPI
  # Note: you need a folder names 'LPI_files' in your working directory
  infile_N <- create_infile(as.data.frame(Pop_N), start_col_name="X1970", end_col_name="X2020", index_vector=index_vector, name="LPI_files/lpi_N")
  lpi_N <- LPIMain(infile_N, REF_YEAR = 1970, PLOT_MAX = 2019, BOOT_STRAP_SIZE = 100, VERBOSE=FALSE, plot_lpi=TRUE)

  # Record the LPI values
  lpi.means[m,] <- lpi_N$LPI_final
  lpi.upper[m,] <- lpi_N$CI_high
  lpi.lower[m,] <- lpi_N$CI_low
}

################################################################################
################################################################################
# The following code generates the figure
# Define the .png file to be saved in your working directory
png(filename="FigureR1.png",width=32,height=16,units="cm",res=300)
# Outlay panels
par(mfrow=c(1,2))
# Define plot pargins
par(mai=c(0.6,0.65,0.4,0.10))

# Make a blank plot
plot(0,0,type="n", ylab="",xlab="",axes=F, ylim=c(0,100),xlim=c(-5,100))

# Define the colours for each combination 
trans <- 0.8 # This is a universal transparency modifier
col.a <- c(rgb(0.2,0.2,1,trans),rgb(1,0,0.6,trans),rgb(1,0.5,0,trans),rgb(0,0.8,0,trans),
           rgb(0,0,0.8,trans),rgb(0.8,0,0.6,trans),rgb(0.75,0.375,0,trans),rgb(0,0.6,0,trans),
           rgb(0,0,0.6,trans),rgb(0.6,0,0.6,trans),rgb(0.5,0.25,0,trans),rgb(0,0.4,0,trans),
           rgb(0,0,0.4,trans),rgb(0.4,0,0.6,trans),rgb(0.25,0.125,0,trans),rgb(0,0.2,0,trans))

# Define the limits for vertical (lower and upper) and horizontal (left and right) margins
vlow <- rep(c(70,50,30,10),each=4)
vup <- rep(c(90,70,50,30),each=4)
hl <- rep(c(10,30,50,70),4)
hr <- rep(c(30,50,70,90),4)

# Draw each combination as a coloured polygon
for (i in 1:16){
    polygon(c(hl[i],hl[i],hr[i],hr[i]),c(vlow[i],vup[i],vup[i],vlow[i]),col=col.a[i])
}

# Label the axes manually
text(20,94,"0", cex=1.1)
text(40,94,"0.01", cex=1.1)
text(60,94,"0.03", cex=1.1)
text(80,94,"0.05", cex=1.1)

text(9,20,"0.05", cex=1.1, pos=2)
text(9,40,"0.03", cex=1.1, pos=2)
text(9,60,"0.01", cex=1.1, pos=2)
text(9,80,"0", cex=1.1, pos=2)

text(50,98, cex=1.3, "Additive stochasticity")
text(-3,50, cex=1.3, "Multiplicative stochasticity", srt=90)

# Label the panel
mtext("a",cex=1.5, side = 3, adj = 0, line = 0.5,font=2)

# This parameter determines how many digit should be included when rounding down values
rd <- 3

# Add the mean LPI values
cex.lab <- 1.1 # Font size modifier
text(20,83,paste(round(lpi.means[1,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(40,83,paste(round(lpi.means[2,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(60,83,paste(round(lpi.means[3,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(80,83,paste(round(lpi.means[4,length(years)],rd)), col="white", font=2, cex=cex.lab)

text(20,63,paste(round(lpi.means[5,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(40,63,paste(round(lpi.means[6,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(60,63,paste(round(lpi.means[7,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(80,63,paste(round(lpi.means[8,length(years)],rd)), col="white", font=2, cex=cex.lab)

text(20,43,paste(round(lpi.means[9,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(40,43,paste(round(lpi.means[10,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(60,43,paste(round(lpi.means[11,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(80,43,paste(round(lpi.means[12,length(years)],rd)), col="white", font=2, cex=cex.lab)

text(20,23,paste(round(lpi.means[13,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(40,23,paste(round(lpi.means[14,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(60,23,paste(round(lpi.means[15,length(years)],rd)), col="white", font=2, cex=cex.lab)
text(80,23,paste(round(lpi.means[16,length(years)],rd)), col="white", font=2, cex=cex.lab)

# Add the percentage declines
text(20,77,paste(round((-1+lpi.means[1,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(40,77,paste(round((-1+lpi.means[2,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(60,77,paste(round((-1+lpi.means[3,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(80,77,paste(round((-1+lpi.means[4,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)

text(20,57,paste(round((-1+lpi.means[5,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(40,57,paste(round((-1+lpi.means[6,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(60,57,paste(round((-1+lpi.means[7,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(80,57,paste(round((-1+lpi.means[8,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)

text(20,37,paste(round((-1+lpi.means[9,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(40,37,paste(round((-1+lpi.means[10,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(60,37,paste(round((-1+lpi.means[11,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(80,37,paste(round((-1+lpi.means[12,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)

text(20,17,paste(round((-1+lpi.means[13,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(40,17,paste(round((-1+lpi.means[14,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(60,17,paste(round((-1+lpi.means[15,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)
text(80,17,paste(round((-1+lpi.means[16,length(years)])*100,rd-2),"%"), col="white", font=1, cex=cex.lab)

# Set the margins for the second panel
par(mai=c(0.75,0.75,0.4,0.10))
# Plot the axes and labels
plot(0,0,type="n",las=1,xlim=c(1970,2020),
     ylim=c(0.9,1.1),ylab="LPI (1970 = 1)",xlab="Year", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.5,0.6,0))

#Add the polygons for the confidence bands
trans <- 0.1 # Global transparency 10%
col.a <- c(rgb(0.2,0.2,1,trans),rgb(1,0,0.6,trans),rgb(1,0.5,0,trans),rgb(0,0.8,0,trans),
           rgb(0,0,0.8,trans),rgb(0.8,0,0.6,trans),rgb(0.75,0.375,0,trans),rgb(0,0.6,0,trans),
           rgb(0,0,0.6,trans),rgb(0.6,0,0.6,trans),rgb(0.5,0.25,0,trans),rgb(0,0.4,0,trans),
           rgb(0,0,0.4,trans),rgb(0.4,0,0.6,trans),rgb(0.25,0.125,0,trans),rgb(0,0.2,0,trans))
for (i in 1:16){
  # Include error bars as polygons
  polygon(c(seq(1970,2020),seq(2020,1970)), c(lpi.lower[i,],rev(lpi.upper[i,])),col=col.a[i],border=NA)
}

# Add mean LPI values
trans <- 1 # Non-transparent
col.a <- c(rgb(0.2,0.2,1,trans),rgb(1,0,0.6,trans),rgb(1,0.5,0,trans),rgb(0,0.8,0,trans),
           rgb(0,0,0.8,trans),rgb(0.8,0,0.6,trans),rgb(0.75,0.375,0,trans),rgb(0,0.6,0,trans),
           rgb(0,0,0.6,trans),rgb(0.6,0,0.6,trans),rgb(0.5,0.25,0,trans),rgb(0,0.4,0,trans),
           rgb(0,0,0.4,trans),rgb(0.4,0,0.6,trans),rgb(0.25,0.125,0,trans),rgb(0,0.2,0,trans))
for (i in 1:16){
  # Add mean line
  lines(c(1970:2020),lpi.means[i,],col=col.a[i],lwd=2)
}

# Label panel
mtext("b",cex=1.5, side = 3, adj = 0, line = 0.5,font=2)

# Close plot device and save to file
dev.off()
