## code to prepare `dur_john` dataset goes here

# Load data #

data <- read.table("data-raw/dur_john_raw.dat")

# Exclude missing variables and oil states #

k <- ncol(data)
indx <- as.matrix(1-(data[,5]== -999))%*%matrix(c(1),1,k)
data <- as.matrix(data[indx>0])
data <- matrix(data,nrow=nrow(data)/k,ncol=k)
indx <- as.matrix(1-(data[,6]== -999))%*%matrix(c(1),1,k)
data <- as.matrix(data[indx>0])
data <- matrix(data,nrow=nrow(data)/k,ncol=k)
indx <- as.matrix(1-(data[,10]== -999))%*%matrix(c(1),1,k)
data <- as.matrix(data[indx>0])
data <- matrix(data,nrow=nrow(data)/k,ncol=k)
indx <- as.matrix(1-(data[,11]== -999))%*%matrix(c(1),1,k)
data <- as.matrix(data[indx>0])
data <- matrix(data,nrow=nrow(data)/k,ncol=k)
indx <- as.matrix(data[,2]== 1)%*%matrix(c(1),1,k)
data <- as.matrix(data[indx>0])
data <- matrix(data,nrow=nrow(data)/k,ncol=k)

# Make data transformations #

diff <- log(data[,6])-log(data[,5])
q <- data[,5]
gdp60 <- log(q)
iony <- log(data[,9]/100)
pgro <- log(data[,8]/100+.05)
sch <- log((data[,10])/100)
lit <- data[,11]
dat <- cbind(diff, gdp60, iony, pgro, sch, q, lit)

# Convert to data frame #

df <- as.data.frame(dat)
colnames(df) <- c("GDPGwth", "LogGDP1960", "LogInvGDP", "LogPopGwth", "LogSchool",
                  "GDP1960", "Literacy")

# Export to .rda file
# (name of object passed to use_data() determines name of .rda file)

dur_john <- df
use_data(dur_john, overwrite = T)
