library(thrreg)

## Load data
data <- dur_john

## Define variable names for output (optional)
var_names <- c("GDP\\_Gwth", "Log\\_GDP\\_1960", "Log\\_Inv/GDP",
               "Log\\_Pop\\_Gwth", "Log\\_School", "GDP\\_1960", "Literacy")
xi <- 2:5
h <- 1


## First Level
sink("output/output.tex")
cat("\\section{First Sample Split}", "\n\n")

## Test for a Sample Split
# test_hom <- thr_test_hom(data, "GDPGwth", xi, "GDP1960", var_names)
test_1_gdp <- thr_test_het(data, "GDPGwth", xi, "GDP1960", var_names)
test_1_lit <- thr_test_het(data, "GDPGwth", xi, "Literacy", var_names)

## Estimate Sample Split
qhat_1 <- thr_est(data, "GDPGwth", xi, "GDP1960", h, test_1_gdp$p_value,
                  var_names, digits.thr = 0, header = var_names[1])
sink()


## Second Level
dat_2l <- subset(data, GDP1960 <= qhat_1)
dat_2u <- subset(data, GDP1960 > qhat_1)

sink("output/output.tex", append = T)
cat(paste0("\\section{Second Sample Split: Subsample, Incomes below ", qhat_1,
           "}", "\n\n"))

## Test for a Sample Split
test_2l_gdp <- thr_test_het(dat_2l, "GDPGwth", xi, "GDP1960", var_names)
test_2l_lit <- thr_test_het(dat_2l, "GDPGwth", xi, "Literacy", var_names)

cat(paste0("\\section{Second Sample Split: Subsample, Incomes above ", qhat_1,
           "}", "\n\n"))

## Test for a Sample Split
test_2u_gdp <- thr_test_het(dat_2u, "GDPGwth", xi, "GDP1960", var_names)
test_2u_lit <- thr_test_het(dat_2u, "GDPGwth", xi, "Literacy", var_names)

## Estimate Sample Split
qhat_2 <- thr_est(dat_2u, "GDPGwth", xi, "Literacy", h, test_2u_lit$p_value,
                  var_names, digits.thr = 0, header = var_names[1])
sink()


## Third Level
dat_3l <- subset(dat_2u, Literacy <= qhat_2)
dat_3u <- subset(dat_2u, Literacy > qhat_2)

sink("output/output.tex", append = T)
cat(paste0("\\section{Third Sample Split: Subsample, Incomes above ", qhat_1,
           ", Literacy below ", qhat_2, "}", "\n\n"))

## Test for a Sample Split
test_3l_gdp <- thr_test_het(dat_3l, "GDPGwth", xi, "GDP1960")
test_3l_lit <- thr_test_het(dat_3l, "GDPGwth", xi, "Literacy")

cat(paste0("\\section{Third Sample Split: Subsample, Incomes above ", qhat_1,
           ", Literacy above ", qhat_2, "}", "\n\n"))

## Test for a Sample Split
test_3u_gdp <- thr_test_het(dat_3u, "GDPGwth", xi, "GDP1960")
test_3u_lit <- thr_test_het(dat_3u, "GDPGwth", xi, "Literacy")

sink()
