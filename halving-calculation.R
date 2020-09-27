library(ggplot2)
library(beanplot)

# data: contains data for multiple patients
# columns: PatID, TIME, LRATIO, lQL
# PatID: unique ID number for a patient
# TIME: months, time since treatment start
# LRATIO: log10(BCR-ABL/ABL * 100)
# lQL: lower quantification limit
path_to_data <- "data.csv" 
data <- read.csv(path_to_data) 

# summary: shows for each patient treatment outcome
# columns: PatID, relapse
# PatID: unique ID number for a patient
# realpse: TRUE or FALSE
path_to_summary <- "summary.csv"
summary <- read.csv(summary)

# calculates halving time based on two 2D points
# x: TIME
# y: log10(LRATIO)
halving_time  <- function(x1,x2,y1,y2){
    -(x2-x1)/(y2-y1)*log10(2)
} 

# list to include all patients with sufficient data 
good_patients  <- c()

for(patient in unique(data$PatID)){
    pat  <- subset(data[1:3], PatID == patient) # subset for a single patient with lQL column excluded
    pat  <- na.omit(pat) # omit rows with NA in LRATIO or TIME
    pat$TIME <- pat$TIME - min(pat$TIME) # shift TIME variable to zero if there are negative values for it
    if (nrow(pat) < 3){next}
    cond1  <- nrow(subset(pat, TIME >= 0 & TIME < 2)) > 0 # at least 1 point with TIME between 0 and 2 
    cond2  <- nrow(subset(pat, TIME >=2 & TIME < 7)) > 0  # at least 1 point with TIME between 2 and 7
    cond3  <- sum(pat[1:3,]$LRATIO == sort(pat[1:3,]$LRATIO, decreasing = TRUE)) == length(pat[1:3,]$LRATIO) # first 3 points are in strictly decreasing order
    cond4  <- !is.na(subset(summary , PatID == patient)$relapse) # treatment outcome is known for the patient
    if (cond1 & cond2 & cond3 & cond4) {good_patients <- append(good_patients, patient)}
}

# subset of sufficient data and summary
data  <- subset(data, PatID %in% good_patients)
summary  <- subset(summary, PatID %in% good_patients)

# include halftime column in summary
summary$halftime  <- 0

for(patient in unique(data$PatID)){
    pat  <- subset(data, PatID == patient)
    k  <- 1 # defines distance between rows in subset to calculate halving time
    x1 = pat$TIME[1]
    x2 = pat$TIME[1+k]
    y1 = pat$LRATIO[1]
    y2 = pat$LRATIO[1+k]
    h  <- halving_time(x1,x2,y1,y2)
    # it is searching for 2 points for halving time being located between them (excludes 2 initial points that are too close)
    while (x1+h > x2){ 
        k  <- k + 1
        x1 = pat$TIME[1]
        x2 = pat$TIME[1+k]
        y1 = pat$LRATIO[1]
        y2 = pat$LRATIO[1+k]
        h  <- halving_time(x1,x2,y1,y2)
    }
    summary$halftime[summary$PatID == patient]  <- h
}

# rename labels for better perception 
summary$relapse  <- ifelse(summary$relapse == TRUE, 'Molecular relapse', 'Sustained TFR')
summary$relapse <- factor(summary$relapse, levels = c('Sustained TFR', 'Molecular relapse'),ordered = TRUE) #order labels for ggplot (optional)

# halftime is multiplied by (365/12) in order to convert months to days
ggplot(summary, aes(y=(365/12)*halftime, x = relapse)) + 
    geom_boxplot(fill=c('purple', 'blue')) + labs(y = 'halving time (days)', x = 'Status')+
    theme_classic()