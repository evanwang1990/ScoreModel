data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
str(chileancredit) # Quick description of the data
table(chileancredit$FlagGB) # Tabulate target variable

# Training and testing samples (Just some basic formality for Modeling)
chileancredit.train=subset(chileancredit,FlagSample==1 & !is.na(FlagGB))
chileancredit.test=subset(chileancredit,FlagSample==0 & !is.na(FlagGB))

# Package application
result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB") # Run and save result
result$ivtable # Tabulation and Information Value
result$iv # Information value
result$bands # Bins or bands
result$ctree # Decision tree from partykit


round <- function(x, unit)
{
  x_min <- floor(x/unit)*unit
  x_ <- ifelse(x <= x_min + 0.5 * unit, x_min, x_min + unit)
}

result_1 <- collapseLevel(x = chileancredit.train$TOB, y =chileancredit.train$FlagGB, 20, method = 'iv', mode = 'J', sourcefile = 'test.R', sqlfile = 'test.sql')


library(smbinning)
chileancredit.train=subset(chileancredit,FlagSample==1 & !is.na(FlagGB))
freqMatrix <- as.matrix(with(chileancredit.train, table(Performance, FlagGB)))
freqMatrix <- freqMatrix[freqMatrix[,1] + freqMatrix[,2] > 0, ]
CollapseZeroCells(freqMatrix, matrix(NA, 5, 2), "J")



tmp <- sample(20:100, 10, replace = T)
tmp <- matrix(c(tmp, round(tmp * (seq(1, 3, length.out = 10) + 0.5 * runif(1:10)) + sample(1:10, 10))), ncol = 2)
rownames(tmp) <- letters[1:10]



df <- data.frame(x1 = sample(1:1000, 1e4, replace = T),
               x2 = sample(c(letters[1:3], NA), 1e4, replace = T),
                 y  = c(sample(c(1, 1, 0), 5e3, replace = T), sample(c(1, 0, 0, 0), 5e3, replace = T)))
res <- bestCollapse(c('x1', 'x2'), y, df, 20, method = 'iv', mode = 'J', tracefile = 'trace.Rout', sqlfile = 'sql_code.sql')
collapseLevel(x = df$x2, y = df$y, levels = 20, method = 'iv', mode = 'A', minp = 0.05, sourcefile = 'test.R', sqlfile = 'test.sql')