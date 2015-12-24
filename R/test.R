
data(chileancredit, package = 'smbinning') # Load smbinning sample dataset (Chilean Credit)
chileancredit.train=subset(chileancredit,FlagSample==1 & !is.na(FlagGB))
attach(chileancredit.train)
chileancredit.train$FlagGB[FlagGB == 1 & Performance == '70: Never delinquent' & runif(length(FlagGB)) < 0.05] <- 0
chileancredit.train$FlagGB[FlagGB == 1 & Performance == '62: 1 x 1-29' & runif(length(FlagGB)) < 0.05] <- 0
chileancredit.train$FlagGB[FlagGB == 0 & Performance == '20: 1+ x 90+' & runif(length(FlagGB)) < 0.05] <- 1
chileancredit.train$Miss[FlagGB == 1] <- sample(c('a', NA), sum(FlagGB == 1), replace = T)
chileancredit.train$Miss[FlagGB == 0] <- sample(c('a', rep(NA, 10)), sum(FlagGB == 0), replace = T)
detach(chileancredit.train)

#factor
splitWOE_factor <- splitLevel(FlagGB ~ Performance, chileancredit.train)
catLog(splitWOE_factor)
collapseWOE_factor <- collapseLevel(FlagGB ~ Performance, chileancredit.train, method = "max_likehood")
catLog(collapseWOE_factor)

#numeric
splitWOE_numeric <- splitLevel(FlagGB ~ TOB, chileancredit.train)
catLog(splitWOE_numeric)
collapseWOE_numeric <- collapseLevel(FlagGB ~ TOB, chileancredit.train, method = "max_likehood")
catLog(collapseWOE_numeric)

#one level with missing values
splitWOE_onelevel <- splitLevel(FlagGB ~ Miss, chileancredit.train)
catLog(splitWOE_onelevel)
collapseWOE_onelevel <- collapseLevel(FlagGB ~ Miss, chileancredit.train)
catLog(collapseWOE_onelevel)


tmp <- sample(20:100, 10, replace = T)
tmp <- matrix(c(tmp, round(tmp * (seq(1, 3, length.out = 10) + 0.5 * runif(1:10)) + sample(1:10, 10))), ncol = 2)
rownames(tmp) <- letters[1:10]



df <- data.frame(x1 = sample(1:1000, 1e4, replace = T),
               x2 = sample(c(letters[1:3], NA), 1e4, replace = T),
                 y  = c(sample(c(1, 1, 0), 5e3, replace = T), sample(c(1, 0, 0, 0), 5e3, replace = T)))
res <- bestCollapse(c('x1', 'x2'), y, df, 20, method = 'iv', mode = 'J', tracefile = 'trace.Rout', sqlfile = 'sql_code.sql')
collapseLevel(x = df$x2, y = df$y, levels = 20, method = 'iv', mode = 'A', minp = 0.05, sourcefile = 'test.R', sqlfile = 'test.sql')
