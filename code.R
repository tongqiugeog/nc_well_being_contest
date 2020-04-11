library(magrittr)
library(data.table)
library(spBayes)
library(coda)
xdata = read.csv('all_xdata.csv', stringsAsFactors = FALSE) %>% setDT()
ydata = read.csv('WBI-NC-County.csv', stringsAsFactors = FALSE) %>% setDT()

get_county_names <- function(x){
  strsplit(x, ' ', fixed = TRUE)[[1]][1]
}

county_names = sapply(ydata$CNTY_NAME, get_county_names) %>% unname()
ydata[['CNTY_name']] = county_names

all_combine = merge(xdata, ydata, by.x = 'NAME', by.y = 'CNTY_name')
write.csv(all_combine, 'data_set.csv')
predictors <- c('civengidx','civinstidx','comdividx','ecnsegidx','econ2idx',
                'ecostabidx','emp1654pct','emp5575pct','engltwpct','hmv1kpct',
                'homeincmdn','homevalmdn','hvi40uprto','labmarkidx','labrfrcpct',
                'pbn1964pct', 'pbns65upct','phyinacpct','primsecidx','prsvt16pct',
                'pubasspct','sdoh1idx','secmrtgpct','sesclusidx','sknurseprp',
                'skyhouidx','soccomidx','stabastdx','stutch2rto')
response   <- c("Well.Being.Index.Score")


cor_list   <- list()
for (covariates in c(16:108)){
  y        <- all_combine[, response, with = FALSE]
  x        <- all_combine[, colnames(xdata)[covariates], with = FALSE]
  cor_list[[covariates-15]] <- cor.test(as.matrix(x),as.matrix(y))$estimate 
}
names(cor_list)             <- colnames(xdata)[16:108]
cor_list                    <- unlist(cor_list)
predictors                  <- cor_list[abs(cor_list)>0.3] %>% names()

get_predictor_name <- function(x){
  strsplit(x, '.', fixed = TRUE)[[1]][1]
}
predictors                  <- sapply(predictors, get_predictor_name) %>% unname()
predictors <- setdiff(predictors, 'homevalmdn')


formula    <- as.formula(paste(response, paste(predictors, collapse = '+'), sep = '~'))


p         <- length(predictors) + 1
n.samples <- 20000
starting  <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
tuning    <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.1  <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                 "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))
cov.model <- "exponential"
burn.in   <- 0.5*n.samples
set.seed(9025)
calibrate_sites <- sample(all_combine$FID_1, 80)
validate_sites  <- setdiff(all_combine$FID_1, calibrate_sites)
calibrate_data  <- all_combine[FID_1 %in% calibrate_sites]
validate_data   <- all_combine[FID_1 %in% validate_sites]


validate_cov    <- cbind(rep(1, nrow(validate_data)), validate_data[, predictors, with = FALSE]) %>% as.matrix()

m.1       <- spLM(formula, data = calibrate_data, coords = as.matrix(cbind(calibrate_data$lon, calibrate_data$lat)),
                  starting=starting,
                  tuning=tuning, priors=priors.1, cov.model=cov.model,
                  n.samples=n.samples, verbose=TRUE)

m.1.pred  <- spPredict(m.1, start = burn.in, thin = 10, pred.covars = validate_cov,
                       pred.coords = as.matrix(cbind(validate_data$lon, validate_data$lat)))


y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, mean)
quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}
y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, quant)
plot(validate_data$Well.Being.Index.Score, y.hat[2,], pch=19, cex=0.5, xlab="observed y", ylab="predicted y")
R     <- cor.test(validate_data$Well.Being.Index.Score,
                  y.hat[2,])$estimate
m.1       <- spRecover(m.1, start=burn.in, verbose=FALSE)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
table <-   round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
table <- table[2:13,1]
table <- as.data.frame(table)
table$name = rownames(table)
table$sig <- table$table<0

ggplot(table, aes(x=name, y = table, fill = sig)) + geom_bar(stat = 'identity')+ theme_minimal() + scale_fill_manual("legend", values = c("FALSE" = '#92c5de', "TRUE" = '#f4a582')) +ylab('importance')+xlab('Social Determinants of Health')


plot_df   <- as.data.frame(cbind(validate_data$Well.Being.Index.Score, y.hat[2,]))

ggplot(data = plot_df, aes_string(x = "V1", y = "V2")) +
  geom_point(alpha = 0.2, shape = 21, fill = 'blue', color = 'grey', size = 2) +
  
  coord_cartesian(xlim = c(50, 70), ylim = c(55, 70))+
  theme_minimal() +
  theme(axis.text    = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+xlab('observed social well-being')+ylab('predicted social well-being')
