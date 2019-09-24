
library(dplyr)
library(tidyr)
library(ggplot2)
library(feather)
library(randomForest)
library(lubridate)
library(stringr)
library(caret)
library(zoo)
library(ade4)
library(ModelMetrics)

## Reading test data
test = read.csv("Data/Test.csv",h=T,stringsAsFactors=F,na.strings="")
test$time = as.POSIXct(test$time,format ="%d/%m/%Y %H:%M:%S",tz="UTC")


## Reading training data
train = read.csv("Data/Train.csv",h=T,stringsAsFactors=F,na.strings="")
train$time = as.POSIXct(train$time,format ="%d/%m/%Y %H:%M:%S",tz="UTC")
## Filtering out rows where indicator1<0.8
train = filter(train,indicator1<0.8 | is.na(indicator1))

## Fill NA values in test set with either the median value in
test_medianfill or the next value in test_nextfill
test_medianfill = test
for(column in str_subset(names(test),"indicator"))
{
  val = median(test[[column]],na.rm=T)
  test_medianfill[[column]][is.na(test_medianfill[[column]])] = val
}

test_nextfill = test %>% mutate_all(funs(na.locf(.,na.rm=F,fromLast=T)))

## creating a smaller training dataset keeping only rows where the water
quality is not NA
train_small = train[!is.na(train$waterQuality),]


## Plotting data

## Histogram of water quality values
pdf("Plots/Hist_waterQuality.pdf")
ggplot(train_small) + geom_histogram(aes(waterQuality),bins=100) + xlab("Water quality")
dev.off()

## Water quality throughout time
pdf("Plots/WaterQuality_throughout_time.pdf")
ggplot(train_small, aes(time,waterQuality)) + geom_point(alpha=0.5) + geom_smooth(col=2) + xlab("Time") + ylab("Water quality")
dev.off()

## Water quality in weekdays
pdf("Plots/WaterQuality_in_weekdays.pdf")
ggplot(train_small) + geom_boxplot(aes(weekdays(time),waterQuality))
dev.off()

## Divide water quality into 4 categories:
## w<1, 1<=w<2, 2<=w<3, w>=3
## Plot the categories throughout time monthly

## We first add a column with a date corresponding to the first day of the month for each row
train_small$firstDayMonth = as.Date(train_small$time)
mday(train_small$firstDayMonth) = 1

## First plot is with number of monthly values
pdf("Plots/waterQuality_category_monthly.pdf")
train_small %>% mutate(waterCat = case_when(waterQuality<1~"inf1",waterQuality>=1 & waterQuality<2~"sup1_inf2",waterQuality>=2 & waterQuality<3~"sup2_inf3",waterQuality>=3~"sup3")) %>% ggplot(aes(firstDayMonth,fill=waterCat)) + geom_bar()
dev.off()

## Second plot is with proportion of monthly values
pdf("Plots/waterQuality_category_monthly_proportion.pdf")
train_small %>% mutate(waterCat = case_when(waterQuality<1~"inf1",waterQuality>=1 & waterQuality<2~"sup1_inf2",waterQuality>=2 & waterQuality<3~"sup2_inf3",waterQuality>=3~"sup3")) %>% ggplot(aes(firstDayMonth,fill=waterCat)) + geom_bar(position="fill")
dev.off()


## Plot the histograms of all indicators
## We order the names in the process
namesIndicators = subset(names(train_small),"indicator")
namesIndicators = namesIndicators[order(as.numeric(str_extract(namesIndicators,"(?<=indicator\\.)\\d+")))]

pdf(paste0("Plots/Indicators_histograms.pdf"))
for(i in namesIndicators)
{
	print(ggplot() + geom_histogram(aes(na.omit(train[[i]]))) + xlab(i))
}
dev.off()

## Plot the water quality against each indicator
pdf(paste0("Plots/WaterQuality_vs_indicators.pdf"))
for(i in namesIndicators)
{
	print(ggplot() +
	geom_point(aes(x=train_small[[i]],y=train_small[["waterQuality"]]),alpha=0.2,col=4) + xlab(i) + ylab("Water quality") + geom_smooth(aes(x=train_small[[i]],y=train_small[["waterQuality"]]),col=2))
}
dev.off()

## PCA exploration

pca1 = dudi.pca(train_small,nf=3,scann=F)
pdf("Plots/Corcircle_PCA_Axis1-2.pdf")
s.corcircle(pca1$co, cgrid = 0, full = FALSE, clab = 0.6, xax = 1, yax = 2)
dev.off()
pdf("Plots/Corcircle_PCA_Axis1-3.pdf")
s.corcircle(pca1$co, cgrid = 0, full = FALSE, clab = 0.6, xax = 1, yax = 3)
dev.off()
pdf("Plots/Corcircle_PCA_Axis2-3.pdf")
s.corcircle(pca1$co, cgrid = 0, full = FALSE, clab = 0.6, xax = 2, yax = 3)
dev.off()


## Caret training of randomForest
set.seed(28)
fmla = as.formula("waterQuality~.")
grid = data.frame(.mtry = seq(8,20,2))
control = trainControl(method="cv", number=10)
rffit = train(fmla, data = na.omit(train_small), method = 'rf',
trControl = control, importance = T, tuneGrid = grid)

## Best model is for mtry=12
rf1 = randomForest(fmla, train_small, mtry = 12)

## Look at variable importance
varImpPlot(rf1)

## Order indicators according to their importance in rf1
orderedIndicators = row.names(rf1$importance)[order(as.numeric(rf1$importance),decreasing=T)]


## Transform irregular time series into a regular one every 10 minutes
## waterQuality_zoo = zoo(train$waterQuality,order.by=train$time)
## Grid every 10 minutes (600 seconds)
## g = seq(start(waterQuality_zoo),end(waterQuality_zoo),600)
## waterQuality_zoo = merge(g,waterQuality_zoo)


## Calculate error in cross-validation 70/30 done 10 times
n = nrow(train_small)
niter = 10
p = 0.7
err = NA
for(i in 1:niter)
{
  
	id = sample(1:n,size=round(p*n))
	rftemp = randomForest(waterQuality~., train_small[id,])
	pred = predict(rftemp,select(train_small,-waterQuality)[-id,])
	err = c(err,rmse(train_small_simple$Output.Parameter[-id],pred))
}
err = err[-1]
hist(err)

## Do a series of randomForests by each time removing the n least
   important indicators (1<=n<=18) and calculating the median rmse over 10 iterations

err2 = NA
t = Sys.time()

for(i in c(1,10,12,14,18))
{
  fmlatemp = as.formula(paste0("waterQuality~",paste0(head(orderedIndicators,-i),collapse="+")))
  errtemp = NA
  for(j in 1:10)
  {
    id = sample(1:n,size=round(p*n))
    rftemp = randomForest(fmlatemp, train_small[id,])
    predtemp = predict(rftemp,select(train_small,-waterQuality)[-id,])
    errtemp = c(errtemp,rmse(train_small$Output.Parameter[-id],predtemp))    
  }
  err2 = c(err2,median(errtemp,na.rm=T))
}

err2 = err2[-1]

Sys.time()-t

	    


      

      

