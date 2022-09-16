## Step 0: Please install all the packages listed below prior to run the R code
```r
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("caret")
install.packages("catch")
install.packages("patchwork")
install.packages("randomForest")
```

## Step 1: Load the libraries into your R
```r
library("pheatmap")
library("ggplot2")
library("caret")
library("catch")
library("patchwork")
library("randomForest")
```

## Step 2: Source the R script
```r
source("code/catch_demo.R")
```

## Step 3: Load the datasets
```r
tensor=read.csv(file = "data/tensor.csv",check.names = FALSE)
covariate=read.csv(file = "data/covariate.csv")
label=read.csv(file = "data/label.csv", header = T)
pathway=read.csv(file = "data/pathway.csv", row.names = 1,check.names = FALSE)
```

## Step 4: Perform data preprocessing
```r
ppo = preprocessing(x = tensor, y = label$MSIStatus, z = covariate,
                    alpha_group=pathway,testx = tensor, testy = label$MSIStatus, testz = covariate, rescale_x = "standard", baseline.factor = "MSS")
```

## Step 5: Run the catch model
```r
model = catch_model(ppo=ppo)
```

## Step 6: Show the figures and tables
```r
plot_sorted_complicated_heatmap(label, model$heatmap_matrix2) # figure 2
print(plot_boxplot(model$estimate_value, ppo, model$heatmap_matrix)) # figure 3
plot_alpha_heatmap(model$alpha_matrix, pathway) # figure 4
data.frame(Beta=sort(round(model$estimate_value, 2))) # table 1
# It may take a few minutes to finish the task below, so please be patient
suppressWarnings(cvprocess(tensor, label, covariate, ppo)) # table 2
```
