rescale_function = function(x,f){
  xdim = dim(x[[1]])
  rescale_matrix = matrix(NA,ncol = length(x[[1]]),nrow = length(x))
  for (i in 1:length(x)) {
    rescale_matrix[i,] = as.matrix(x[[i]],nrow = 1)
  }
  rescale_matrix = f(rescale_matrix)
  for (i in 1:length(x)) {
    x[[i]] = array(rescale_matrix[i,],dim=xdim)
  }
  return(x)
}
rescale_standard = function(x){
  for (i in 1:length(x[1,])) {
    mean_x = mean(x[,i])
    sd_x = sd(x[,i])
    if(sd_x==0){
      x[,i] = x[,i]-mean_x
    }else{
      x[,i] = (x[,i]-mean_x)/(sd_x)
    }
    
  }
  
  return(x)
}
t_adjten = function (x, z, y, testx = NULL, testz = NULL, is.centered = FALSE)
{
  n = length(x)
  dimen = dim(x[[1]])
  nvars = prod(dimen)
  vec_x = matrix(0, nvars, n)
  for (i in 1:n) {
    vec_x[, i] = matrix(x[[i]], nvars, 1)
  }
  vec_x = t(vec_x)
  if (!is.null(testx)) {
    tesize = length(testx)
    vec_testx = matrix(0, nvars, tesize)
    for (i in 1:tesize) {
      vec_testx[, i] = matrix(testx[[i]], nvars, 1)
    }
    vec_testx = t(vec_testx)
  }
  z = as.matrix(z)
  cz = z
  n = dim(z)[1]
  q = dim(z)[2]
  cvecx = vec_x
  nclass <- as.integer(length(unique(y)))
  if (is.centered == FALSE) {
    for (i in 1:nclass) {
      if (q > 1) {
        cz[y == i, ] = sweep(z[y == i, ], 2, colMeans(z[y == 
                                                          i, ]))
      }
      else {
        cz[y == i] = z[y == i] - mean(z[y == i])
      }
      cvecx[y == i, ] = sweep(vec_x[y == i, ], 2, colMeans(vec_x[y == i, ]))
    }
  }
  c = solve(t(cz) %*% cz) %*% t(cz)
  c = c %*% cvecx
  vec_xres = vec_x - z %*% c
  xres = array(list(), n)
  for (i in 1:n) {
    xres[[i]] = array(vec_xres[i, ], dim = dimen)
  }
  if (!is.null(testx)) {
    vec_testxres = vec_testx - testz %*% c
    testxres = array(list(), tesize)
    for (i in 1:tesize) {
      testxres[[i]] = array(vec_testxres[i, ], dim = dimen)
    }
  }
  else {
    testxres = NULL
  }
  muz = matrix(0, nrow = q, ncol = (nclass - 1))
  for (i in 2:nclass) {
    if (q > 1) {
      muz[, (i - 1)] = apply(z[y == i, ], 2, mean) - apply(z[y == 
                                                               1, ], 2, mean)
    }
    else {
      muz[i - 1] = mean(z[y == i]) - mean(z[y == 1])
    }
  }
  gamma = solve(cov(z)) %*% muz
  outlist = list(gamma = gamma, xres = xres, testxres = testxres,alpha = c,vec_xres=vec_xres)
  outlist
}
predict.catch.prob = function (object, newx, z = NULL, ztest = NULL, gamma = NULL) 
{
  if (is.null(gamma)) {
    pred = predict.tsda(object, newx)
  }
  else {
    thetatm <- object$beta
    theta = array(list(), length(thetatm))
    mu <- object$mu
    prior <- object$prior
    nclass <- length(prior)
    dimen <- dim(newx[[1]])
    nvars <- prod(dimen)
    nlambda <- length(theta)
    gamma = as.matrix(gamma)
    q = dim(gamma)[1]
    z = as.matrix(z)
    ztest = as.matrix(ztest)
    for (i in 1:nlambda) {
      theta[[i]] = matrix(0, nrow = dim(thetatm[[i]])[1] + 
                            q, ncol = nclass - 1)
      for (j in 1:nclass - 1) {
        theta[[i]][1:nvars, j] = matrix(thetatm[[i]][, 
                                                     j], ncol = 1)
        for (qq in 1:q) {
          theta[[i]][nvars + qq, j] = gamma[qq, j]
        }
      }
    }
    mubar = matrix(list(), nclass - 1, 1)
    for (i in 1:(nclass - 1)) {
      mubar[[i]] = (mu[[i + 1]] + mu[[1]])/2
    }
    n <- length(newx)
    nn <- length(object$x)
    x.train <- object$x
    vecx.train = matrix(0, ncol = nn, nrow = nvars + q)
    vecnewx = matrix(0, ncol = n, nrow = nvars + q)
    for (i in 1:nn) {
      vecx.train[1:nvars, i] <- matrix(x.train[[i]], ncol = 1)
      for (qq in 1:q) {
        vecx.train[nvars + qq, i] = z[i, qq]
      }
    }
    vecx.train = t(vecx.train)
    for (i in 1:(length(newx))) {
      vecnewx[1:nvars, i] <- matrix(newx[[i]], ncol = 1)
      for (qq in 1:q) {
        vecnewx[nvars + qq, i] = ztest[i, qq]
      }
    }
    vecnewx = t(vecnewx)
    y.train <- object$y
    pred <- matrix(0, n, nlambda)
    pred[1] <- which.max(prior)
    for (i in 1:nlambda) {
      nz <- sum(theta[[i]][, 1] != 0)
      if (nz == 0) {
        pred[, i] <- which.max(prior)
      }
      else {
        xfit <- vecx.train %*% theta[[i]][, 1:(min(nclass - 
                                                     1, nz)), drop = FALSE]
        xfit.sd <- matrix(0, nclass, ncol(xfit))
        for (j in 1:nclass) {
          xfit.sd[j, ] <- apply(xfit[y.train == j, , 
                                     drop = FALSE], 2, sd)
        }
        xfit.sd <- apply(xfit.sd, 2, min)
        if (min(xfit.sd) < 1e-04) {
          print("xfit.sd less then 1e-04,can't output probability")
          pred = NULL
        }
        else {
          l <- lda(xfit, y.train)
          pred <- predict(l, vecnewx %*% theta[[i]][, 
                                                         1:(min(nclass - 1, nz))])$posterior
        }
      }
    }
  }
  pred
}
preprocessing = function(x,y,z=NULL,alpha_group,testx=NULL,testy=NULL,testz=NULL,
                         baseline.factor=NULL,rescale_x=NULL,rescale_z=NULL,
                         tensor_namelist=NULL,cov_matrix_namelist=NULL){
  x_train_as_test = is.null(testx)
  y_train_as_test = is.null(testy)
  z_train_as_test = is.null(testz)
  train_as_test=FALSE
  if(is.null(z)){
    if(x_train_as_test&y_train_as_test){
      train_as_test = TRUE
    }else if(!((!x_train_as_test)&(!y_train_as_test))){
      stop("either testx and testy must be exist or not exist")
    }
  }else{
    if(x_train_as_test&y_train_as_test&z_train_as_test){
      train_as_test = TRUE
    }else if(!((!x_train_as_test)&(!y_train_as_test)&(!z_train_as_test))){
      stop("either testx and testy must be exist or not exist")
    }
  }
  if(is.numeric(y)){
    orig_y = sort(unique(y))
  }else{
    orig_y = unique(y)
  }
  if(!is.null(baseline.factor)){
    if(baseline.factor %in% orig_y){
      orig_y = orig_y[orig_y!=baseline.factor]
      orig_y = c(baseline.factor,orig_y)
    }else{
      warning("baseline.factor not in y,ignore baseline.factor.")
    }
  }
  trans_y = c(1:length(orig_y))
  new_y = rep(NA,length(y))
  for (i in 1:length(orig_y)) {
    new_y[y==orig_y[i]] = trans_y[i]
  }
  #transform testy.
  if(!is.null(testy)){
    new_test_y = rep(NA,length(testy))
    for (i in 1:length(orig_y)) {
      new_test_y[testy==orig_y[i]] = trans_y[i]
    }
  }else if(train_as_test){
    testy = y
    new_test_y = new_y
  }
  if(class(x)=="data.frame"){
    y_n = length(new_y)
    actual_x_length = length(colnames(x))
    matrix_n = ceiling(length(x[1,])^(1/2))
    model_x_length = matrix_n*matrix_n
    x_dim = c(matrix_n,matrix_n)
    while(length(x[1,]) < (matrix_n*matrix_n)){x = cbind(x,0)}
    tensor.array = array(list(),y_n)
    for (i in 1:y_n){
      tensor.array[[i]] =  array(as.numeric(x[i,]),dim=x_dim)
    }
    if(is.null(tensor_namelist)){
      tensor_namelist = colnames(x)
    }else if (length(tensor_namelist) != colnames(x)){
      tensor_namelist = colnames(x)
      warning("length of tensor_namelist not equal to length of x ,ignore tensor_namelist.")
    }
    name_array = NULL
    if(train_as_test){
      test.tensor.array = tensor.array
    }else{
      test_y_n = length(new_test_y)
      if(class(testx) != "data.frame"){
        stop("x is dataframe but testx isn't a dataframe.")
      }else{
        while(length(testx[1,]) < (matrix_n*matrix_n)){
          testx = cbind(testx,0)
        }
        test.tensor.array = array(list(),test_y_n)
        for (i in 1:test_y_n){
          test.tensor.array[[i]] =  array(as.numeric(testx[i,]),dim=x_dim)
        }
      }
    }
    
  }else{
    y_n = length(new_y)
    actual_x_length = length(x[[1]])
    model_x_length = length(x[[1]])
    x_dim = dim(x[[1]])
    tensor.array = x
    if(train_as_test){
      test.tensor.array = tensor.array
    }else{
      test.tensor.array = testx
    }
    if(is.null(tensor_namelist)){
      tensor_namelist = paste0("tensor-",1:actual_x_length)
    }else if(length(tensor_namelist) != actual_x_length){
      tensor_namelist = paste0("tensor-",1:actual_x_length)
      warning("length of tensor_namelist not equal to length of x ,ignore tensor_namelist.")
    }
    name_array = array(tensor_namelist,dim=x_dim)
  }
  
  if(is.null(rescale_x)){
  }else if(class(rescale_x)=="character"){
    if(rescale_x=="standard"){
      
      tensor.array = rescale_function(tensor.array,rescale_standard)
      if(train_as_test){
        test.tensor.array = tensor.array
      }else{
        test.tensor.array = rescale_function(test.tensor.array,rescale_standard)
      }
      
    }else if(rescale_x=="01"){
      
      tensor.array = rescale_function(tensor.array,rescale_01)
      if(train_as_test){
        test.tensor.array = tensor.array
      }else{
        test.tensor.array = rescale_function(test.tensor.array,rescale_01)
      }
  }else{stop("rescale_x must be NULL \"01\" \"standard\" or a rescale function.")}
  
    
  }else if(class(rescale_x)=="function"){
    
    tensor.array = rescale_function(tensor.array,rescale_x)
    if(train_as_test){
      test.tensor.array = tensor.array
    }else{
      test.tensor.array = rescale_function(test.tensor.array,rescale_x)
    }
    
  }else{stop("rescale_x must be NULL , \"01\" , \"standard\" or a rescale function.")}
  
  if((!is.null(z))&is.null(cov_matrix_namelist)){
    if(is.null(colnames(z))){
      cov_matrix_namelist = paste0("cov-",1:length(z[1,]))
    }else{cov_matrix_namelist = colnames(z)}
  }else if(length(cov_matrix_namelist) != length(z[1,])){
    warning("length of cov_matrix_namelist can't fit with matrix z,ignore cov_matrix_namelist.")
    if(is.null(colnames(z))){
      cov_matrix_namelist = paste0("cov-",1:length(z[1,]))
    }else{cov_matrix_namelist = colnames(z)}
  }
  if(!is.null(z)){
    z = as.matrix(z)
  }
  if(!is.null(testz)){
    testz = as.matrix(testz)
  }
  if(train_as_test){
    testz = z
  }
  if(is.null(rescale_z)){
  }else if(class(rescale_z)=="character"){
    if(rescale_z=="standard"){
      
      z = rescale_standard(z)
      if(train_as_test){
        testz = z
      }else{
        testz = rescale_standard(testz)
      }
      
    }else if(rescale_z=="01"){
      
      z = rescale_01(z)
      if(train_as_test){
        testz = z
      }else{
        testz = rescale_01(testz)
      }
      
    }else{stop("rescale_z must be NULL , \"01\" , \"standard\" or a rescale function.")}
  }else if(class(rescale_z)=="function"){
    
    z = rescale_z(z)
    if(train_as_test){
      testz = z
    }else{
      testz = rescale_z(testz)
    }
    
  }else{stop("rescale_z must be NULL , \"01\" , \"standard\" or a rescale function.")}
  
  
  alpha_group$`Metabolic Pathway`=factor(alpha_group$`Metabolic Pathway`)
  output  = list()
  output$alpha_group=alpha_group
  output$actual_x_length = actual_x_length
  output$x = x
  output$tensor.array = tensor.array
  output$testx = testx
  output$test.tensor.array = test.tensor.array
  
  output$trans_y = trans_y
  output$orig_y = orig_y
  output$y = y
  output$new_y = new_y
  output$testy = testy
  
  output$z = z
  output$testz = testz
  
  output$name_array = name_array
  output$tensor_namelist = tensor_namelist
  output$cov_matrix_namelist = cov_matrix_namelist
  return(output)
}


catch_model = function(ppo,lambda=NULL,max_p=15,min_p=5,positive.factor=NULL,cv.seed=1,annotation_colors = NA,alpha_show_name=F,min_abs_alpha=0.1){
  positive.factor="MSI"
  alpha_group=ppo$alpha_group
  if(max_p<min_p){
    stop("max_p should't smaller then min_p")
  }
  if(!("catch" %in% (.packages()))){
    library("catch")
  }
  if(!("caret" %in% (.packages()))){
    library("caret")
  }
  if(!("pheatmap" %in% (.packages()))){
    library("pheatmap")
  }
  if(!("MASS" %in% (.packages()))){
    library("MASS")
  }
  in_x = c(rep(TRUE,ppo$actual_x_length),rep(FALSE,(length(ppo$tensor.array[[1]])-ppo$actual_x_length)))
  y_n = length(ppo$new_y)
  test_y_n = length(ppo$testy)
  if(!is.null(positive.factor)){
    if(!(positive.factor %in% ppo$orig_y)){
      positive.factor = NULL
      warning("positive.factor not in y,ignore baseline.factor.")
    }
  }
  
  set.seed(cv.seed)
  cv_model = cv.catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y)
  estimate_mcve = cv_model$cvm
  cv_model2 = catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y,
                    testx=ppo$tensor.array,testz=ppo$z)
  estimate_significant_n = cv_model2$df
  estimate_lambda = cv_model2$lambda
  estimate_y = cv_model2$pred
  estimate_precision = c()
  for(i in 1:length(estimate_y[1,])){
    pre = sum(estimate_y[,i]==ppo$new_y)/length(estimate_y[,i])
    estimate_precision = c(estimate_precision,pre)
  }
  if(sum(estimate_significant_n<max_p)==0){
    max_p = min(estimate_significant_n)
    warning("max_p is smaller then all of cross validation estimate,ignore it information.")
  }
  if(sum(estimate_significant_n>min_p)==0){
    min_p = max(estimate_significant_n)
    warning("min_p is larger then all of cross validation estimate,ignore it information.")
  }  
  if(!is.null(lambda)){
    best.lambda = lambda
  }else{
    
    bool1 = estimate_significant_n<=max_p
    bool2 = estimate_significant_n>=min_p
    bool3 = bool1&bool2
    
    
    if(sum(bool3)>0){
      best.lambda = estimate_lambda[(estimate_mcve == min(estimate_mcve[bool3]))&bool3]
      best.lambda = best.lambda[length(best.lambda)]
    }else{
      up_find = FALSE
      down_find = FALSE
      while(!up_find){
        max_p = max_p+1
        if (max_p %in% estimate_significant_n) {
          up_lambda = estimate_lambda[estimate_significant_n==max_p]
          up_find = TRUE
        }
      }
      while(!down_find){
        min_p = min_p-1
        if (min_p %in% estimate_significant_n) {
          down_lambda = estimate_lambda[estimate_significant_n==min_p]
          down_find = TRUE
        }
      }
      best.lambda = (up_lambda+down_lambda)/2
    }
    
  }
  if(!is.null(ppo$z)){
    adj_model = t_adjten(x = ppo$tensor.array,z=ppo$z,y=ppo$new_y,testx=ppo$test.tensor.array,testz=ppo$testz)
    alpha = as.matrix(adj_model$alpha)
    colnames(alpha) = ppo$tensor_namelist
    row.names(alpha) = ppo$cov_matrix_namelist
    gamma = adj_model$gamma
    
    testxres = adj_model$testxres
    adj = adj_model$xres
    if(class(ppo$x)=="data.frame"){
      matrix_n = ceiling(length(ppo$x[1,])^(1/2))
      adj.matrix = data.frame(array(dim=c(length(ppo$testy),matrix_n*matrix_n)))
      colnames(adj.matrix) = ppo$tensor_namelist
      for (i in 1:test_y_n) {
        adj.matrix[i,] = array(testxres[[i]])
      }
      adj.matrix = adj.matrix[,1:ppo$actual_x_length]
    }else{adj.matrix = testxres}
  }
  model = catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y ,
                lambda=best.lambda,testx = ppo$test.tensor.array,testz = ppo$testz)
  
  model_pred = model$pred[,1]
  model_pred_prob = predict.catch.prob(model,adj_model$testxres,ppo$z,ppo$testz,adj_model$gamma)
  model_pred2 = rep(NA,length(model_pred))
  for (i in 1:length(ppo$orig_y)) {
    model_pred2[model_pred==ppo$trans_y[i]] = ppo$orig_y[i]
  }
  if(is.null(positive.factor)){
    pred_table = confusionMatrix(factor(model_pred2),factor(ppo$testy),positive=levels(factor(ppo$testy))[length(levels(factor(ppo$testy)))])
  }else{
    pred_table = confusionMatrix(factor(model_pred2),factor(ppo$testy),positive=positive.factor)
  }
  
  if(length(unique(ppo$y))==2){
    pred_table = list(pred_table,pred_table$byClass[7])
  }
  beta = model$beta$`1`
  if(is.null(beta)){
    beta = model$beta[[1]]
  }
  rownames(beta) = ppo$tensor_namelist
  significant = rep(FALSE,length(beta[,1]))
  for (i in 1:(length(ppo$orig_y)-1)) {
    significant = significant|(beta[,i] !=0)
  }
  

  output=list()

  if(!is.null(ppo$z)){
    output$adj.array = adj.matrix
  }
  if(sum(significant)!=0){
    output$estimate_value = beta[significant&in_x,]
  }
  output$pred = model_pred2
  output$pred_table = pred_table
  # output$precision_plot = p1
  # output$mcve_plot = p2
  if(sum(significant)>1){
    heatmap_matrix = array(dim=c(test_y_n,length(testxres[[1]])))
    for (i in 1:test_y_n) {
      heatmap_matrix[i,] = testxres[[i]][1:length(testxres[[1]])]
    }
    colnames(heatmap_matrix) = ppo$tensor_namelist[1:length(ppo$test.tensor.array[[1]])]
    row.names(heatmap_matrix) = 1:length(heatmap_matrix[,1])
    heatmap_beta = beta[significant&in_x,1]
    heatmap_matrix = heatmap_matrix[,significant&in_x]
    heatmap_matrix = heatmap_matrix[,order(heatmap_beta,decreasing=TRUE)]
    #========
    prob = matrix(ncol=length(ppo$trans_y),nrow = length(ppo$test.tensor.array))
    prob[,1] = 0
    phi = matrix(ncol=length(ppo$trans_y),nrow = length(ppo$z[1,]))
    mu = matrix(ncol=length(ppo$trans_y),nrow = length(matrix(model$mu[1,1][[1]],ncol=1)))
    for (i in 1:length(ppo$trans_y)) {
      if(!is.null(ppo$z)){
        phi[,i] = colMeans(x=ppo$z[ppo$new_y==i,])
        mu[,i] = matrix(model$mu[i,1][[1]],ncol=1)        
      }else{mu[,i] = matrix(model$mu[i,1][[1]],ncol=1)  }

    }
    f_alpha = c()
    
    for (i in 2:length(ppo$trans_y)) {
      if(!is.null(ppo$testz)){
        f_alpha = c(f_alpha,log(model$prior[i]/model$prior[1]) - ((phi[,i]+phi[,1]) %*% gamma[,(i-1)])/2 - sum((mu[,i]+mu[,1])/2*beta[,i-1]))
      }else{
        f_alpha = c(f_alpha,log(model$prior[i]/model$prior[1])-sum((mu[,i]+mu[,1])/2*beta[,i-1]))
        gamma = matrix(0,nrow=i-1)
      }
    }
    
    for (i in 2:length(ppo$trans_y)) {
      for (j in 1:length(ppo$test.tensor.array)) {
        if(!is.null(ppo$z)){
          prob[j,i] = f_alpha[i-1]+sum(ppo$testz[j,]*gamma[,i-1])+sum(matrix(testxres[[j]],ncol=1)*beta[,i-1])
        }else{prob[j,i] = f_alpha[i-1]+sum(matrix(ppo$test.tensor.array[[j]],ncol=1)*beta[,i-1])}
        
      }
    }
    if(length(ppo$trans_y)==2){
      p0 = ppo$testy==ppo$orig_y[1]
      p1 = ppo$testy==ppo$orig_y[2]
      mt0 = heatmap_matrix[p0,,drop=F]
      mt0 = mt0[order(model_pred_prob[p0,2]),,drop=F]
      mt1 = heatmap_matrix[p1,,drop=F]
      mt1 = mt1[order(model_pred_prob[p1,2]),,drop=F]
      heatmap_matrix = rbind(mt0,mt1)
    }else if (length(ppo$trans_y)==3){
      
      pp0 = model_pred2==ppo$orig_y[1]
      pp1 = model_pred2==ppo$orig_y[2]
      pp2 = model_pred2==ppo$orig_y[3]
      p0 = ppo$testy==ppo$orig_y[1]
      p1 = ppo$testy==ppo$orig_y[2]
      p2 = ppo$testy==ppo$orig_y[3]
      
      mt0_0 = heatmap_matrix[p0&pp0,,drop=F]
      mt0_0 = mt0_0[order(model_pred_prob[p0&pp0,1],decreasing=T),,drop=F]
      mt0_12 = heatmap_matrix[p0&(!pp0),,drop=F]
      mt0_12 = mt0_12[order(model_pred_prob[p0&(!pp0),1],decreasing=T),,drop=F]
      
      mt1_0 = heatmap_matrix[p1&pp0,,drop=F]
      mt1_0 = mt1_0[order(model_pred_prob[p1&pp0,2]),,drop=F]
      mt1_1 = heatmap_matrix[p1&pp1,,drop=F]
      mt1_1 = mt1_1[order(model_pred_prob[p1&pp1,2]),,drop=F]
      mt1_2 = heatmap_matrix[p1&pp2,,drop=F]
      mt1_2 = mt1_2[order(model_pred_prob[p1&pp2,2]),,drop=F]
      
      mt2_01 = heatmap_matrix[p2&(!pp2),,drop=F]
      mt2_01 = mt2_01[order(model_pred_prob[p2&(!pp2),3]),,drop=F]
      mt2_2 = heatmap_matrix[p2&pp2,,drop=F]
      mt2_2 = mt2_2[order(model_pred_prob[p2&pp2,3]),,drop=F]
      
      heatmap_matrix = rbind(mt0_0,mt0_12,mt1_0,mt1_1,mt1_2,mt2_01,mt2_2)
    }
    
    output$prob = prob
    output$prob2 = model_pred_prob
    output$heatmap_matrix = heatmap_matrix
    
    for (i in 1:length(heatmap_matrix[1,])) {
      if(sd(heatmap_matrix[,i])==0){
        heatmap_matrix[,i] = rep(0,length(heatmap_matrix[,i]))
      }else{
        heatmap_matrix[,i] = (heatmap_matrix[,i]-mean(heatmap_matrix[,i]))/sd(heatmap_matrix[,i])
      }
      
    }
    output$heatmap_matrix2 = heatmap_matrix
    


    if(!is.null(ppo$z)){
      heatmap_matrix2 = alpha[,significant&in_x]
      heatmap_matrix2 = heatmap_matrix2[,order(heatmap_beta,decreasing=TRUE)]
      alpha_group=ppo$alpha_group
      if(!is.null(alpha_group)){
        mt = heatmap_matrix2[alpha_group[,1]==unique(alpha_group[,1])[1],]
        bool = rowSums(abs(mt)>min_abs_alpha)>0
        if((sum(bool)!=0)&(sum(!bool)!=0)){
          mtT = mt[bool,,drop=F]
          mtF = mt[!bool,,drop=F]
          mt = rbind(mtF,mtT)
        }

        for (i in 2:length(unique(alpha_group[,1]))) {
          mt0 = heatmap_matrix2[alpha_group[,1]==unique(alpha_group[,1])[i],]
          bool = rowSums(abs(mt0)>min_abs_alpha)>0
          if((sum(bool)!=0)&(sum(!bool)!=0)){
          mtT = mt0[bool,,drop=F]
          mtF = mt0[!bool,,drop=F]
          mt0 = rbind(mtF,mtT)
          }
          mt = rbind(mt,mt0)
        }
        heatmap_matrix2 = mt
      }
    output$alpha_matrix = t(heatmap_matrix2)
    #output$alpha_heatmap = pheatmap(t(heatmap_matrix2),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=alpha_group,show_colnames = alpha_show_name)
    }
  }
  else{
    warning("Less then 2 significant tensor,no estimate_value???tensor_heatmap???alpha_heatmap")
    }
  return(output)
}
plot_alpha_heatmap=function(heatmap_matrix2,alpha_group){
  heatmap_matrix2=heatmap_matrix2[rev(rownames(heatmap_matrix2)),]
  return(pheatmap(heatmap_matrix2,cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=alpha_group,show_colnames = F))
}
plot_sorted_complicated_heatmap=function(label2,heatmap_matrix){
  heatmap_matrix=heatmap_matrix[,rev(colnames(heatmap_matrix))]
  col.matrix = label2[,c(2,3,4,5)]
  row.names(col.matrix) = seq(1,300)#
  colnames(col.matrix)[1] = "MS"
  for (i in 1:4) {
    col.matrix[,i] = factor(col.matrix[,i])
  }
  ann_colors = list(
    TP53 = c("M"="#d41515", "W"="#bcc213"),
    APC = c("M"="#7d0c0c", "W"="#7d7d0c")
  )
  a = heatmap_matrix
  col.matrix = col.matrix[row.names(heatmap_matrix),]
  a = a[order(col.matrix[,4]),]
  col.matrix = col.matrix[order(col.matrix[,4]),]
  a = a[order(col.matrix[,3]),]
  col.matrix = col.matrix[order(col.matrix[,3]),]
  a = a[order(col.matrix[,2]),]
  col.matrix = col.matrix[order(col.matrix[,2]),]
  a = a[order(col.matrix[,1],decreasing = FALSE),]
  col.matrix = col.matrix[order(col.matrix[,1]),]
  return(pheatmap(t(a),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=col.matrix,show_colnames=F,annotation_colors=ann_colors))
}### -> put in source.R
plot_complicated_heatmap=function(label2,heatmap_matrix){
  col.matrix = label2[,c(2,3,4,5)]
  row.names(col.matrix) = seq(1,300)#
  colnames(col.matrix)[1] = "MS"
  for (i in 1:4) {
    col.matrix[,i] = factor(col.matrix[,i])
  }
  ann_colors = list(
    TP53 = c("0"="#d41515", "1"="#bcc213"),
    APC = c("0"="#7d0c0c", "1"="#7d7d0c")
  )
  b = heatmap_matrix
  col.matrix = col.matrix[row.names(heatmap_matrix),]
  b = b[order(col.matrix[,1],decreasing = FALSE),]
  col.matrix = col.matrix[order(col.matrix[,1]),]
  return(pheatmap(t(b),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=col.matrix,show_colnames=F,annotation_colors=ann_colors))
}### -> put in source.R
plot_heatmap=function(label2,heatmap_matrix){
  ann_colors = list(
    TP53 = c("0"="#d41515", "1"="#bcc213"),
    APC = c("0"="#7d0c0c", "1"="#7d7d0c")
  )
  col.matrix = as.data.frame(label2[,c(2)],byrow=TRUE)
  row.names(col.matrix) = seq(1,300)
  col.matrix[,1] = factor(col.matrix[,1])
  a = heatmap_matrix
  col.matrix = col.matrix[row.names(heatmap_matrix),,drop=FALSE]
  a = a[order(col.matrix[,1],decreasing = FALSE),]
  col.matrix = col.matrix[order(col.matrix[,1],decreasing=FALSE),]
  col.matrix=as.data.frame(col.matrix,byrow=TRUE)
  colnames(col.matrix)[1] = "y_reference"
  return(pheatmap(t(a),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=col.matrix,show_colnames=F,annotation_colors=ann_colors))
}### -> put in source.R
plot_boxplot=function(est_value,ppo2,heatmap_matrix){
  n = 75
  outputbox=list()
  for (i in row.names(data.frame(est_value))) {
    set.seed(1)
    target = ppo2$x[,colnames(ppo2$x)==i]
    p.value = round(wilcox.test(target[1:n],target[(1+n):(n*4)])$p.value,3)
    j=i
    j=gsub("X","",j)
    j=gsub("\\.","-",j)
    if(j=="C18-0-LPE"){j="LPE,C18:0"}
    if(j=="C36-1-PC"){j="PC,C36:1"}
    if(j=="C14-0-CE"){j="CE,C14:0"}
    if(p.value==0){p.value="<0.001"}
    group1 = factor(c(rep("MSI",n),rep("MSS",n*3)))
    p1 = ggplot(data.frame(target), aes(x = group1 ,y = target ,fill  = group1))+
      geom_boxplot(outlier.alpha=0)+
      geom_jitter(color="black", size=0.1, alpha=0.8)+
      scale_x_discrete("")+scale_y_continuous(j)+ labs(title = "(Non-adjusted)", subtitle = paste0(expression("p")," = ",p.value))+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
    target2 = (target-mean(target))/sd(target)
    p.value2 = round(wilcox.test(target2[1:n],target2[(1+n):(n*4)])$p.value,3)
    if(p.value2==0){p.value2="<0.001"}
    p2 = ggplot(data.frame(target2), aes(x = group1 ,y = target2 ,fill  = group1))+
      geom_boxplot(outlier.alpha=0)+
      geom_jitter(color="black", size=0.1, alpha=0.8)+
      scale_x_discrete("")+scale_y_continuous(j)+ labs(title = "(Standardized)", subtitle = paste0(expression("p")," = ",p.value2))+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
    
    target3 = heatmap_matrix[,colnames(heatmap_matrix)==i]
    p.value3 = round(wilcox.test(target3[1:(n*3)],target3[(1+n*3):(n*4)])$p.value,3)
    if(p.value3==0){p.value3="<0.001"}
    group = factor(c(rep("MSS",n*3),rep("MSI",n)))
    p3 = ggplot(data.frame(target3), aes(x = group ,y = target3 ,fill  = group))+
      geom_boxplot(outlier.alpha=0)+
      geom_jitter(color="black", size=0.1, alpha=0.8)+
      scale_x_discrete("")+scale_y_continuous(j)+ labs(title = "(Catch-adjusted)", subtitle = paste0(expression("p")," = ",p.value3))+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
    
    p4 = ggplot()+theme(panel.background=element_rect(fill='#ffffff', color="#ffffff"))
    ggp_all <- p1 + p2+p3
    outputbox[[i]]=ggp_all
  }
  return(outputbox)
}
cvprocess=function(tensor2,label2,covariate2,ppo2){
  library("randomForest")
  total_size = 300
  MSI_size = 75
  groupn = 9
  MIS_multi = 3
  test_Accuracy = c()
  test_Sensitivity = c()
  test_Specificity = c()
  test_Precision = c()
  test_recall = c()
  test_F1 = c()
  for (i in 1:100) {

    set.seed(i*10)
    rs1 = sample(1:MSI_size,groupn,replace=F)
    set.seed(i*10)
    rs2 = sample((MSI_size+1):total_size,groupn*MIS_multi,replace=F)
    rs = c(rs1,rs2)
    train_tensor = tensor2[-rs,]
    valid_tensor = tensor2[rs,]
    train_label = label2$MSIStatus[-rs]
    valid_label = label2$MSIStatus[rs]
    train_cov_matrix = as.matrix(covariate2[-rs,])
    valid_cov_matrix = as.matrix(covariate2[rs,])

    ppo2_2 = preprocessing(x=train_tensor,y=train_label,z=train_cov_matrix,
                           alpha_group = ppo2$alpha_group,testx=valid_tensor,testy=valid_label,testz=valid_cov_matrix,rescale_x = "standard",baseline.factor = "NOT MSI")
    invisible(capture.output(model2_2 <- catch_model(ppo=ppo2_2,max_p=15,min_p=5,cv.seed=1,positive.factor = "MSI")))
    test_Sensitivity = c(test_Sensitivity,model2_2$pred_table[[1]]$byClass[1])
    test_Specificity = c(test_Specificity,model2_2$pred_table[[1]]$byClass[2])
    test_Accuracy = c(test_Accuracy,model2_2$pred_table[[1]]$overall[1])
    test_Precision = c(test_Precision,model2_2$pred_table[[1]]$byClass[5])
    test_recall = c(test_recall,model2_2$pred_table[[1]]$byClass[6])
    test_F1 = c(test_F1,model2_2$pred_table[[1]]$byClass[7])
  }

  a2 = c(mean(test_Accuracy,na.rm = T),mean(test_Sensitivity,na.rm = T),
         mean(test_Specificity,na.rm = T),mean(test_Precision,na.rm = T),
         mean(test_recall,na.rm = T),mean(test_F1,na.rm = T))
  a2 = round(a2,3)
  test_Accuracy = c()
  test_Sensitivity = c()
  test_Specificity = c()
  test_Precision = c()
  test_recall = c()
  test_F1 = c()
  for (i in 1:100) {
    set.seed(i*10)
    rs1 = sample(1:MSI_size,groupn,replace=F)
    set.seed(i*10)
    rs2 = sample((MSI_size+1):total_size,groupn*MIS_multi,replace=F)
    rs = c(rs1,rs2)
    names(ppo2$x)[names(ppo2$x)=="3-phosphoglycerate"]="X3.phosphoglycerate"
    names(ppo2$x)[names(ppo2$x)=="Glutathione reduced"]="Glutathione.reduced"
    names(ppo2$x)[names(ppo2$x)=="6-phosphogluconate"]="X6.phosphogluconate"
    names(ppo2$x)[names(ppo2$x)=="C18:0 lysophosphatidylethanolamine (LPE)"]="C18.0.LPE"
    names(ppo2$x)[names(ppo2$x)=="C36:1 phosphatidylcholine (PC)"]="C36.1.PC"
    names(ppo2$x)[names(ppo2$x)=="C14:0 cholesterol ester (CE)"]="C14.0.CE"
    train_tensor = rescale_standard(ppo2$x)[-rs,]
    valid_tensor = rescale_standard(ppo2$x)[rs,]
    train_label = label2$MSIStatus[-rs]
    valid_label = label2$MSIStatus[rs]
    train_cov_matrix = as.matrix(ppo2$z[-rs,])
    valid_cov_matrix = as.matrix(ppo2$z[rs,])
    train_tensor = cbind(train_tensor,train_cov_matrix)
    valid_tensor = cbind(valid_tensor,valid_cov_matrix)
    train_tensor$label = factor(train_label)
    valid_tensor$label = factor(valid_label)
    rf.model = suppressWarnings(randomForest(label~.,data=train_tensor))
    
    cm = confusionMatrix(predict(rf.model,valid_tensor),valid_tensor$label)
    
    test_Sensitivity = c(test_Sensitivity,cm$byClass[1])
    test_Specificity = c(test_Specificity,cm$byClass[2])
    test_Accuracy = c(test_Accuracy,cm$overall[1])
    test_Precision = c(test_Precision,cm$byClass[5])
    test_recall = c(test_recall,cm$byClass[6])
    test_F1 = c(test_F1,cm$byClass[7])
  }
  #===
  a5 = c(mean(test_Accuracy,na.rm = T),mean(test_Sensitivity,na.rm = T),
         mean(test_Specificity,na.rm = T),mean(test_Precision,na.rm = T),
         mean(test_recall,na.rm = T),mean(test_F1,na.rm = T))
  a5 = round(a5,3)
  b5=a5[-5]
  b2=a2[-5]
  table2=rbind(b5,b2)
  rownames(table2)=c("Random Forest","Catch")
  colnames(table2)=c("Accuracy","Sensitivity","Specificity","Precision","F1")
  return(table2)
}
associated=function(alpha_matrix,pathway2){
  pathname=c("amino acids and derivatives","amino acids and derivatives","lipids","carbohydrates","lipids","lipids","carbohydrates","amino acids and derivatives")
  sigfeature=c()
  for(i in 1:8){
    headf=rownames(data.frame(head(sort(alpha_matrix[i,colnames(alpha_matrix[i,,drop=FALSE]) %in% row.names(pathway2)[pathway2$pathway==pathname[i]]]),3)))
    tailf=rownames(data.frame(tail(sort(alpha_matrix[i,colnames(alpha_matrix[i,,drop=FALSE]) %in% row.names(pathway2)[pathway2$pathway==pathname[i]]]),3)))
    sigfeature[[i]]=c(tailf,headf)
  }
  nf=rownames(alpha_matrix)
  testdf=as.data.frame(pathname,nf)
  featurestring=c()
  for(i in 1:8){
    featurestring[i]=paste(sigfeature[[i]][3],sigfeature[[i]][2],sigfeature[[i]][1],sigfeature[[i]][4],sigfeature[[i]][5],sigfeature[[i]][6],sep="/")
  }
  testdf=cbind(testdf,featurestring)
}