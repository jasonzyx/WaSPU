impute = function(X){
    
    for(k in 1:ncol(X)) {
        a = X[,k]
        a[is.na(a)]<- mean(a,na.rm=T)
        X[,k] <- a
    }
    
    return(X)
}
