
## @knitr setup 

data=read.table("Givens_Data/easysmooth.dat",header=TRUE)
x=data[,1]
y=data[,2]

## @knitr plotdata

plot(data,ylab="Response",xlab="Predictor",pch=19)
#lines(x,x^3*sin((x+3.4)/2),col="green")

## @knitr rmimpl

RunningMean=function(k,data){
  #data is a n by 2 matrix with x_i in the 1st column
  #and y_i in the 2nd column. x_i's are assumed sorted
  n=length(data[,1])
  Y=data[,2]
  #build the span matrix row by row.
  S=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    #get indices of the symmetric nearest neighbors
    left=max(1,(i-(k-1)/2));right=min(n,(i+(k-1)/2))
    len=right-left+1
    S[i,left:right]=rep(1/len,len)
  }
  S%*%Y #Return the running average
}

## @knitr prints
n=8
k=5
#build the span matrix row by row.
S=matrix(0,ncol=n,nrow=n)
for(i in 1:n){
  #get indices of the symmetric nearest neighbors
  left=max(1,(i-(k-1)/2));right=min(n,(i+(k-1)/2))
  len=right-left+1
  S[i,left:right]=rep(1/len,len)
}
print(S)


## @knitr rm13plot
k13=RunningMean(13,data)
plot(data,ylab="Response",xlab="Predictor",pch=19, main="Running Mean with k=13")
x=data[,1]
lines(x,k13,col="blue")
#plot the true curve
lines(x,x^3*sin((x+3.4)/2),col="green")
legend(2,-5,legend=c("k=13","True curve"),lty=c(1,1),col=c("blue","green"))

## @knitr comparespan

k3=RunningMean(3,data)
k50=RunningMean(49,data)
plot(data,ylab="Response",xlab="Predictor",pch=19, main="Running Mean Comparison")
x=data[,1]
lines(x,k3,col="blue")
lines(x,k50,col="red")
#plot the true curve
lines(x,x^3*sin((x+3.4)/2),col="green")
legend(2,-5,legend=c("k=3","k=49","True curve"),lty=c(1,1,1),col=c("blue","red","green"))

## @knitr bvplot
plot(data,ylab="Response",xlab="Predictor",pch=19, xlim=c(-3,-1.5), ylim=c(-7.2,-4))
k3=RunningMean(3,data)
k50=RunningMean(49,data)
lines(x,k3,col="blue")
lines(x,k50,col="red")
#plot the true curve
lines(x,x^3*sin((x+3.4)/2),col="green")
legend(-1.8,-6.5,legend=c("k=3","k=49","True curve"),lty=c(1,1,1),col=c("blue","red","green"))

## @knitr docvcode
DoCrossValidate=function(k,data){
  n=length(data[,1])
  S=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    #get indices of the symmetric nearest neighbors
    left=max(1,(i-(k-1)/2));right=min(n,(i+(k-1)/2))
    len=right-left+1
    S[i,left:right]=rep(1/len,len)
  }
  Y=data[,2]
  mean(((Y-S%*%Y)/(1-diag(S)))^2)
}

## @knitr cvplot
DoCrossValidate=function(k,data){
  n=length(data[,1])
  S=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    #get indices of the symmetric nearest neighbors
    left=max(1,(i-(k-1)/2));right=min(n,(i+(k-1)/2))
    len=right-left+1
    S[i,left:right]=rep(1/len,len)
  }
  Y=data[,2]
  mean(((Y-S%*%Y)/(1-diag(S)))^2)
}
ks=seq(3,49,by=2)
cvs=numeric(length(ks))
for(i in 1:length(ks)){
  k=ks[i]
  cvs[i]=DoCrossValidate(k,data)
}
plot(ks,cvs,type='b',pch=19,
     xlab="k",ylab="CVRSS",main="CVRSS vs. k")


## @knitr rlimpl
RunLine=function (k,data){
  x=data[,1];y=data[,2];n=length(x)
  vals=matrix(nrow=n,ncol=1)
  S=matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    #get indices of the symmetric nearest neighbors
    left=max(1,(i-(k-1)/2));right=min(n,(i+(k-1)/2))
    len=right-left+1
    Xi=matrix(1,ncol=2,nrow=len)
    Xi[,2]=x[left:right]
    Yi=y[left:right];
    Hi=Xi%*%solve(t(Xi)%*%Xi)%*%t(Xi)
    S[i,left:right]=Hi[i-left+1,] #Select the row
    #corresponding x_i
  }
  vals=S%*%y
}

## @knitr rlplot
y=RunLine(23,data)
plot(data,ylab="Response",xlab="Predictor",pch=19)
lines(data[,1],y,col="blue")
lines(x,x^3*sin((x+3.4)/2),col="green")
legend(1,-5,legend=c("Running Line (k=23)","True curve"),lty=c(1,1),col=c("blue","green"))

## @knitr rlcvplot
RunLineCV=function (k,data){
  x=data[,1];y=data[,2];n=length(x)
  vals=matrix(nrow=n,ncol=1)
  S=matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    #get indices of the symmetric nearest neighbors
    left=max(1,(i-(k-1)/2));right=min(n,(i+(k-1)/2))
    len=right-left+1
    Xi=matrix(1,ncol=2,nrow=len)
    Xi[,2]=x[left:right]
    Yi=y[left:right];
    Hi=Xi%*%solve(t(Xi)%*%Xi)%*%t(Xi)
    S[i,left:right]=Hi[i-left+1,]
  }
  vals=mean(((y-S%*%y)/(1-diag(S)))^2)
}
ks=seq(3,45,by=2)
cv=lapply(ks,function(x){RunLineCV(x,data)})
plot(ks,cv,type="b", pch=19)

## @knitr csplot
getCubicSmoothMat = function(lambda,data){
  x=data[,1];y=data[,2];n=length(x)
  h=diff(x)
  W=diag(n)
  Tb=matrix(0,ncol=n-2,nrow=n-2)
  Q=matrix(0,ncol=n,nrow=n-2)
  for(i in 1:(n-2)){
    Tb[i,i]=(h[i]+h[i+1])/3
    if(i>1){Tb[i-1,i]=h[i]/6;Tb[i,i-1]=h[i]/6}
    Q[i,i]=1/h[i]
    if(i<n-2) Q[i,i+1]=-(1/h[i]+1/h[i+1])
    if(i<n-3)Q[i,i+2]=1/h[i+1]
  }
  
  K=t(Q)%*%solve(Tb)%*%Q
  S=solve(diag(n)+lambda*K)
}

doCubicSpline=function(lambda,data){
  S=getCubicSmoothMat(lambda,data);
  S%*%y
}
yhat=doCubicSpline(.066,data)
plot(data,ylab="Response",xlab="Predictor",pch=19)
lines(x,yhat,col="blue")
lines(x,x^3*sin((x+3.4)/2),col="green")
legend(2,-5,legend=c(expression(lambda==0.066),"True curve"),lty=c(1,1,1),col=c("blue","green"))

## @knitr plotcscv
cubicSplineCV=function (lambdas,data){
  cvs=numeric(length(lambdas))
  i=1
  Y=data[,2]
  n=length(Y)
  for(l in lambdas){
    S=getCubicSmoothMat(l,data)
    S=S[,1:(n-2)]
    Y=Y[1:(n-2)]
    cvs[i]=mean(((Y-S%*%Y)/(1-diag(S)))^2)
    i=i+1
  }
  cvs
}
lambdas=seq(.005,1,by=.005)
cvs=cubicSplineCV(lambdas,data)
plot(lambdas,cvs)