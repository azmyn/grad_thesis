df=read.table("/Users/azumi/Dropbox/Labs/8semester/grad_thesis/grad_source/lasso/crime.txt", sep="")
X = scale(df[,3:7])
y = scale(df[,1])

S.lambda=function(x,lambda)max(abs(x)-lambda,0)*sign(x)
inner.prod=function(x,y)sum(x*y)

linear.reg=function(x,y,lambda){
  q=ncol(x)
  z=diag(t(x)%*%x)
  beta.hat=array(0, dim=q)
  beta.hat.old=array(10, dim=q)
  beta.hat.0=0
  eps=1
  while(eps>0.01){
    for(j in 1:q){
      r=y-beta.hat.0-x[,-j]%*%beta.hat[-j]
      beta.hat[j]=S.lambda(inner.prod(r,x[,j]),lambda)/z[j]
    }
    beta.hat.0=sum(y)/n-sum(x%*%beta.hat)/n
    eps=max((beta.hat-beta.hat.old)/beta.hat.old)
    beta.hat.old=beta.hat
  }
  return(list(beta=beta.hat,beta.0=beta.hat.0))
}


range=c(0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0, 20, 25, 30)
p=ncol(X)
n=nrow(X)
q=length(range)
tab=array(dim=c(p,q))
for(i in 1:q)tab[,i]=linear.reg(X,y,range[i])$beta

cols <- c("funding", "hs", "not-hs", "college", "college4")
plot(0,xlim=c(min(range),max(range)), xlab="lambda", ylab="coefficient", ylim=c(-0.6,0.6), type="n")
for(j in 1:p)lines(range,tab[j,],col=j)
legend("topright", legend=cols, col=1:p, lty=1)
