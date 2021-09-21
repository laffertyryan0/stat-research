library(arrangements)
library(invgamma)
library(rjson)
C = function(j, r, h, Upaj,Upajj){
  s = 2*rep(h[j],length(Upajj))
  Upajm1 = solve(Upaj)
  term = function(nu){
    out = (-1)^(sum(abs(nu)))
    for(k in 1:length(nu)){
      out = out*choose(s[k],nu[k])
    }
    fp = t(s/2 - nu)%*%Upajm1%*%(s/2-nu)
    sp = t(s/2-nu)%*%Upajm1%*%Upajj
    denom = factorial(r)*factorial((sum(abs(s))-2*r))
    out  = out * (fp^r * sp^(sum(abs(s))-2*r) / denom)
    return(out)
  }
  out2 = 0
  perm = permutations(x = 0:s[1],k=s[1],replace=TRUE)
  for(i in 1:length(perm[,1])){
    out2 = out2 + term(perm[i,])
  }
  return(out2[1,1])
}

sampDmarginal = function(j,Ujj,Upaj,Upajj,h,alpha){
  s = 2*rep(h[j],length(Upajj))
  maxr = floor(sum(abs(s))/2)
  probs = c()
  shapes = c()
  Upajm1 = solve(Upaj)
  mj = length(s)
  for(k in 1:(maxr+1)){
    r = k-1
    probs[k] = C(j,r,h,Upaj,Upajj)
    shapes[k] = -1 - (r + mj/2 - alpha[j]/2)
  }
  betaj = c((1/2)*(Ujj - t(Upajj)%*%Upajm1%*%Upajj))
  probs = probs/sum(probs)
  rchoice = sample(0:maxr,size = 1,prob = probs)
  #print(shapes[rchoice])
  samp = rinvgamma(1,shapes[rchoice],rate = betaj)
  return(samp)
}

#sample from f(x) = (x+a)^k * standardnormalpdf(x)
args <- fromJSON(file = ".json/configTable.json")
aMin = args$aMin
aMax = args$aMax
aStep = args$aStep
xMin = args$xMin
xMax = args$xMax
xStep = args$xStep
sampF = function(k,a){
  u = runif(1,0,1)
  path = paste("./table/tablekequals",k,".csv",sep="")
  line = read.csv(file = path,
                  header = TRUE, 
                  nrows = 1,
                  check.names = FALSE,
                  row.names = 1,
                  skip = as.integer((a-aMin)/aStep))
  leftendpoint = findInterval(u,line)
  rightendpoint = leftendpoint + 1
  left_xval = xMin + xStep*(leftendpoint - 1)
  right_xval = xMin + xStep*leftendpoint
  left_Fval = if(leftendpoint != 0) line[[leftendpoint]] else 0
  right_Fval = if(rightendpoint != (length(line)+1)) line[[rightendpoint]] else 1
  slope = (right_Fval-left_Fval)/(right_xval-left_xval)
  approx_inv = if(slope!=0) left_xval + (u - left_Fval)/slope else left_xval
  return(approx_inv)
}
pp = c();
system.time(for(i in 1:1000){pp = append(pp,sampF(4,.3))})

Upajj = c(0, 1)
Upaj = matrix(c(2,0,0,2),nrow = 2)
h = matrix(c(1,0,0,1),nrow = 2)
alpha = 10*matrix(c(1,1,1,1),nrow = 2)


C(1,1,h,Upaj,Upajj)
samps = c()
for(i in 1:1000){
  samps[i] = sampDmarginal(1,10,Upaj,Upajj,h,alpha)
}
hist(samps,breaks = 100)
mean(samps)
mode(samps)

# - a - 1 = pow -> -a = pow + 1 -> a = - pow - 1



