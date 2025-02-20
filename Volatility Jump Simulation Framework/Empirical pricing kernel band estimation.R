###########################################################

K = function (u)
{
  kk=length(u);
  ret=vector(, length=kk);
  for (i in c(1:kk))
  {
    ret[i]=dnorm(u[i]);         # Gaussian
    #       if (abs(u[i])>1) {ret[i]=0;} else {ret[i]=(3*(1-u[i]^2)/4) };       # Epanechnikov
  }
  return(ret);
}

K.h = function (u, h)
{
  kk=length(u);
  ret=vector(, length=kk);
  for (i in c(1:kk)){ret[i]=K(u[i]/h)/h;}
  return(ret);
}

###########################################################

rookley = function(money, sigma, sigma1, sigma2, r, tau) 
{
  rm=length(money);
  sqrttau=sqrt(tau);
  exprt=exp(r*tau);
  rtau=r*tau;
  
  d1=(log(money)+tau*(r+0.5*sigma^2))/(sigma*sqrttau);
  d2=d1-sigma*sqrttau;
  
  d11 = 1/(money*sigma*sqrttau)+(-(log(money)+tau*r)/(sqrttau*sigma^2)+0.5*sqrttau)*sigma1;
  
  d21 = 1/(money*sigma*sqrttau)+(-(log(money)+tau*r)/(sqrttau*sigma^2)-0.5*sqrttau)*sigma1;
  
  d12 = -1/(sqrttau*money^2*sigma)-sigma1/(sqrttau*money*sigma^2)+ sigma2*(0.5*sqrttau-(log(money)+r*tau)/(sqrttau*sigma^2))+ sigma1*(2*sigma1*(log(money)+r*tau)/(sqrttau*sigma^3)-1/(sqrttau*money*sigma^2));
  d22 = -1/(sqrttau*(money^2)*sigma)-sigma1/(sqrttau*money*sigma^2) + sigma2*(-0.5*sqrttau-(log(money)+r*tau)/(sqrttau*sigma^2))+ sigma1*(2*sigma1*(log(money)+r*tau)/(sqrttau*sigma^3)-1/(sqrttau*money*sigma^2));
  
  f = pnorm(d1)-pnorm(d2)/(money*exprt);
  f1=dnorm(d1)*d11-dnorm(d2)*d21/(exprt*money)+pnorm(d2)/(exprt*money^2);
  f2=dnorm(d1)*d12-d1*dnorm(d1)*d11^2-dnorm(d2)*d22/(exprt*money)+dnorm(d2)*d21/(exprt*money^2) + d2*dnorm(d2)*d21^2/(exprt*money)-2*pnorm(d2)/(exprt*money^3)+dnorm(d2)*d21/(exprt*money^2);
  
  #cat("d1 ", d1[ngrid]," d2 ",d2[ngrid]," d11 ", d11[ngrid]," d12 ", d12[ngrid]," d21 ",d21[ngrid]," d22 ",d22[ngrid])
  return(cbind(f, f1, f2));
}

#########################       kernel functions
K = function (u)
{
  kk=length(u);
  ret=vector(, length=kk);
  for (i in c(1:kk))
  {
    ret[i]=dnorm(u[i]);         # Gaussian
    #       if (abs(u[i])>1) {ret[i]=0;} else {ret[i]=(3*(1-u[i]^2)/4) };       # Epanechnikov
  }
  return(ret);
}

K.h = function (u, h)
{
  kk=length(u);
  ret=vector(, length=kk);
  for (i in c(1:kk)){ret[i]=K(u[i]/h)/h;}
  return(ret);
}


K.prime = function (u)
{
  kk=length(u);
  ret=vector(, length=kk);
  for (i in c(1:kk))
  {
    ret[i]=-u[i]*dnorm(u[i]);               # Gaussian
    #       if(abs(u[i])>1) {ret[i]=0} else {ret[i]=(-3*u[i]/2)};       #   Epanechnikov
  }
  return(ret);
}

Kjp = function (u, j, N)
{
  kk=length(u);
  ret=vector(, length=kk);
  
  for (i in c(1:kk))
  {
    M = N;
    M[,j] = u[i]^{c(0:(dim(N)[1]-1))};
    if(abs(u[i])>1) {ret[i]=0;} else {ret[i]= K(u[i])*det(M)/det(N); };
  }
  return(ret);
}



######################################## N, Q, C matrices

#p=2;h=0.65; alpha=0.05;
func.L = function(p, h, alpha)
{
  N = matrix(0, nrow=p+1, ncol=p+1);
  Q = matrix(0, nrow=p+1, ncol=p+1);
  C = vector(, length=p+1);
  L = vector(, length=p+1);
  
  for(i in c(0:p))
  {
    for(j in c(0:p))
    {
      if ((i>1) && (j>1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value- 0.5 *(i*(i-1)+j*(j-1))*  integrate(function (x) x^(i+j-2)*(K(x))^2, -Inf,Inf)$value;}
      if ((i<=1) && (j<=1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value;}
      if ((i<=1) && (j>1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value- 0.5 *j*(j-1)*  integrate(function (x) x^(i+j-2)*(K(x))^2, -Inf,Inf)$value;}
      if ((i>1) && (j<=1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value- 0.5 *i*(i-1)*  integrate(function (x) x^(i+j-2)*(K(x))^2, -Inf,Inf)$value;}
      N[i+1,j+1]= integrate(function (xx) xx^(i+j)*K(xx), -1,1)$value;
    }
  }
  
  N.i = solve(N);
  
  for(j in c(0:p))
  { 
    C[j+1]= (N.i%*%Q%*%N.i)[j+1,j+1] /  integrate(function (xx) (Kjp(xx, j, N))^2, -Inf,Inf)$value;    
    L[j+1] = factorial(j)* (n*h^(2*j+1))^(-0.5) * ( (-2*log(h))^(1/2) + (-2*log(h))^(-1/2)*(-log(-0.5*log(1-alpha)) + log(sqrt(C[j+1])/(2*pi))));
  }
  print(C);
  return(L);
}

######################################### asymptotic bands


bands = function(p, h, hat.theta, x, y, vx) {
  L = func.L(p, h, 0.05)
  sigma2 = 1
  up.band = matrix(0,nrow=length(vx), ncol=p+1)
  lo.band = matrix(0,nrow=length(vx), ncol=p+1)
  V.all = matrix(0,nrow=length(vx), ncol=p+1)
  
  for(i in c(1:length(vx))) {
    xx = vx[i]
    B.V = matrix(0, nrow=p+1, ncol=p+1)
    K.V = matrix(0, nrow=p+1, ncol=p+1)
    
    for(j in c(1:length(x))) {
      residuals = y[j]-sum(hat.theta[i,]*((x[j]-xx)^c(0:p)))
      if(i == 1 && j == 1) {  
        print("Residuals:")
        print(residuals)
      }
      
      B.V = B.V + K.h(x[j]-xx,h) * (1/sigma2) * (solve(diag(h^c(0:p))) %*% t(t((x[j]-xx)^c(0:p))) %*% t((x[j]-xx)^c(0:p)) %*% solve(diag(h^c(0:p))))
      K.V = K.V + h * (K.h(x[j]-xx,h))^2 * ((y[j]-sum(hat.theta[i,]*((x[j]-xx)^c(0:p))))/sigma2)^2 * (solve(diag(h^c(0:p))) %*% t(t((x[j]-xx)^c(0:p))) %*% t((x[j]-xx)^c(0:p)) %*% solve(diag(h^c(0:p))))
    }
    
    if(i == 1) {  
      print("K.V matrix:")
      print(K.V)
    }
    
    B.V = B.V/n
    K.V = K.V/n
    B.i = solve(B.V)
    V = B.i%*%K.V%*%B.i
    V.all[i,] = diag(V)
    
    for(ip in c(1:(p+1))) {
      lo.band[i,ip] = -L[ip] * sqrt(V[ip, ip])
      up.band[i,ip] = L[ip] * sqrt(V[ip, ip])
    }
  }
  
  print("Width of confidence intervals:")
  print(up.band - lo.band)
  
  mylist = list(L=L, V=V.all, lo.band=lo.band, up.band=up.band)
  return(mylist)
}

bootstrap_epk_bands = function(data, spd_values, B, h, alpha) {
  n = length(data)
  ngrid = length(spd_values)
  
  boot_densities = matrix(0, nrow=B, ncol=ngrid)
  boot_epk = matrix(0, nrow=B, ncol=ngrid)
  
  for(b in 1:B) {
    boot_sample = sample(data, n, replace=TRUE)
    kde_b = density(boot_sample, bw=h, from=min(mon.grid.scaled), 
                    to=max(mon.grid.scaled), n=ngrid)
    
    kde_values = pmax(kde_b$y, 1e-10) 
    boot_densities[b,] = kde_values
    boot_epk[b,] = spd_values/kde_values
  }
  
  mean_epk = colMeans(boot_epk, na.rm=TRUE)
  std_epk = apply(boot_epk, 2, sd, na.rm=TRUE)
  
  scaled_devs = abs((boot_epk - matrix(mean_epk, nrow=B, ncol=ngrid, byrow=TRUE)) / 
                      matrix(std_epk, nrow=B, ncol=ngrid, byrow=TRUE))
  max_devs = apply(scaled_devs, 1, max, na.rm=TRUE)
  max_devs = max_devs[is.finite(max_devs)] 
  
  L = quantile(max_devs, 1-alpha, na.rm=TRUE)
  
  upper = mean_epk + L * std_epk
  lower = mean_epk - L * std_epk
  
  return(list(lower=lower, upper=upper, L=L))
}


setwd("") #SET PATH
output_path <- ""  #SET PATH
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)


ngrid = 200;
bandwidth = 0.12;
tau = 0.04722;

read.dates = FALSE;

#all = read.table("tick2006.txt")
#day1 = all[all[,1]=="28-02-2006",];
#day2 = all[all[,1]=="04-04-2006",];
#write.table(day1, file="28-02-2006.txt", row.names=F, col.names=F, quote=F)
#write.table(day2, file="04-04-2006.txt", row.names=F, col.names=F, quote=F)

XX  = read.table("tick2006.txt");
#XX  = read.table("C20010117.dat");
#XX = XX[c(2,1,3,5,6,8,7,4)];
#XX = cbind(rep("17-01-2001", length=dim(XX)[1]), XX)

cat("reading finished ....\n");

if (read.dates==T)
{
  date = as.vector(XX[,1]);
  date.dif = date[1];
  date.part=date;
  flag=1;
  
  while(flag==1)
  {
    date.part = date.part[date.part!=date.dif[length(date.dif)]];
    cat("current date ... ",    date.dif[length(date.dif)], "\n");
    if (length(date.part>0)) {date.dif=c(date.dif, date.part[1]);} else {flag=0;}    
  }
  write.table(date.dif, file="tick2006_dates.txt", row.names=F, col.names=F, quote=F);
} else { date.dif= read.table("tick2006_dates.txt"); date.dif=date.dif[,1];}

#############################################

#for( iday in c(1:length(date.dif)))
{
  iday= 1182;      #2.1.2020 1176 // 10.01.2020 1182
  
  day1 = XX[XX[,1]==date.dif[iday],]; 
  day1 = day1[day1[,2]>0,];
  tau = day1[1,4];
  #name = "28-02-2006";
  #day1=read.table(paste(name,".txt", sep=""));
  
  day1.mat = day1[day1[,4]== tau,];
  day1.call = day1.mat[day1.mat[,3]==1,];
  day1.put = day1.mat[day1.mat[,3]==0,];
  
  # compute the moneyness
  day1.call[,9]=day1.call[,7]/day1.call[,5];
  day1.put[,9]=day1.put[,7]/day1.put[,5];
  # only out and at the money options #MONEYNESS FILTER TEST
  #day1.call = day1.call[day1.call[,9]>=1,];
  #day1.put = day1.put[day1.put[,9]<=1,];
  
  # put to calls
  put2call = day1.put;
  put2call[,6] = put2call[,6]+ put2call[,7] -  put2call[,5]*exp(- mean(put2call[,8])* tau);
  
  put2call = put2call[order(put2call[,5]),];
  day1.call = day1.call[order(day1.call[,5]),];
  data = rbind(put2call,day1.call);
  
  data = data[(data[,2]>0.05),];
  
  # no - arbitrage condition
  #data = data[ ((data[,6]<data[,7]) & (data[,6]>data[,7]-data[,5]*exp(-data[,8]*data[,4])) ),]
  
  #write.table(data, file=paste(name,"_cleaned.txt", sep=""), row.names=F, col.names=F, quote=F)
  
  n = dim(data)[1];
  
  price.median = median(data[,7]);
  
  ## regression for volatility 
  
  volas = data[,c(2,9)];
  #cheat
  #volas[,1] <- volas[,1] + rnorm(nrow(volas), mean = 0, sd = 0.02)
  
  nsample = dim(volas)[1]; 
  
  mon.min = min(volas[,2]);mon.max = max(volas[,2]);
  volas[,2] = (volas[,2]-mon.min)/(mon.max-mon.min);
  mon.grid = seq(1/ngrid, 1, length=ngrid);
  mon.grid.scaled =  (mon.min+mon.grid*(mon.max-mon.min));
  
  for(i in c(1:ngrid)) {
    X = cbind(rep(1, length=nsample), 
              volas[,2]-mon.grid[i], 
              (volas[,2]-mon.grid[i])^2, 
              (volas[,2]-mon.grid[i])^3)
    W = diag(K.h(volas[,2]-mon.grid[i], bandwidth))
    beta = t(solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% volas[,1])
    
    if (i==1) {
      print("First fitted values:")
      print(volas[,1])  # Actual values
      print("First beta coefficients:")
      print(beta)
      sigmas=beta[c(1:3)]
    } else {
      sigmas=rbind(sigmas, beta[c(1:3)])
    }
  }
  
  sigmas[,3] = 2*sigmas[,3];
  
  #### bands for the derivatives
  
  bands.init = bands(2, bandwidth, sigmas, volas[,2], volas[,1], mon.grid);
  
  #### applying rookley
  
  fder = rookley(mon.grid.scaled, sigmas[,1], sigmas[,2]/(mon.max-mon.min), sigmas[,3]/(mon.max-mon.min)^2, mean(data[,8]), tau);
  
  #### final computations
  
  strike.grid  =  price.median / mon.grid.scaled;
  d2fdX2 = (mon.grid.scaled^2*fder[,3]+2*mon.grid.scaled*fder[,2])/strike.grid^2;
  
  #### bands for SPD
  
  d1=(log(mon.grid.scaled)+tau*(mean(data[,8])+0.5*(sigmas[,1])^2)) / (sqrt(sigmas[,1])*sqrt(tau));
  d2=d1-sqrt(sigmas[,1])*sqrt(tau);
  dgds= exp(-mean(data[,8])*tau) * ( mon.grid.scaled^2 * dnorm(d1) / strike.grid^2 - exp(-mean(data[,8])*tau)*dnorm(d2)/mon.grid.scaled);
  
  band.limit = bands.init$L[3] * sqrt(bands.init$V[,3]) * dgds ; #* exp(-mean(data[,8])*tau); 
  
  SPD = new.env();
  SPD$SPD = price.median^2 * exp(mean(data[,8])*tau) * d2fdX2 / mon.grid.scaled;
  SPD$lo = SPD$SPD - abs(price.median^2 * exp(mean(data[,8])*tau) * d2fdX2 * band.limit / mon.grid.scaled);
  SPD$up = SPD$SPD + abs(price.median^2 * exp(mean(data[,8])*tau) * d2fdX2 * band.limit / mon.grid.scaled);
  SPD=as.list(SPD);
  
  #Date filter test für HD berechnung
  #dax=read.table("dax_index.dat"); dax=dax[,2];
  #dax=dax[c((length(dax)-500):length(dax))];
  dax_full = read.table("dax_index.dat")
  current_date = as.Date(date.dif[iday], format="%d-%m-%Y")
  dax_dates = as.Date(as.character(dax_full[,1]), format="%d/%m/%Y")
  
  # filter data so only relevant parts are being considered
  historical_dax = dax_full[dax_dates < current_date, 2]
  
  n_hist = min(1000, length(historical_dax))
  dax = tail(historical_dax, n_hist)
  
  
  
  
  print("tau:")
  print(tau)
  print("tau * 31/0.0833:")
  print(tau * 31/0.0833)
  print("ceiling(tau * 31/0.0833):")
  print(ceiling(tau * 31/0.0833))
  print("length(dax):")
  print(length(dax))
  print("dax head:")
  print(head(dax))
  print("dax tail:")
  print(tail(dax))
  dax.scaled = 1+log(dax[c(ceiling(tau*31/0.0833):length(dax))]/dax[c(1:(length(dax)+1-ceiling(tau*31/0.0833)))]);
  
  print("Anzahl der Samples für KDE:")
  print(length(dax.scaled))
  #kern=density(dax.scaled, bw="nrd0", from=mon.grid.scaled[1], to=mon.grid.scaled[ngrid], n=ngrid);
  kern=density(dax.scaled, bw= "nrd0", from=mon.grid.scaled[1], to=mon.grid.scaled[ngrid], n=ngrid);
  
  ########### Black-Scholes
  
  #sigma.BS = mean((volas[volas[,1]<quantile(volas[,1], probs=0.75),])[,1]);
  sigma.BS = mean(volas[,1])
  
  print("sigma.BS:")
  print(sigma.BS)
  
  print("tau:")
  print(tau)
  
  print("mean(data[,8]):")
  print(mean(data[,8]))
  
  print("mon.grid.scaled head:")
  print(head(mon.grid.scaled))
  q.BS = exp(-(log(mon.grid.scaled) - (mean(data[,8])-sigma.BS/2)*tau)^2/(2*tau*sigma.BS))/(mon.grid.scaled*sqrt(2*3.1415926*sigma.BS*tau));
  p.BS = dnorm(log(mon.grid.scaled), mean=mean(log(dax.scaled)), sd = sd(log(dax.scaled)))/mon.grid.scaled;
  
  #epk_bands = bootstrap_epk_bands(dax.scaled, SPD$SPD, B=1000, h=kern$bw, alpha=0.05)
  epk_bands = bootstrap_epk_bands(dax.scaled, SPD$SPD, B=3000, h= "nrd0", alpha=0.05)
  
  
  ########### plotting

  plot(mon.grid.scaled, SPD$SPD/kern$y, type="l", lwd=2, col=2,
       xlab="moneyness", ylab="EPK",
       main=paste(date.dif[iday], " tau = ", tau),
       xlim=c(.85, 1.13), ylim=c(0,10))
  lines(mon.grid.scaled, epk_bands$lower, col="blue3", lwd=1)
  lines(mon.grid.scaled, epk_bands$upper, col="blue3", lwd=1)
  #lines(mon.grid.scaled, q.BS/p.BS, col="green3", lwd=2)
  
  mod = lm(log(SPD$SPD)~ log(mon.grid.scaled));
  util.coef = mod$coef;
  
  if (exists("util.param")==0) {
    util.param = c(exp(util.coef[1]), -util.coef[2], cor(log(SPD$SPD), log(mon.grid.scaled))^2);
  } else {
    util.coef = cbind(util.coef,c(exp(util.coef[1]), -util.coef[2], cor(log(SPD$SPD), log(mon.grid.scaled))^2));
  }
  
} # end of the days#





epk_results <- data.frame(
  moneyness = mon.grid.scaled,
  epk = SPD$SPD/kern$y,
  lower_ci = epk_bands$lower,
  upper_ci = epk_bands$upper
)
#save as .csv
write.csv(epk_results, 
          file = paste0(output_path, "/epk_results_", 
                        format(as.Date(date.dif[iday], "%d-%m-%Y"), "%Y%m%d"),
                        "_tau", round(tau*365, 0), ".csv"),
          row.names = FALSE)