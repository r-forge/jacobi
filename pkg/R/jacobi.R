#
#   jacobi package
#   Copyright (C) 2008  Jan de Leeuw <deleeuw@stat.ucla.edu>
#   UCLA Department of Statistics, Box 951554, Los Angeles, CA 90095-1554
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################
#
# version 0.0.1, 2008-12-05   Eigen and SVD
# version 0.1.0, 2008-12-06   Simultaneous Diagonalization
# version 0.2.0, 2008-12-08   Simultaneous SVD
# version 0.2.1, 2008-12-08   Various Bugfixes
# version 0.2.2, 2008-12-08   Small efficiency gains
# version 0.3.0, 2008-12-10   PREHOM in jMCA added
# version 1.0.0, 2008-12-10   Tucker3 Methods added
#

jEigen<-function(a,eps1=1e-6,eps2=1e-10,itmax=100,vectors=TRUE,verbose=FALSE) {
n<-nrow(a); k<-diag(n); itel<-1; mx<-0; saa<-sum(a^2)
repeat {
    for (i in 1:(n-1)) for (j in (i+1):n) {
        aij<-a[i,j]; bij<-abs(aij);
        if (bij < eps1) next()
        mx<-max(mx,bij)
        am<-(a[i,i]-a[j,j])/2
        u<-c(aij,-am); u<-u/sqrt(sum(u^2))
        c<-sqrt((1+u[2])/2); s<-sign(u[1])*sqrt((1-u[2])/2) 
        ss<-s^2; cc<-c^2; sc<-s*c
        ai<-a[i,]; aj<-a[j,]
        aii<-a[i,i]; ajj<-a[j,j]
        a[i,]<-a[,i]<-c*ai-s*aj
        a[j,]<-a[,j]<-s*ai+c*aj
        a[i,j]<-a[j,i]<-0
        a[i,i]<-aii*cc+ajj*ss-2*sc*aij
        a[j,j]<-ajj*cc+aii*ss+2*sc*aij 
        if (vectors) {
            ki<-k[,i]; kj<-k[,j]
            k[,i]<-c*ki-s*kj
            k[,j]<-s*ki+c*kj
            }
        }
    ff<-sqrt(saa-sum(diag(a)^2))
    if (verbose)
        cat("Iteration ",formatC(itel,digits=4),"maxel ",formatC(mx,width=10),"loss ",formatC(ff,width=10),"\n")
    if ((mx < eps1) || (ff < eps2) || (itel == itmax)) break()
    itel<-itel+1; mx<-0
    }
d<-diag(a); o<-order(d,decreasing=TRUE)
if (vectors) return(list(values=d[o],vectors=k[,o]))
    else return(values=d[o])
}

jSVD<-function(x,eps1=1e-6,eps2=1e-6,itmax=1000,vectors=TRUE,verbose=FALSE) {
n<-nrow(x); m<-ncol(x); itel<-1; mx<-0
kkk<-diag(n); lll<-diag(m);
sxx<-sum(x^2); sxm<-sqrt(sxx/(n*m))
repeat {
    for (i in 1:(n-1)) {
		if (i > m) next()		
		for (j in (i+1):n) {
	        xi<-x[i,]; xj<-x[j,]        
	        xij<-ifelse(j > m,0,x[i,j])
	        xjj<-ifelse(j > m,0,x[j,j])
			xii<-x[i,i]; xji<-x[j,i]
	        mx<-max(mx,abs(xij)/sxm,abs(xji)/sxm)
	        v<-matrix(0,2,2)
	        v[1,1]<-xij^2+xji^2
	        v[1,2]<-v[2,1]<-xii*xji-xjj*xij
	        v[2,2]<-xii^2+xjj^2
	        u<-eigen(v)$vectors[,1]
	        x[i,]<-u[2]*xi+u[1]*xj
	        x[j,]<-u[2]*xj-u[1]*xi
	        if (vectors) {
	            ki<-kkk[i,]; kj<-kkk[j,]
	            kkk[i,]<-u[2]*ki+u[1]*kj
	            kkk[j,]<-u[2]*kj-u[1]*ki
	            }
	        }
		}
    ff<-sqrt((sxx-sum(diag(x)^2))/sxx)
    if (verbose)
        cat(" Left iteration ",formatC(itel,digits=4),"maxel ",formatC(mx,width=10),"loss ",formatC(ff,width=10),"\n")
    for (k in 1:(m-1)) {
    	if (k > n) next()
		for (l in (k+1):m) {
	        xk<-x[,k]; xl<-x[,l]
	        xlk<-ifelse(l > n,0,x[l,k])
	        xll<-ifelse(l > n,0,x[l,l])
	        xkk<-x[k,k]; xkl<-x[k,l]
	        mx<-max(mx,abs(xkl)/sxm,abs(xlk)/sxm)
	        v<-matrix(0,2,2)
	        v[1,1]<-xkl^2+xlk^2
	        v[1,2]<-v[2,1]<-xll*xlk-xkk*xkl
	        v[2,2]<-xkk^2+xll^2
	        u<-eigen(v)$vectors[,1]
	        x[,k]<-u[2]*xk-u[1]*xl
	        x[,l]<-u[1]*xk+u[2]*xl
	        if (vectors) {
	            lk<-lll[,k]; ll<-lll[,l]
	            lll[,k]<-u[2]*lk-u[1]*ll
	            lll[,l]<-u[1]*lk+u[2]*ll
	            }
	        }
		}
    ff<-sqrt((sxx-sum(diag(x)^2))/sxx)
    if (verbose)
        cat("Right iteration ",formatC(itel,digits=4),"maxel ",formatC(mx,width=10),"loss ",formatC(ff,width=10),"\n")
    if ((mx < eps1) || (ff < eps2) || (itel == itmax)) break()
    itel<-itel+1; mx<-0
    }
return(list(d=diag(x),u=t(kkk),v=lll))
}

jSimDiag<-function(a,eps=1e-10,itmax=100,vectors=TRUE,verbose=FALSE) {
n<-dim(a)[1]; kk<-diag(n); m<-dim(a)[3]; itel<-1; saa<-sum(a^2)
fold<-saa-sum(apply(a,3,function(x) sum(diag(x^2))))
repeat {
    for (i in 1:(n-1)) for (j in (i+1):n) {
    ad<-(a[i,i,]-a[j,j,])/2
    av<-a[i,j,]
    v<-matrix(0,2,2)
    v[1,1]<-sum(ad^2)
    v[1,2]<-v[2,1]<-sum(av*ad)
    v[2,2]<-sum(av^2)
    u<-eigen(v)$vectors[,2]
    c<-sqrt((1+u[2])/2); s<-sign(u[1])*sqrt((1-u[2])/2) 
    for (k in 1:m) {
            ss<-s^2; cc<-c^2; sc<-s*c
            ai<-a[i,,k]; aj<-a[j,,k]
            aii<-a[i,i,k]; ajj<-a[j,j,k]; aij<-a[i,j,k]
            a[i,,k]<-a[,i,k]<-c*ai-s*aj
            a[j,,k]<-a[,j,k]<-s*ai+c*aj
            a[i,j,k]<-a[j,i,k]<-u[1]*(aii-ajj)/2+u[2]*aij
            a[i,i,k]<-aii*cc+ajj*ss-2*sc*aij
            a[j,j,k]<-ajj*cc+aii*ss+2*sc*aij 
            }
    if (vectors) {
        ki<-kk[,i]; kj<-kk[,j]
        kk[,i]<-c*ki-s*kj
        kk[,j]<-s*ki+c*kj
        }
    }
    fnew<-saa-sum(apply(a,3,function(x) sum(diag(x^2))))
    if (verbose)
        cat("Iteration ",formatC(itel,digits=4),"old loss ",formatC(fold,width=10),"new loss ",formatC(fnew,width=10),"\n")
    if (((fold-fnew) < eps) || (itel == itmax)) break()
    itel<-itel+1; fold<-fnew
    }
return(list(a=a,d<-apply(a,3,diag),k=kk))
}

jSimSVD<-function(x,eps=1e-6,itmax=1000,vectors=TRUE,verbose=FALSE) {
n<-dim(x)[1]; m<-dim(x)[2]; nmat<-dim(x)[3]; itel<-1
kkk<-diag(n); lll<-diag(m); sxx<-sum(x^2); fold<-Inf
print(dim(x))
repeat {
    for (i in 1:(n-1)) {
    	if (i > m) next()
		for (j in (i+1):n) {
	        v<-matrix(0,2,2)
	        for (imat in 1:nmat) {
	            xij<-ifelse(j > m,0,x[i,j,imat])
	            xjj<-ifelse(j > m,0,x[j,j,imat])
				xii<-x[i,i,imat]; xji<-x[j,i,imat]
	            v[1,1]<-v[1,1]+(xij^2+xji^2)
	            v[1,2]<-v[2,1]<-v[1,2]+(xii*xji-xjj*xij)
	            v[2,2]<-v[2,2]+(xii^2+xjj^2)
	            }
	        u<-eigen(v)$vectors[,1]
	        for (imat in 1:nmat) {
	            xi<-x[i,,imat]; xj<-x[j,,imat]
	            x[i,,imat]<-u[2]*xi+u[1]*xj
	            x[j,,imat]<-u[2]*xj-u[1]*xi
	            }
	        if (vectors) {
	            ki<-kkk[i,]; kj<-kkk[j,]
	            kkk[i,]<-u[2]*ki+u[1]*kj
	            kkk[j,]<-u[2]*kj-u[1]*ki
	            }
	        }
	}
    ss<-sum(apply(x,3,diag)^2); fnew<-sqrt((sxx-ss)/sxx)
    if (verbose)
        cat(" Left iteration ",formatC(itel,digits=4),"loss ",formatC(fnew,digits=6,width=10),"\n")
    for (k in 1:(m-1)) {
    	if (k > n) next()
		for (l in (k+1):m) {
	        v<-matrix(0,2,2)
	        for (imat in 1:nmat) {
	            xlk<-ifelse(l > n,0,x[l,k,imat])
	            xll<-ifelse(l > n,0,x[l,l,imat])
	            xkl<-x[k,l,imat]; xkk<-x[k,k,imat]
	            v[1,1]<-v[1,1]+(xkl^2+xlk^2)
	            v[1,2]<-v[2,1]<-v[1,2]+(xll*xlk-xkk*xkl)
	            v[2,2]<-v[2,2]+(xkk^2+xll^2)
	            }
	        u<-eigen(v)$vectors[,1]
	        for (imat in 1:nmat) {
	            xk<-x[,k,imat]; xl<-x[,l,imat]
	            x[,k,imat]<-u[2]*xk-u[1]*xl
	            x[,l,imat]<-u[1]*xk+u[2]*xl
	            }
	        if (vectors) {
	            lk<-lll[,k]; ll<-lll[,l]
	            lll[,k]<-u[2]*lk-u[1]*ll
	            lll[,l]<-u[1]*lk+u[2]*ll
	            }
	        }
	}
    ss<-sum(apply(x,3,diag)^2); fnew<-sqrt((sxx-ss)/sxx)
    if (verbose)
        cat("Right iteration ",formatC(itel,digits=4),"loss ",formatC(fnew,digits=6,width=10),"\n")
    if (((fold - fnew) < eps) || (itel == itmax)) break()
    itel<-itel+1; fold<-fnew
    }
return(list(d=apply(x,3,diag),u=t(kkk),v=lll))
}

jTucker3Diag<-function(a,eps=1e-6,itmax=100,vectors=TRUE,verbose=TRUE) {
n<-dim(a)[1]; m<-dim(a)[2]; k<-dim(a)[3]; nmk<-min(n,m,k)
kn<-diag(n); km<-diag(m); kk<-diag(k); ossq<-0; itel<-1
repeat {
	for (i in 1:(n-1)) for (j in (i+1):n) {
		ai<-a[i,,]; aj<-a[j,,]
		acc<-ass<-asc<-0
		if (i <= min(m,k)) {
			acc<-acc+a[i,i,i]^2
			ass<-ass+a[j,i,i]^2
			asc<-asc+a[i,i,i]*a[j,i,i]
			}
		if (j <= min(m,k)) {
			acc<-acc+a[j,j,j]^2
			ass<-ass+a[i,j,j]^2
			asc<-asc-a[j,j,j]*a[i,j,j]
			}
		u<-eigen(matrix(c(acc,asc,asc,ass),2,2))$vectors[,1]
		c<-u[1]; s<-u[2]
		a[i,,]<-c*ai+s*aj
		a[j,,]<-c*aj-s*ai
		if (vectors) {
		    ki<-kn[i,]; kj<-kn[j,]
		    kn[i,]<-c*ki+s*kj
		    kn[j,]<-c*kj-s*ki
			}
		}
	for (i in 1:(m-1)) for (j in (i+1):m) {
		ai<-a[,i,]; aj<-a[,j,]
		acc<-ass<-asc<-0
		if (i <= min(n,k)) {
			acc<-acc+a[i,i,i]^2
			ass<-ass+a[i,j,i]^2
			asc<-asc+a[i,i,i]*a[i,j,i]
			}
		if (j <= min(n,k)) {
			acc<-acc+a[j,j,j]^2
			ass<-ass+a[j,i,j]^2
			asc<-asc-a[j,j,j]*a[j,i,j]
			}
		u<-eigen(matrix(c(acc,asc,asc,ass),2,2))$vectors[,1]
		c<-u[1]; s<-u[2]
		a[,i,]<-c*ai+s*aj
		a[,j,]<-c*aj-s*ai
		if (vectors) {
		    ki<-km[i,]; kj<-km[j,]
		    km[i,]<-c*ki+s*kj
		    km[j,]<-c*kj-s*ki
			}
		}
	for (i in 1:(k-1)) for (j in (i+1):k) {
		ai<-a[,,i]; aj<-a[,,j]
		acc<-ass<-asc<-0
		if (i <= min(n,m)) {
			acc<-acc+a[i,i,i]^2
			ass<-ass+a[i,i,j]^2
			asc<-asc+a[i,i,i]*a[i,i,j]
			}
		if (j <= min(n,m)) {
			acc<-acc+a[j,j,j]^2
			ass<-ass+a[j,j,i]^2
			asc<-asc-a[j,j,j]*a[j,j,i]
			}
		u<-eigen(matrix(c(acc,asc,asc,ass),2,2))$vectors[,1]
		c<-u[1]; s<-u[2]
		a[,,i]<-c*ai+s*aj
		a[,,j]<-c*aj-s*ai
		if (vectors) {
		    ki<-kk[i,]; kj<-kk[j,]
		    kk[i,]<-c*ki+s*kj
		    kk[j,]<-c*kj-s*ki
			}
		}
	nssq<-0; for (v in 1:nmk) nssq<-nssq+a[v,v,v]^2
    if (verbose)
        cat("Iteration ",formatC(itel,digits=4),"ssq ",formatC(nssq,digits=10,width=15),"\n")
    if (((nssq - ossq) < eps) || (itel == itmax)) break()
    itel<-itel+1; ossq<-nssq
    }
d<-rep(0,nmk); for (v in 1:nmk) d[v]<-a[v,v,v]
return(list(a=a,d=d,kn=kn,km=km,kk=kk))
}

jTucker3Block<-function(a,dims,eps=1e-6,itmax=100,vectors=TRUE,verbose=TRUE) {
n<-dim(a)[1]; m<-dim(a)[2]; k<-dim(a)[3]; nmk<-min(n,m,k)
p<-dims[1]; q<-dims[2]; r<-dims[3]
kn<-diag(n); km<-diag(m); kk<-diag(k); ossq<-0; itel<-1
repeat {
	for (i in 1:(n-1)) for (j in (i+1):n) {
		ai<-a[i,,]; aj<-a[j,,]
		acc<-ass<-asc<-0
		if (i <= p) 
			for (u in 1:q) for (v in 1:r) {
				acc<-acc+a[i,u,v]^2
				ass<-ass+a[j,u,v]^2
				asc<-asc+a[i,u,v]*a[j,u,v]
				}
		if (j <= p)
			for (u in 1:q) for (v in 1:r) {
				acc<-acc+a[j,u,v]^2
				ass<-ass+a[i,u,v]^2
				asc<-asc-a[j,u,v]*a[i,u,v]
				}
		u<-eigen(matrix(c(acc,asc,asc,ass),2,2))$vectors[,1]
		c<-u[1]; s<-u[2]
		a[i,,]<-c*ai+s*aj
		a[j,,]<-c*aj-s*ai
		if (vectors) {
		    ki<-kn[i,]; kj<-kn[j,]
		    kn[i,]<-c*ki+s*kj
		    kn[j,]<-c*kj-s*ki
			}
		}
	for (i in 1:(m-1)) for (j in (i+1):m) {
		ai<-a[,i,]; aj<-a[,j,]
		acc<-ass<-asc<-0
		if (i <= q) 
			for (u in 1:p) for (v in 1:r) {
				acc<-acc+a[u,i,v]^2
				ass<-ass+a[u,j,v]^2
				asc<-asc+a[u,i,v]*a[u,j,v]
				}
		if (j <= q)
			for (u in 1:p) for (v in 1:r) {
				acc<-acc+a[u,j,v]^2
				ass<-ass+a[u,i,v]^2
				asc<-asc-a[u,i,v]*a[u,j,v]
				}
		u<-eigen(matrix(c(acc,asc,asc,ass),2,2))$vectors[,1]
		c<-u[1]; s<-u[2]
		a[,i,]<-c*ai+s*aj
		a[,j,]<-c*aj-s*ai
		if (vectors) {
		    ki<-km[i,]; kj<-km[j,]
		    km[i,]<-c*ki+s*kj
		    km[j,]<-c*kj-s*ki
			}
		}
	for (i in 1:(k-1)) for (j in (i+1):k) {
		ai<-a[,,i]; aj<-a[,,j]
		acc<-ass<-asc<-0
		if (i <= r) 
			for (u in 1:p) for (v in 1:q) {
				acc<-acc+a[u,v,i]^2
				ass<-ass+a[u,v,j]^2
				asc<-asc+a[u,v,i]*a[u,v,j]
				}
		if (j <= r)
			for (u in 1:p) for (v in 1:q) {
				acc<-acc+a[u,v,j]^2
				ass<-ass+a[u,v,i]^2
				asc<-asc-a[u,v,i]*a[u,v,j]
				}
		u<-eigen(matrix(c(acc,asc,asc,ass),2,2))$vectors[,1]
		c<-u[1]; s<-u[2]
		a[,,i]<-c*ai+s*aj
		a[,,j]<-c*aj-s*ai
		if (vectors) {
		    ki<-kk[i,]; kj<-kk[j,]
		    kk[i,]<-c*ki+s*kj
		    kk[j,]<-c*kj-s*ki
			}
		}
	nssq<-0; for (i in 1:p) for (j in 1:q) for (l in 1:r) nssq<-nssq+a[i,j,l]^2
    if (verbose)
        cat("Iteration ",formatC(itel,digits=4),"ssq ",formatC(nssq,digits=10,width=15),"\n")
    if (((nssq - ossq) < eps) || (itel == itmax)) break()
    itel<-itel+1; ossq<-nssq
    }
d<-a[1:p,1:q,1:r]
return(list(a=a,d=d,kn=kn,km=km,kk=kk))
}

jMCA<-function(burt,k,eps=1e-6,itmax=500,verbose=TRUE,vectors=TRUE) {
m<-length(k); burt<-m*m*burt/sum(burt); sk<-sum(k)
db<-diag(burt); ll<-kk<-ww<-diag(sk); itel<-1; ossq<-0
klw<-1+cumsum(c(0,k))[1:m]; kup<-cumsum(k)
ind<-lapply(1:m,function(i) klw[i]:kup[i])
sburt<-burt/sqrt(outer(db,db))
for (i in 1:m) 
    kk[ind[[i]],ind[[i]]]<-t(svd(sburt[ind[[i]],])$u)
kbk<-kk%*%sburt%*%t(kk)
for (i in 1:m) for (j in 1:m) 
    ww[ind[[i]],ind[[j]]]<-ifelse(outer(1:k[i],1:k[j],"=="),1,0)
repeat {
    for (l in 1:m) {
        if (k[l] == 2) next()
        li<-ind[[l]]
        for (i in (klw[l]+1):(kup[l]-1)) for (j in (i+1):kup[l]) {
            bi<-kbk[i,-li]; bj<-kbk[j,-li]
            wi<-ww[i,-li]; wj<-ww[j,-li]
            acc<-sum(wi*bi^2)+sum(wj*bj^2)
            acs<-sum((wi-wj)*bi*bj)
            ass<-sum(wi*bj^2)+sum(wj*bi^2)
            u<-eigen(matrix(c(acc,acs,acs,ass),2,2))$vectors[,1]
            c<-u[1]; s<-u[2]
            kbk[-li,i]<-kbk[i,-li]<-c*bi+s*bj
            kbk[-li,j]<-kbk[j,-li]<-c*bj-s*bi
            if (vectors) {
                ki<-kk[i,li]; kj<-kk[j,li]
                kk[i,li]<-c*ki+s*kj
                kk[j,li]<-c*kj-s*ki
                }
            }
        }
    nssq<-sum(ww*kbk^2)-sum(diag(kbk)^2)
    if (verbose)
        cat("Iteration ",formatC(itel,digits=4),"ssq ",formatC(nssq,digits=10,width=15),"\n")
    if (((nssq - ossq) < eps) || (itel == itmax)) break()
    itel<-itel+1; ossq<-nssq
    }
kl<-unlist(sapply(k,function(i) 1:i))
pp<-ifelse(outer(1:sk,order(kl),"=="),1,0)
pkbkp<-t(pp)%*%kbk%*%pp
pk<-t(pp)%*%kk
km<-as.vector(table(kl)); nm<-length(km)
klw<-1+cumsum(c(0,km))[1:nm]; kup<-cumsum(km)
for (i in 1:length(km)) {
	if (km[i]==1) next()
	ind<-klw[i]:kup[i]
	ll[ind,ind]<-eigen(pkbkp[ind,ind])$vectors
	}
lpkbkpl<-t(ll)%*%pkbkp%*%ll
lpk<-t(ll)%*%pk
return(list(kbk=kbk,pkbkp=pkbkp,lpkbkpl=lpkbkpl,kk=t(kk),pp=pp,ll=ll,kpl=t(lpk)))
}
