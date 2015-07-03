interp2D.fun=function(x,y,obj){
    xobj = obj$x
    yobj = obj$y
    zobj = obj$z
    nx = length(xobj)
    ny = length(yobj)
    lx = approx(xobj, 1:nx, x)$y
    ly = approx(yobj, 1:ny, y)$y
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    ey = ly - ly1
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    return(
	zobj[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
	zobj[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
	zobj[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
	zobj[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}

interp3D.fun=function(x,y,z,obj){
    xobj = obj$x
    yobj = obj$y
    zobj = obj$z
    tobj = obj$t
    nx = length(xobj)
    ny = length(yobj)
    nz = length(zobj)
    lx = approx(xobj, 1:nx, x)$y
    ly = approx(yobj, 1:ny, y)$y
    lz = approx(zobj, 1:nz, z)$y
    lx1 = floor(lx)
    ly1 = floor(ly)
    lz1 = floor(lz)
    ex = lx - lx1
    ey = ly - ly1
    ez = lz - lz1
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    ez[lz1 == nz] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    lz1[lz1 == nz] = nz - 1
    return(
	tobj[cbind(lx1, ly1, lz1)] * (1 - ex) * (1 - ey) * (1 - ez) +
	tobj[cbind(lx1 + 1, ly1, lz1)] * ex * (1 - ey) * (1 - ez) +
	tobj[cbind(lx1, ly1 + 1, lz1)] * (1 - ex) * ey * (1 - ez) +
	tobj[cbind(lx1 + 1, ly1 + 1, lz1)] * ex * ey * (1 - ez) +
	tobj[cbind(lx1, ly1, lz1 + 1)] * (1 - ex) * (1 - ey) * ez +
	tobj[cbind(lx1 + 1, ly1, lz1 + 1)] * ex * (1 - ey) * ez +
	tobj[cbind(lx1, ly1 + 1, lz1 + 1)] * (1 - ex) * ey * ez +
	tobj[cbind(lx1 + 1, ly1 + 1, lz1 + 1)] * ex * ey * ez
    )
}

findiff=function(arr,bins){

dimarr=dim(arr)
dimx=dimarr[1]
dimy=dimarr[2]
dimz=dimarr[3]

Mx=1:(dimx-2)
x=2:(dimx-1)
Px=3:dimx
My=1:(dimy-2)
y=2:(dimy-1)
Py=3:dimy
Mz=1:(dimz-2)
z=2:(dimz-1)
Pz=3:dimz

f=array(NA,dim=dimarr)

fx=array(NA,dim=dimarr)
fy=array(NA,dim=dimarr)
fz=array(NA,dim=dimarr)

fxx=array(NA,dim=dimarr)
fyy=array(NA,dim=dimarr)
fzz=array(NA,dim=dimarr)

fxy=array(NA,dim=dimarr)
fxz=array(NA,dim=dimarr)
fyz=array(NA,dim=dimarr)

f[x,y,z]=arr[x,y,z]

fx[x,y,z]=(arr[Px,y,z]-arr[Mx,y,z])/(2*bins[1])
fy[x,y,z]=(arr[x,Py,z]-arr[x,My,z])/(2*bins[2])
fz[x,y,z]=(arr[x,y,Pz]-arr[x,y,Mz])/(2*bins[3])

fxx[x,y,z]=(arr[Px,y,z]-2*arr[x,y,z]+arr[Mx,y,z])/(bins[1]^2)
fyy[x,y,z]=(arr[x,Py,z]-2*arr[x,y,z]+arr[x,My,z])/(bins[2]^2)
fzz[x,y,z]=(arr[x,y,Pz]-2*arr[x,y,z]+arr[x,y,Mz])/(bins[3]^2)

fxy[x,y,z]=(arr[Px,Py,z]-arr[Px,My,z]-arr[Mx,Py,z]+arr[Mx,My,z])/(4*bins[1]*bins[2])
fxz[x,y,z]=(arr[Px,y,Pz]-arr[Px,y,Mz]-arr[Mx,y,Pz]+arr[Mx,y,Mz])/(4*bins[1]*bins[3])
fyz[x,y,z]=(arr[x,Py,Pz]-arr[x,Py,Mz]-arr[x,My,Pz]+arr[x,My,Mz])/(4*bins[2]*bins[3])

return=list(f=f,fx=fx,fy=fy,fz=fz,fxx=fxx,fyy=fyy,fzz=fzz,fxy=fxy,fxz=fxz,fyz=fyz)
}

findeig=function(i,j,k,data){
	temp=matrix(c(data$fxx[i,j,k], data$fxy[i,j,k], data$fxz[i,j,k], data$fxy[i,j,k], data$fyy[i,j,k], data$fyz[i,j,k], data$fxz[i,j,k], data$fyz[i,j,k], data$fzz[i,j,k]),3)
	temp2=eigen(temp)
	return(list(eigval=temp2$values,eigvec=temp2$vectors,hess=temp))
}

alleig=function(data){

dimarr=dim(data$f)
dimx=dimarr[1]
dimy=dimarr[2]
dimz=dimarr[3]

eigarr1=array(NA,dim=dimarr)
eigarr2=array(NA,dim=dimarr)
eigarr3=array(NA,dim=dimarr)

for(i in 2:(dimx-1)){
    for(j in 2:(dimy-1)){
        for(k in 2:(dimz-1)){

            tempeig=findeig(i,j,k,data)
		
            eigarr1[i,j,k]=tempeig$eigval[1]
            eigarr2[i,j,k]=tempeig$eigval[2]
            eigarr3[i,j,k]=tempeig$eigval[3]

        }
    }
}

return=list(eig1=eigarr1,eig2=eigarr2,eig3=eigarr3)

}
