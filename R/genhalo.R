genhalo <-
function(x=0,y=0,z=0,vx=0,vy=0,vz=0,Npart=1000,Mpart=1,conc=10,iseed=87618820,Munit=1e10,Lunit=1e3,Vunit=1,IDoffset=0){
if(length(x)==6){y=x[2];z=x[3];vx=x[4];vy=x[5];vz=x[6];x=x[1]}

out=.Fortran("genhalofortran",
npart=as.integer(Npart),
x=as.double(rep(NA,2*Npart)),
y=as.double(rep(NA,2*Npart)),
z=as.double(rep(NA,2*Npart)),
vx=as.double(rep(NA,2*Npart)),
vy=as.double(rep(NA,2*Npart)),
vz=as.double(rep(NA,2*Npart)),
mpart=as.double(Mpart*Munit/1e10),
sigma=as.double(NA),
conc=as.double(conc),
iseed=as.integer(iseed),
NAOK=TRUE,PACKAGE='simpack')

select= ! is.na(out$x)
Ntot=length(which(select))
return=list(part=data.frame(ID=1:Ntot+IDoffset,x=(out$x[select]*1e6/Lunit)+x, y=(out$y[select]*1e6/Lunit)+y, z=(out$z[select]*1e6/Lunit)+z, vx=(out$vx[select]+vx)*1/Vunit, vy=(out$vy[select]+vy)*1/Vunit, vz=(out$vz[select]+vz)*1/Vunit, Mass=out$mpart*1e10/Munit), params=c(Npart=Npart, Mpart=out$mpart*1e10, Mvir=Npart*out$mpart*1e10, Mtot=Ntot*out$mpart*1e10, sigma=out$sigma, conc=out$conc, iseed=iseed))
}
