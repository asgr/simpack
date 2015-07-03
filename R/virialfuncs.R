MvirToSigma=function(Mvir=1e12,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1){
H0=h*100
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
rho_crit=(3*H0^2*(kms_to_ms/(1e6*pc_to_m))^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit)
return=(Mvir*sqrt(32*pi*g^3*DeltaVir*rho_crit/3))^(1/3)
}

SigmaToMvir=function(Sigma=230,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1){
H0=h*100
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return=Sigma^3*sqrt(3/(32*pi*g^3*DeltaVir*rho_crit))
}

MvirToRvir=function(Mvir=1e12,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1){
H0=h*100
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return=((3*Mvir)/(4*pi*DeltaVir*rho_crit))^(1/3)
}

RvirToMvir=function(Rvir=162.635,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1){
H0=h*100
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return=(4*pi/3)*DeltaVir*rho_crit*Rvir^3
}

SigmaToRvir=function(Sigma=230,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1){
H0=h*100
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return=Sigma*(27/(512*pi^3*g^3*DeltaVir^3*rho_crit^3))^(1/6)
}

RvirToSigma=function(Rvir=162.635,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1){
H0=h*100
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return=Rvir*((512*pi^3*g^3*DeltaVir^3*rho_crit^3)/27)^(1/6)
}

Mpart=function(VirialParam=1e12,Npart=1000,type='mass',Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1,int=FALSE){
if(type=='mass'){temp=VirialParam/Npart}
if(type=='sig'){temp=SigmaToMvir(VirialParam,Munit=Munit,Lunit=Lunit,Vunit=Vunit,DeltaVir=DeltaVir,h=h)/Npart}
if(type=='rad'){temp=RvirToMvir(VirialParam,Munit=Munit,Lunit=Lunit,Vunit=Vunit,DeltaVir=DeltaVir,h=h)/Npart}
if(int){temp=as.integer(temp)}
return=temp
}

Npart=function(VirialParam=1e12,Mpart=1e9,type='mass',Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1,int=FALSE){
if(type=='mass'){temp=VirialParam/Mpart}
if(type=='sig'){temp=SigmaToMvir(VirialParam,Munit=Munit,Lunit=Lunit,Vunit=Vunit,DeltaVir=DeltaVir,h=h)/Mpart}
if(type=='rad'){temp=RvirToMvir(VirialParam,Munit=Munit,Lunit=Lunit,Vunit=Vunit,DeltaVir=DeltaVir,h=h)/Mpart}
if(int){temp=as.integer(temp)}
return=temp
}

FreeFall=function(M1=1e12,M2=1,Rad=162.635,Munit=1,Lunit=1e3,Vunit=1){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
return=(pi/sqrt(g*(M1+M2)))*(Rad/2)^(3/2)
}

VisViva=function(M=1e12,Rad=162.635,PeriToApo=1,Munit=1,Lunit=1e3,Vunit=1){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
RadMaj=Rad*(1+1/PeriToApo)/2
return=sqrt(g*M*(2/Rad-1/RadMaj))
}

RocheRad=function(M1=1e12,M2=1e10,Size=35.03865,Rfac=2.44){
return=Rfac*Size*(M1/M2)^(1/3)
}

RocheSize=function(M1=1e12,M2=1e10,Rad=396.8294,Rfac=2.44){
return=(Rad/Rfac)*(M2/M1)^(1/3)
}

