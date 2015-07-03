snapboundR=function(snapfile,maxN=1e4,soft=0){
temp=snapread(snapfile)
CenHalo=which(temp$part[,'ID']<length(temp$part[,'ID'])/2)
SatHalo=which(temp$part[,'ID']>length(temp$part[,'ID'])/2)
if(length(CenHalo)>maxN){CenHalo=sample(CenHalo,maxN)}
if(length(SatHalo)>maxN){SatHalo=sample(SatHalo,maxN)}
tempCen=orbit(pos=temp$part[CenHalo,c('x','y','z')],vel=temp$part[CenHalo,c('vx','vy','vz')],mass=temp$part[CenHalo,'Mass'],Lunit=1e3,Munit=1,Nstep=0,soft=soft)
CenID=temp$part[CenHalo[which.min(tempCen$gpe/temp$part[CenHalo,'Mass']+tempCen$ke/temp$part[CenHalo,'Mass'])],'ID']
tempSat=orbit(pos=temp$part[SatHalo,c('x','y','z')],vel=temp$part[SatHalo,c('vx','vy','vz')],mass=temp$part[SatHalo,'Mass'],Lunit=1e3,Munit=1,Nstep=0,soft=soft)
SatID=temp$part[SatHalo[which.min(tempSat$gpe/temp$part[SatHalo,'Mass']+tempSat$ke/temp$part[SatHalo,'Mass'])],'ID']
return(c(CenID=CenID,SatID=SatID))
}

snappot=function(snap,grid=seq(-300,300,by=10),doeig=FALSE,doplot=TRUE){

bins=grid[2]-grid[1]
snappotden=potden(pos=snap$part[,2:4], mass=snap$part[,8], grid=grid, bins=bins)

if(doeig){
	snapfindiff=findiff(snappotden$arrpot$t,bins=snappotden$bins)
	snapeig=alleig(snapfindiff)
	allpos=which(snapeig$eig1>0 & snapeig$eig2>0 & snapeig$eig3>0 ,arr.ind=T)
	
	if(doplot){
		image(snappotden$arrpot$x,snappotden$arrpot$y,magmap(snapeig$eig1[,,ceiling(length(grid)/2)],log=T,lo=0.1,hi=1,flip=T)$map,col=rainbow(1e4,end=2/3),asp=1)
		points(snap$part[snap$part[,4]> -bins*1.5 & snap$part[,4]<bins*1.5,2:3],pch='.',col=hsv(v=0,alpha=0.2))
	}
temp=list(snap=snap, potden=snappotden, diff= snapfindiff, eig= snapeig, allpos=allpos)
}else{temp=list(snap=snap, potden=snappotden)}
return=temp
}

snapbound=function(snapfile){
temp=snapread(snapfile)
temp=addenergy(temp)
CenHalo=temp$part[,'ID']<length(temp$part[,'ID'])/2
SatHalo=temp$part[,'ID']>length(temp$part[,'ID'])/2

CenID=temp$part[which(CenHalo)[which.min(temp$part[CenHalo,'gpe']/temp$part[CenHalo,'Mass']+temp$part[CenHalo,'ke']/temp$part[CenHalo,'Mass'])],'ID']
SatID=temp$part[which(SatHalo)[which.min(temp$part[SatHalo,'gpe']/temp$part[SatHalo,'Mass']+temp$part[SatHalo,'ke']/temp$part[SatHalo,'Mass'])],'ID']

return(c(CenID=CenID,SatID=SatID))
}
