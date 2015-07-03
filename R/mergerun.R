mergerun=function(CenMass=1e12,RelMass=2^(-4),RelRad=2^(-2),PeriToApo=2^(-1),N1=1e3,N2=1e3,C1=10,C2=10,gadgetloc='/Applications/GadgetRun/',gadgetbin='Gadget2noNFW',cpus=4,mpi='/opt/local/bin/openmpirun',outputloc='/Users/aaron/DarkMatter/',soft=1,Tmax=14,Tgap=0.01,Tag='',dosim=TRUE,dopic=TRUE,domovie=TRUE,lim=2,alpha=1/3,pad=4,denylim=c(0,3),bw=0.05,Tconv=0.9778139,...){

name=paste('CM',round(log10(CenMass),2),'RM',round(log2(RelMass),2),'RR',round(log2(RelRad),2),'PA',round(log2(PeriToApo),2),Tag,sep='')

Rvir=MvirToRvir(CenMass)

if(dosim){

Tconv=0.9778139

dir.create(paste(outputloc,name,'/',sep=''),showWarnings=F)
if(dopic){dir.create(paste(outputloc,name,'/pic/',sep=''),showWarnings=F)}

Rad=Rvir*RelRad

Vel=VisViva(CenMass,Rad=Rad,PeriToApo=PeriToApo)

hern1=genhalo(Npart=N1,Mpart=CenMass/N1,conc=C1,Munit=1,Lunit=1e3)
hern2=genhalo(x=Rad,vy=Vel,Npart=N2,Mpart=RelMass*CenMass/N2,conc=C2,Munit=1,Lunit=1e3,IDoffset=max(hern1$part[,1]))
herncomb=rbind(hern1$part,hern2$part)
herncomb=addhead(part=herncomb)
snapwrite(part=herncomb$part,head=herncomb$head,file=paste(outputloc,name,'/hernmerg.gdt',sep=''))

genparam(
ParamBase = paste(outputloc,'/',name,sep=''),
InitCondFile = paste(outputloc,'/',name,'/hernmerg.gdt',sep=''),
OutputDir = paste(outputloc,'/',name,sep=''),
TimeMax = Tmax,
TimeBetSnapshot = Tgap,
UnitMass_in_g = 1.989e+33,
SofteningGas = soft,
SofteningHalo = soft,
SofteningGasMaxPhys = soft,
SofteningHaloMaxPhys = soft,
...
)

system(command=paste(mpi,' -np ',cpus,' ',gadgetloc,'/',gadgetbin,' ',outputloc,'/',name,'/galaxy.param',sep=''))

}

tempfiles=list.files(paste(outputloc,'/',name,'/',sep=''))
tempfiles=tempfiles[grep('snapshot',tempfiles)]
resort={}
for(i in 1:length(tempfiles)){resort=c(resort,as.numeric(strsplit(tempfiles[i],'_')[[1]][2]))}
tempfiles=tempfiles[order(resort)]

if(dopic){

if(require(doParallel)==FALSE){stop('Need package doParallel')}
if(require(foreach)==FALSE){stop('Need package foreach')}
if(require(Cairo)==FALSE){stop('Need package Cairo')}

registerDoParallel(cores=6)

boundIDs=snapboundR(paste(outputloc,name,'/',formatC(tempfiles[1],flag=0,width=3),sep=''),maxN=2e4)

seps={}
for(i in 1:length(tempfiles)){
temp=snapread(paste(outputloc,name,'/',formatC(tempfiles[i],flag=0,width=3),sep=''))
time=temp$head$Time*Tconv
rsep=(1/Rvir)*sqrt((temp$part[temp$part[,'ID']==boundIDs[1],'x']-temp$part[temp$part[,'ID']==boundIDs[2],'x'])^2+(temp$part[temp$part[,'ID']==boundIDs[1],'y']-temp$part[temp$part[,'ID']==boundIDs[2],'y'])^2+(temp$part[temp$part[,'ID']==boundIDs[1],'z']-temp$part[temp$part[,'ID']==boundIDs[2],'z'])^2)
vsep=(Tconv/Rvir)*sqrt((temp$part[temp$part[,'ID']==boundIDs[1],'vx']-temp$part[temp$part[,'ID']==boundIDs[2],'vx'])^2+(temp$part[temp$part[,'ID']==boundIDs[1],'vy']-temp$part[temp$part[,'ID']==boundIDs[2],'vy'])^2+(temp$part[temp$part[,'ID']==boundIDs[1],'vz']-temp$part[temp$part[,'ID']==boundIDs[2],'vz'])^2)
seps=rbind(seps,c(time,rsep,vsep))
}

foreach(i = 1:length(tempfiles)) %dopar%{

snapim(snapfile=paste(outputloc,name,'/',formatC(tempfiles[i],flag=0,width=3),sep=''),lim=lim,outfile=paste(outputloc,name,'/pic/pic',formatC(i,flag=0,width=pad),'.png',sep=''),alpha=alpha,denylim=denylim,bw=bw,dobound=TRUE,CenBound=boundIDs[1],SatBound=boundIDs[2],addbound=T,rseppast=rbind(seps[1:i,c(1,2)]),vseppast=rbind(seps[1:i,c(1,3)]))

#temp=snapread(paste(outputloc,name,'/',formatC(tempfiles[i],flag=0,width=3),sep=''))
#png(file=paste(outputloc,name,'/pic/pic',formatC(tempfiles[i],flag=0,width=3),'.png',sep=''),width=400,height=400)
#plot.new()
#par(mar=c(0,0,0,0))
#plot.window(xlim=c(-300,300), ylim=c(-300,300),asp=1)
#points(temp$part[,2:3], pch='.', col=hsv(magmap(temp$part[,1])$map, alpha=1/3),cex=3)
#draw.circle(0,0,radius=c(100,200,300),lty=3)
#abline(h=0,lty=2)
#abline(v=0,lty=2)
#temptime=as.numeric(strsplit(tempfiles[i],'snapshot_')[[1]][2])
#legend('topleft',legend=paste('Time =',round(as.numeric(temptime)*Tgap*Tconv,4),'Gyrs'))
#dev.off()

}

if(domovie){
invisible(suppressWarnings(file.remove(paste(outputloc,name,'/pic/merge.mov',sep=''))))
MakeMovie(imdir=paste(outputloc,name,'/pic/',sep=''),imhead='pic',pad=pad,bitrate='1280k',outmov='merge.mov')
}

}

}

