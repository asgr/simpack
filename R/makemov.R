MakeMovie=function(imdir,imhead,outmov='out.mov',ffmpeg='/usr/bin/ffmpeg',framerate=30,bitrate='640k',pad=3,type='png'){
system(paste(ffmpeg,' -r ',framerate,' -i ',imdir,'/',imhead,'%0',pad,'d.',type,' -c:v libx264 -b:v ',bitrate,' -pix_fmt yuv420p ',imdir,'/',outmov,sep=''))
}

