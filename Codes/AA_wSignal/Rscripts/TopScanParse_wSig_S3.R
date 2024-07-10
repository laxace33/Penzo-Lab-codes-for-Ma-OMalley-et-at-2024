root="C:/Users/Jun/Documents/R/AA_wSignal/"
max_cue_dur=15 #in seconds

require(readxl)
require("changepoint")


funcdir=paste(root,"track_func/",sep="") # r function folder
tdir=paste(root,"rawtrack",sep="")
tfiles=list.files(tdir)

tdirttl=paste(root,"rawtrack_ttl",sep="")
ttlfl=list.files(tdirttl)

formatch=gsub(".csv","",ttlfl)
rout=paste(root,"rtrack",sep="")

load (paste(funcdir,"cleantrack.rFunc",sep=""))

for (i in (1:length(tfiles))) ({

tfname=paste(tdir,tfiles[i],sep="/")

gethead=read.table (tfname, sep = '\t',nrows=50, header = F,strip.white=F,blank.lines.skip=F,fill=F,stringsAsFactors=F)
startdat=grep("Format:",gethead[,1])
#track=read.table (tfname,skip=startdat-4, sep = '\t', header = F,strip.white=F,blank.lines.skip=F,fill=F,stringsAsFactors=F)
trackidx=read.table (tfname, sep = '\t',skip=startdat, header = F,strip.white=F,blank.lines.skip=F,fill=F,stringsAsFactors=F)

realstart=startdat-(trackidx[1,1]-1)


track=try (read.table (tfname,skip=realstart, sep = '\t', header = F,strip.white=F,blank.lines.skip=T,fill=F,stringsAsFactors=F),silent=T)

if (is.null (dim (str(track)))==T)({
track=read.table (tfname,skip=startdat-6, sep = '\t', header = F,strip.white=F,blank.lines.skip=T,fill=F,stringsAsFactors=F)
})

if (track[1,1]=="Format:FrameNum") ({
#colnames(track)=track[1,]
track=track[-1,]
row.names(track)=seq(1,nrow(track))
})


forcolmatch=c("Format:FrameNum", "CenterX(mm)","CenterY(mm)","NoseX(mm)","NoseY(mm)","EventRule1","EventRule2")
namematch=gethead [seq (startdat,max (grep("EventRule",gethead[,1]))),1]

track=track[,match(forcolmatch,namematch)]
colnames(track)=c("time","posx","posy","xnose","ynose","left","right") #,"NA1","NA2","NA3")

track$fnum=track$time


getinf=read.table (tfname,nrows=22, sep = '\t', header = F,,strip.white=T,blank.lines.skip=T,fill=F,stringsAsFactors=F)
fps=as.numeric (strsplit (getinf[7,],":")[[1]][2])
starttime=strsplit (getinf[10,],":")[[1]][2]
stime=gsub("[[:space:]]", "",starttime)
stime=as.numeric (gsub("\\(s)", "",stime))

splitinf=strsplit (getinf[1,],"\\\\")[[1]]
inf=gsub (".BMP","",splitinf[length(splitinf)])
inf=gsub(".wmv","",inf)
inf=gsub(".avi","",inf)

track$time=as.numeric (track$time) /fps

dat_ttl=read.csv (paste (tdirttl,ttlfl [match(inf,formatch)],sep="/"))
dat_ttl=as.data.frame(dat_ttl,stringsAsFactors=F)

noshock_idx=match("TTL.to.opto.active",colnames(dat_ttl))

if (is.na(noshock_idx)==F) ({
dat_ttl$shock.to.photo.active=0
getshock=rle(dat_ttl$cues.to.photo.active)
forshock=cumsum(getshock$lengths)
endcue=forshock+1
onlycues=which(getshock$values==1)
cuestop=dat_ttl$Time [endcue[onlycues]]
cuestart=dat_ttl$Time [forshock[onlycues-1] +1]
cuedur=cuestop-cuestart
strials=which(cuedur>=max_cue_dur)

if (length(strials)>0) ({
dat_ttl$shock.to.photo.active [match (cuestop[strials],dat_ttl$Time)]=1
})

dat_ttl=dat_ttl[,-noshock_idx]

})

cidx=c(grep ("cues.to.photo.active",colnames(dat_ttl)),grep ("shock.to.photo.active",colnames(dat_ttl)))

dat_ttl=dat_ttl[,c(1,cidx)]

colnames(dat_ttl)=c("time","cues","shock")


getframe.rfunc=function(x,y,s) ({  #x==dat_ttl, y==track
cmatch=grep(s,colnames(x))
rlecues=rle(x[,cmatch])
cumlength=cumsum(rlecues$lengths)
idx=which(rlecues$values==1)
starts=cumlength[idx-1] +1
stops=cumlength[idx]

forbind=cbind(x$time [starts],x$time [stops])
out=unlist (apply(forbind,1,function(p) (which(y$time>=p[1]) [which(y$time>=p[1]) %in% which(y$time<=p[2])])))

return(out)
})


track$cues=0
track$cues [getframe.rfunc(dat_ttl,track,"cues")]=1
track$shock=0
track$shock [getframe.rfunc(dat_ttl,track,"shock")]=1

track=cleantrack.rfunc(track,"center")
track=cleantrack.rfunc(track,"nose")

dat=track

save(dat,file=paste(rout,paste(inf,".rdat",sep=""),sep="/"))

})
























