require(readxl)
require("changepoint")

root="//eu.boehringer.com/users/bib/users1/duhoffma/Desktop/AA_pipeline/aa_test/" 

freeze_rule=1 #length in seconds for a bout of immobility to trigger freezing
itilength=45 #in seconds
max_shock_time=15

tdir=paste(root,"rtrack",sep="")
tfiles=list.files(tdir)

csvdir=paste(root,"csvout",sep="")
pointer=read.csv (paste(root,"pointer/pointer.csv",sep=""),stringsAsFactors=F)

### this code may be specific to Jun's experiments but I left them here just in case
pointer$day=as.vector (unlist (sapply (pointer$Experiment_Subgroup,function(x) (substr(x,nchar(x),nchar(x))))))
pointer$stype=rep(NA,nrow(pointer))
pointer$stype [grep("Train",pointer$Experiment_Subgroup)]="train"
pointer$stype [grep("Test",pointer$Experiment_Subgroup)]="test"
pointer$stype=paste (pointer$stype,pointer$day,sep="_")
pointer$stype=paste (pointer$stype,tolower (pointer$Light_Status),sep="_")

#pointer$stype=paste (pointer$stype,unlist (lapply (strsplit (pointer$Experiment_Subgroup,"_"),function(x) (paste(x[1],x[2],sep="")))),sep="_")

##################################################

arena=c(21,47,335,178)    #xleft, ybottom, xright, ytop

cueout=list()
fout=list()


for (i in (1:length(tfiles))) ({

load (paste(tdir,tfiles[i],sep="/"))
dat=as.data.frame(dat,stringsAsFactors=F)

dat$posx=as.numeric(dat$posx)
dat$posy=as.numeric(dat$posy)
dat$xnose=as.numeric(dat$xnose)
dat$ynose=as.numeric(dat$ynose)
dat$left=as.numeric(dat$left)
dat$ynose=as.numeric(dat$right)
dat$fnum=as.numeric(dat$fnum)

dat$otime=dat$time
fps=c(0,diff(dat$time))


inf=gsub(".rdat","",tfiles[i])
splitinf=strsplit(inf," ")[[1]]
didx=splitinf[1]

dates=as.Date (didx,format="%m%d%Y")
day=splitinf[3]
id=paste(splitinf[4],splitinf[5],sep="_")

infout=data.frame (cbind(id,dates,day),stringsAsFactors=F)

infout$dates=as.Date (didx,format="%m%d%Y")

dat$zone="NA"
dat$event="NA"

dat$zone[which(dat$right==1)]="right"
dat$zone[which(dat$left==1)]="left"
dat$event[which(dat$cues==1)]="cue"

dat$event[which(dat$shock==1)]=NA
tmatch=gsub (" ","",gsub(".rdat","",tfiles[i]))


dat$zchange="NA"
zonerle=rle(dat$zone)
zchange=cumsum (zonerle$lengths)+1

if(max(zchange)>nrow(dat)) ({
zchange=zchange[-length(zchange)]
})

dat$zchange[zchange]="yes"
cuerle=rle(dat$event)
cuesum=cumsum(cuerle$lengths)
cuestarts=cuesum [which(cuerle$values=="cue")-1] +1
cuestops=cuesum [which(cuerle$values=="cue")] 

cuedfn=data.frame (cbind(cuestarts,cuestops),stringsAsFactors=F)

cuedf=data.frame (apply(cuedfn,2,function(x) (dat$time[x])),stringsAsFactors=F)
colnames(cuedf)=c("start","stop")

cuedf$start=round(cuedf$start,digits=3)
cuedf$stop=round(cuedf$stop,digits=3)

cuedf$f_start=cuedfn[,1]
cuedf$f_stop=cuedfn[,2]
cuedf$dur=cuedf[,2]-cuedf[,1]
cuedf$shock="no"
cuedf$shock_start=NA
cuedf$shock_end=NA

cuedf$f_shock_start=NA
cuedf$f_shock_end=NA
cuedf$shock_dur=NA

if (length(which(cuedf$dur>14.8))>0) ({
addtoshock=which(cuedf$dur>14.8)
cuedf$shock[addtoshock]="shock"
addtoshocki=dat[cuedf$f_stop[addtoshock],]
fornewshock=sapply (addtoshocki$fnum,function(x) (which(dat$fnum>x) [which (which(dat$fnum>x) %in% which(dat$zchange=="yes"))][1]))

getnoesc=dat$time[fornewshock]-addtoshocki$time

if (length(getnoesc>max_shock_time)>0) ({
getnoesci=which(getnoesc>max_shock_time)
if (length (getnoesci)!=0) ({
fornewshock[getnoesci]=findInterval (addtoshocki$time[getnoesci]+max_shock_time,dat$time)
})

if (length (getnoesci)==0) ({
fornewshock=findInterval (addtoshocki$time+max_shock_time,dat$time)
})

if (length(which(is.na(fornewshock)))>0) ({
fornewshock=fornewshock [-which(is.na(fornewshock))]})
dat$event[unlist (sapply(seq(1,length(fornewshock)),function(x) (seq(addtoshocki$fnum[x]+1,fornewshock[x]))))]="shock"
})
})

nidx=as.list (seq(1,nrow(cuedf)))
newcue=unlist (lapply(nidx,function(x) (which(dat$time>=cuedf[x,1]) [which(dat$time>=cuedf[x,1]) %in% which (dat$time<=cuedf[x,2])])))


dat$event[newcue]="cue"

allshock=which(dat$event=="shock")

if (length(allshock)>0) ({
strials=findInterval(dat$time[allshock],cuedf$start)
cuedf$shock [unique(strials)]="yes"

sidx=lapply (as.list (unique(strials)),function(x) (which(strials==x)))
stimes=lapply(sidx,function(x) (dat$time[allshock][x]))
reals=do.call ("rbind",lapply(stimes,function(x) (cbind (x[1],x[length(x)]))))
reals=round(reals,digits=3)

sframes=lapply(sidx,function(x) (allshock[x]))
sframe=do.call ("rbind",lapply(sframes,function(x) (cbind (x[1],x[length(x)]))))

cuedf$shock_start [which(cuedf$shock=="yes")]=reals[,1]
cuedf$shock_end [which(cuedf$shock=="yes")] =reals[,2]

cuedf$f_shock_start [which(cuedf$shock=="yes")]=sframe[,1]
cuedf$f_shock_end [which(cuedf$shock=="yes")] =sframe[,2]

cuedf$shock_dur [which(cuedf$shock=="yes")]=as.numeric (round (reals[,2]-reals[,1],digits=2))

})


dat$tint=findInterval(dat$time,cuedf$start)
cuedf$tint=seq(1,nrow(cuedf))

dists=unlist (lapply (as.list (seq(1,nrow(dat)-1)),function(x) (sqrt ( ((dat$posx[x+1]-dat$posx[x])^2)+((dat$posy[x+1]-dat$posy[x])^2)))))
dat$dists=c(NA,dists)

dists_nose=unlist (lapply (as.list (seq(1,nrow(dat)-1)),function(x) (sqrt ( ((dat$xnose[x+1]-dat$xnose[x])^2)+((dat$ynose[x+1]-dat$ynose[x])^2)))))
dat$dists_nose=c(NA,dists_nose)

zdivide=mean(c(arena[1],arena[3]))
dat$edist=unlist (lapply (as.list (seq(1,nrow(dat))),function(x) (sqrt ( ((dat$posx[x]-zdivide)^2)+((dat$posy[x]-dat$posy[x])^2)))))


dat$esc=NA
eidx=match (cuedf$tint,dat$tint[which(dat$zchange=="yes")])

dat$esc[match (dat$time[which(dat$zchange=="yes")][eidx],dat$time)]="esc"
findmin=lapply (as.list (which(dat$esc=="esc")),function(x) (dat[seq((x-10),x),]))
newesc=unlist (lapply(findmin,function(x) (x$fnum [which.min(x$edist)])))
allz=newesc

noesc=which (is.na (dat$esc[which(dat$zchange=="yes")]))
if (length(noesc)>0) ({
allz=c(noesc,allz)
})

changez=unlist (lapply(findmin,function(x) (x$fnum [seq(which.min(x$edist)+1,length(x$edist))])))

dat$zchange [which(dat$zchange=="yes")]=NA
dat$esc [which(dat$esc=="esc")]=NA

dat$zchange [allz]="yes"
dat$esc [newesc]="esc"

cuedf$esc_start=NA
esc_idx=findInterval(round (dat$time[which(dat$esc=="esc")],digits=3),cuedf$start)
cuedf$esc_start[esc_idx]=round (dat$time[which(dat$esc=="esc")],digits=3)
cuedf$esc_lat=cuedf$esc_start-cuedf$start
cuedf$esc_frame=NA
cuedf$esc_frame[esc_idx]=which(dat$esc=="esc")



dat=dat[,-match(c("right","left","cues","shock"),colnames(dat))]

ends=unlist (lapply(findmin,function(x) (x$tint[1])))
splitends=split(dat,dat$tint)

forcheck=splitends[match (ends,as.numeric (names(splitends)))]
endcues=as.vector (unlist (lapply(forcheck,function(x) (x$fnum[nrow(x)]))))


if (is.null(endcues)!=T) ({
nocues=unlist (lapply (as.list (seq(1,length(endcues))),function(x) (seq(newesc[x],endcues[x]))))
})

dat$event[nocues [which (dat$event[nocues]=="cue")] ]="NA"

splitcue=split(dat,dat$tint)
firstcues=as.numeric (unlist (lapply(splitcue,function(x) (row.names(x[nrow(x),])))))
dat$tint[firstcues]=dat$tint[firstcues]+1

getcp.rfunc=function(x,pv) ({
penval=as.numeric (pv)

realindex=row.names(x)
realindex=as.numeric(realindex)

if (nrow(x)>3) ({
cptan=cpt.meanvar(x$edist,method="PELT",penalty="CROPS",minseglen=3,Q=100,pen.value=penval)

pbind=pen.value.full(cptan)
sbind=cpts.full (cptan)

if(nrow(sbind)!=length(pbind)) ({
pbind=pbind[1]
})

cpsum=data.frame (cbind (pbind,sbind),stringsAsFactors=F)
colnames(cpsum)=c("pen",paste("cp_",seq(1,ncol(cpsum)-1),sep=""))

######this removes any NA rows from output of change point analysis
if(is.na (cpsum[nrow(cpsum),2])==T) ({
cpsum=cpsum[-nrow(cpsum),]
})

mdatrm.rfunc=function(y) ({
isna=which (is.na(y))
if (length(isna) !=0) ({
y=y [-which (is.na(y))]
})
return(y)
})

#################this looks for large changes in penalty that flag overly strict analyses
pdiff=diff (cpsum$pen)
penthresh=10

if (length(which(pdiff>penthresh))>0) ({
idx=which(pdiff>penthresh)
idx=idx[1]
cpts=cpsum[idx,]

cpts=mdatrm.rfunc(cpts)

}) else({
cpts=cpsum[nrow(cpsum),]

})

cpts=mdatrm.rfunc(cpts)

if(nrow(cpts)==2) ({
cpts=cpsum[nrow(cpsum)-1,]
cpts=mdatrm.rfunc(cpts)
})

cpidx=as.vector (unlist (cpts[,-1]))
cpidx=(cpidx-1)+realindex[1]

out=list()
out[[1]]=cpts[,1]
out[[2]]=cpidx

}) else ({

out=list()
out[[1]]=10
out[[2]]=as.numeric (row.names(x) [nrow(x)-1])
})


return(out)

})

onlyt=dat[which (dat$tint %in% cuedf$tint),]
splitt=split(onlyt,onlyt$tint)

splitt=splitt[esc_idx]

maxsess=dat$time[nrow(dat)]

if((maxsess-cuedf$start[nrow(cuedf)])<18) ({
splitt=splitt[-length(splitt)]
esc_idx=esc_idx[-length(esc_idx)]
})

if (i==345) ({
splitt=splitt[-length(splitt)]
esc_idx=esc_idx[-length(esc_idx)]
})


newt=lapply(splitt,function(x) (x[1:which(x$esc=="esc")[1],]))
ichange=lapply (newt,function(q) (getcp.rfunc(q,c(10,30))))

cps=as.vector (unlist (lapply(ichange,function(x) (x[[2]]))))

esc=match(cuedf$esc_start,round (dat$time,digits=3))[esc_idx]
idx=lapply (as.list (esc),function(x) (cps-x))
estarts=unlist (lapply (lapply(idx,function(x) (which(x<0))),function(q) (q[length(q)])))

dat$e_loco=NA
dat$e_loco[cps[estarts]]="start"

cuedf$f_emov=NA
cuedf$emov=NA
cuedf$emov_cuelat=NA


cuedf$f_emov [dat$tint [cps[estarts]]]=cps[estarts]
cuedf$emov [dat$tint [cps[estarts]]]=round (dat$time [cps[estarts]],digits=3)

cuedf$emov_cuelat=round (cuedf$emov-cuedf$start,digits=3)
cuedf$zone=dat$zone [cuedf$f_start]

nose_thresh=as.numeric (quantile (dat$dists_nose,probs=0.6,na.rm=T))
body_thresh=as.numeric (quantile (dat$dists,probs=0.6,na.rm=T))

test=which (dat$dists <= body_thresh) [which (dat$dists <= body_thresh) %in%  which (dat$dists_nose<=nose_thresh)]
newtest=rep ("no",length (dat$dists))

newtest[test]="yes"
newrle=rle(newtest)

f_thresh=as.numeric (quantile (rle(newtest)$lengths [which (rle(newtest)$values=="yes")],probs=0.995))

realframesper=round (1/fps[10])

if (f_thresh> (3*realframesper)) ({
f_thresh=3*realframesper
})


newstop=cumsum(newrle$lengths)
newstart=(newstop-newrle$lengths)+1

forr=data.frame (cbind (newstart,newstop),stringsAsFactors=F)
colnames(forr)=c("start","stop")
forr$lengths=newrle$lengths
forr$vals=newrle$values
forr$thresh=NA
forr$freeze=NA
forr$freeze_end=NA
forr$freeze_end_real=NA

forr$freeze[which(forr$lengths>=f_thresh)]="yes"

forr$freeze[which(forr$vals=="no")]="no"

yandn=newtest
yandn [which(yandn=="yes")]=1
yandn [which(yandn=="no")]=0
yandn=as.numeric(yandn)


cptan=cpt.meanvar(yandn,method="PELT",penalty="CROPS",minseglen=30,Q=100,pen.value=c(10,30))
pbind=pen.value.full(cptan)
sbind=cpts.full (cptan)

forr$cumsum=cumsum(forr$lengths)
forr$change=NA
forr$change [match (unique (as.vector(sbind)),forr$cumsum)]="yes"



if(is.na(forr$change[1])==F) ({
forr$change[1]="no"
})


diffidx=c(1,which(forr$change=="yes"),nrow(forr)+1)
realint=sapply(seq(1,length(diffidx)-1),function(p) (seq(diffidx[p],diffidx[p+1]-1)))

ys=unlist (lapply(realint,function(p) (sum (forr$lengths[p] [which(forr$vals[p]=="yes")]))))
ns=unlist (lapply(realint,function(p) (sum (forr$lengths[p] [which(forr$vals[p]=="no")]))))

forr$realint=unlist (sapply (seq(1,length(realint)),function(p) (rep(p,length(realint[[p]])))))

forr$yn=NA
forr$yn[match(unique(forr$realint),forr$realint)]=paste(ys,ns,sep="_")

splitc=split(forr,forr$realint)

getit.rfunc=function(x,y,z) ({
getit=as.numeric (strsplit (x$yn[1],"_")[[1]])
out=vector()

if (getit[1]>y) ({
if (getit[2]<y) ({
newthresh=sum(x$lengths[which(x$vals=="yes")] [which (which(x$vals=="yes")%in%  which(x$lengths<=2))])
if (getit[2]*3<(getit[1]-newthresh)) ({
starts=x$start[1]
stops=x$stop[nrow(x)]

if (x$vals[1]=="no") (starts=x$start[2])
if (x$vals[nrow(x)]=="no") (stops=x$stop[nrow(x)-1])

oyes=x[which(x$vals=="yes"),]
if (median(oyes$lengths)<5) ({
newstop=match (oyes$start[max (which(oyes$lengths>=5))],x$start)
stops=x$stop[newstop]
})


secs=c(z$time[starts],z$time[stops])
mins= secs/60
startst <- as.POSIXct(60 * mins, origin = "1970-01-01", tz = "UTC") # POSIXct
#return(as.character (starts))
out=c(starts,stops,as.character (startst))

return(out)
})
})

}) 


if (nrow(x)==1 & x$vals[1]=="no") ({
out=rep(NA,4)
return(out)
})
 

if (length(out)==0) ({
out=rep(NA,4)
return(out)
})

})

newfreeze=lapply(splitc,function(x) (getit.rfunc(x,realframesper,dat)))

newyo=do.call("rbind",newfreeze)
newyo=data.frame(newyo,stringsAsFactors=F)

newyo$flengths=as.vector (unlist (lapply(splitc,FUN=nrow)))

newyo[,1]=as.numeric(newyo[,1])
newyo[,2]=as.numeric(newyo[,2])

newyo=newyo[which (complete.cases(newyo)),]

fdiffs=sapply (seq(2,nrow(newyo)),function(x) (newyo[x,1]-newyo[x-1,2]))


newyo$fdiffs=c(NA,fdiffs)
newyo$smooth=NA
newyo$smooth[which(newyo$fdiffs<=0.1*realframesper)]="yes"

if(length(which(newyo$smooth=="yes"))>0) ({
smidx=which(newyo$smooth=="yes")
goods=which (is.na (newyo$smooth))
repidx=unlist (lapply (as.list (smidx),function(x) (goods [max (which (goods<x))])))

newyo[repidx,2]=newyo[smidx,2]
newyo[repidx,4]=newyo[smidx,4]

newyo=newyo[-which(newyo$smooth=="yes"),]
})
fdat_fnum=newyo[,c(1,2)]
fdat_vidtime=newyo[,c(3,4)]


fvidtime=cbind (unlist (lapply (strsplit(as.character (fdat_vidtime[,c(1)])," "),function(x) (x[2]))),unlist (lapply (strsplit(as.character (fdat_vidtime[,c(2)])," "),function(x) (x[2]))))
fdat=data.frame (cbind(dat$time[fdat_fnum[,1]],dat$time[fdat_fnum[,2]]),stringsAsFactors=F)

if (length (which (is.na(fdat[,1])))>0) ({
fdat=fdat[-which (is.na(fdat[,1])),]
fdat_fnum=fdat_fnum[-which (is.na(fdat_fnum[,1])),]
fvidtime=fvidtime[-which (is.na(fvidtime[,1])),]
})


fdat$f_start=fdat_fnum[,1]
fdat$f_stop=fdat_fnum[,2]

applythresh=fdat[,2]-fdat[,1]
subthresh_freeze=NA

if (length(which(applythresh<=freeze_rule))>0) ({
subthresh_freeze=fdat[which(applythresh<=freeze_rule),]
fdat=fdat[-which(applythresh<=freeze_rule),]
fdat_fnum=fdat_fnum[-which(applythresh<=freeze_rule),]
fvidtime=fvidtime[-which(applythresh<=freeze_rule),]
})

dat$freeze="NA"
forthis=apply(fdat_fnum,1,function(x) (seq(x[1],x[2])))
forthisu=as.vector (unlist(forthis))
dat$freeze [forthisu]="yes"

dat$fint=NA

dat$fint [forthisu]=unlist (lapply(as.list(seq(1,nrow(fdat_fnum))),function(x) (rep(x,length(forthis[[x]])))))
fmatch=match (unique (dat$fint [which (dat$freeze=="yes")]),dat$fint)

fdat$tint=NA
fdat$tint=dat$tint[fmatch]

colnames(fdat)=c("start","stop","start_f","stop_f","tint")#,"cue","precue")

fdat$cue="NA"
fdat$precue="no"

cue_freeze=which(dat$freeze=="yes") [which (which(dat$freeze=="yes") %in% which(dat$event=="cue"))]

iidx=as.numeric (unique (dat$fint[cue_freeze]))
fdat$cue[iidx]="yes" 

precueidx=cuedf$f_start-1
fpre=which (dat$freeze[precueidx]=="yes")

if (length(fpre)>0) ({
fdat$precue [as.numeric (dat$fint[precueidx][fpre])]="yes"


fdat$tint [which (fdat$precue=="yes")]=as.numeric (fdat$tint [which (fdat$precue=="yes")])+1

if (length(which(fdat$tint==0))>0) ({ 
fvidtime=fvidtime[-which(fdat$tint==0),]
fdat=fdat[-which(fdat$tint==0),]})
fdat$cuetime=cuedf$start [as.numeric (fdat$tint)]
})

fdat$tint=as.numeric(fdat$tint)
fdat$zone=dat$zone [match (fdat$start,dat$time)]


if ( length(which(fdat$tint==0))>0) ({
fdat$cuezone=c(rep (NA,length(which(fdat$tint==0))),cuedf$zone [fdat$tint])
fdat$tint[which(fdat$tint==0)]=NA
}) else ({
fdat$cuezone=cuedf$zone [fdat$tint]
})

fdat$ftype="pre_e"

fdat$ftype [which (fdat$zone!=fdat$cuezone)]="post_e"

fdat$dur=fdat$stop-fdat$start
fdat$precue_dur=0

if (length (which(fdat$precue=="yes"))>0) ({
pcidx=which(fdat$precue=="yes")
fdat$precue_dur[pcidx]=round (fdat$cuetime[pcidx]-fdat$start[pcidx],digits=3)
})


if ( length(which(is.na (fdat$tint)))>0) ({
fdat$tint [which(is.na (fdat$tint))]=NA
})


fdat$cuedur=NA
fdat$cue_stop=cuedf$stop[fdat$tint]

if(length (which (fdat$cue_stop>fdat$stop))>0) ({
fdat$cue_stop [which (fdat$cue_stop>fdat$stop)]=fdat$stop [which (fdat$cue_stop>fdat$stop)]
})

newcuedurs=fdat$cue_stop- fdat$start-fdat$precue_dur

if (length(which(fdat$cue=="yes"))>0) ({
fdat$cuedur[which(fdat$cue=="yes")]= round (newcuedurs [which(fdat$cue=="yes")],digits=3)
})

fdat$vidstart=fvidtime[,1]
fdat$vidstop=fvidtime[,2]

if (length(which (fdat$cuedur<0))>0) ({
fixidx=which (fdat$cuedur<0)
fdat$cuetime [fixidx]=cuedf$start [fdat$tint [fixidx]-1]
fdat$precue_dur [fixidx]=fdat$cuetime [fixidx]-fdat$start [fixidx]
fdat$cuedur [fixidx]=fdat$dur [fixidx]-fdat$precue_dur [fixidx]
})


##########################################################this is where i stopped vetting
freeze_mov=fdat$stop-cuedf$emov[fdat$tint]
fdat$freeze_mov=round (freeze_mov,digits=3)

if (length (which(fdat$cue=="yes"))>0) ({
fdat$freeze_mov [which(fdat$cue!="yes")]="NA"
})

fdat$ttype=cuedf$shock [fdat$tint]

eidx=cbind (cuedf$f_emov,cuedf$esc_frame)

if (length (which (is.na(eidx[,1])==TRUE))>0) ({
eidx=eidx [-which (is.na(eidx[,1])),]
})


cuedf$real_path=NA

fdat$ttype [which(fdat$ttype=="no")]="esc"
fdat$ttype [which(fdat$ttype=="yes")]="shock"
fdat$freeze_mov[which(fdat$ftype=="post_e")]=NA

fdat$freeze_mov[which(fdat$ttype=="shock")]=NA

formax=cbind (cuedf$f_start,cuedf$esc_frame)
formax=formax[esc_idx,]

if (length(formax)!=0) ({

if (is.null (dim(formax))==F) ({
maxidx=apply(formax,1,function(x) (dat[x[1]:x[2],]))
}) else ({
maxidx=list()
maxidx[[1]]=dat[formax[1]:formax[2],]
})



cuedf$maxframe=NA
cuedf$maxframe[esc_idx]=unlist (lapply(maxidx,function(x) (x$fnum [which.max(x$dists)])))
cuedf$maxtime=round (dat$time[cuedf$maxframe],digits=3)
cuedf$lat_tomax=round (cuedf$maxtime-cuedf$start,digits=3)
}) else ({

cuedf$maxframe=NA
cuedf$maxframe=NA
cuedf$maxtime=NA
cuedf$lat_tomax=NA
})


###################################apply freezing criterion

itistart=cuedf$esc_start[1:nrow(cuedf)-1]
itiend=cuedf$start[2:nrow(cuedf)]-0.001
if(is.na (itistart[1])==TRUE) (itistart[1]=itiend[1]-30)

itibind=cbind(itistart,itiend)


if (length (which (is.na(itibind[,1])))>0) ({
nostart=which (is.na(itibind[,1]))
itibind [nostart,1]= cuedf$stop[nostart-1]
})

dat$iti="NA"


ititrack=apply(itibind,1,function(x) (which(dat$time>x[1]) [which(dat$time>x[1]) %in% which(dat$time<x[2])]))


dat$iti [unlist(ititrack)]="iti"
dat$iti [which(dat$tint==0)]="iti"


cuedf$t_end=NA
cuedf$t_end_f=NA
gettend=split(dat,dat$tint)
getend=lapply (gettend [match (cuedf$tint,as.numeric (names(gettend)))],function(x) (x[nrow(x),]))
## get max velocity during escape

cuedf$t_end=as.vector (unlist (lapply(getend,function(x) (x$time))))
cuedf$t_end_f=as.vector (unlist (lapply(getend,function(x) (x$fnum))))

formax=cbind (cuedf$f_start,cuedf$t_end_f)

dat$sub_thresh_freeze=NA
dat$sub_thresh_freeze [apply(formax,1,function(x) (dat$fnum[x[1]:x[2]] [which.max (dat$dists[x[1]:x[2]])]))]="freeze"

csvname=paste (inf,".csv",sep="")

write.csv(cuedf,file=paste (csvdir,paste("cues",csvname,sep="_"),sep="/"),row.names=F)
write.csv(fdat,file=paste (csvdir,paste("freeze",csvname,sep="_"),sep="/"),row.names=F)

cueout[[i]]=cuedf
fout[[i]]=fdat

})

save (fout,file=paste (paste (root,"rout",sep=""),"fout.rdat",sep="/"))
save (cueout,file=paste (paste (root,"rout",sep=""),"cueout.rdat",sep="/"))









