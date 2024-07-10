root="C:/Users/Jun/Documents/R/AA_wSignal/"

pre=2  #for peri-event histogram in seconds
post=2 #for peri-event histogram  in seconds
preheat=10 #for heatmap in seconds
postheat=30 #for heatmap in seconds

sps=0.1258291 #signals/second

tdir=paste(root,"rtrack",sep="")
tfiles=list.files(tdir)

outpath=paste(root,"output_dump/",sep="")

pointer=read.csv (paste(root,"pointer/finalpointer.csv",sep=""),stringsAsFactors=F)

load (paste (paste (root,"rout",sep=""),"cueout.rdat",sep="/"))
load (paste (paste (root,"rout",sep=""),"fout.rdat",sep="/"))
load (paste (paste (root,"rout",sep=""),"sigout.rdat",sep="/"))


lidx=match (gsub(".rdat","",tfiles),gsub (".wmv","",pointer$Video_Name))

binidx=round (seq(1,16)*sps,digits=2)
bins=c(rev(-binidx),0, binidx)
bins=sapply(seq(1,length(bins)-1),function(s) (paste (bins[s],bins[s+1],sep="_")))


binidxheat1=round (seq(1,80)*sps,digits=2)
binidxheat2=round (seq(1,240)*sps,digits=2)
binsheat=c(rev(-binidxheat1),0, binidxheat2)
binsheat=sapply(seq(1,length(binsheat)-1),function(s) (paste (binsheat[s],binsheat[s+1],sep="_")))


bigb=list()
bigsig=list()
forheatshock=list()
forheatcue=list()
allcue_onset_out=list()
esc_mov_onset_out=list()
at_esc_out=list()
at_maxvel_out=list()
cue_freeze_out=list()
post_esc_freeze_out=list()
fast_esc_cue_freeze_out=list()
slow_esc_cue_freeze_out=list()
fast_esc_freeze_out=list()
slow_esc_freeze_out=list()
allfreeze_out=list()


for (i in (1:length(cueout))) ({
finf=pointer[lidx [i],]

cuedf=cueout[[i]]
fdat=fout[[i]]
sig=sigout[[i]]


inf=gsub(".rdat","",finf$Video_Name)
splitinf=strsplit(inf," ")[[1]]
didx=splitinf[1]
id=paste(splitinf[length(splitinf)],sep="_")
outname=paste("ZS",didx,id,finf$stype,sep="_")

bindstype=rep ("shock",nrow(cuedf))
bindstype [which(cuedf$shock=="no")]="esc"
bindout=data.frame (cbind (inf,finf$stype,cuedf$tint,bindstype),stringsAsFactors=F)
colnames(bindout)=c("id","etype","tint","ttype")

sps=mean(diff(sig$time))
spshist=round (1/sps)*pre
spshist_post=round (1/sps)*post
spsheat_pre=round (1/sps)*preheat
spsheat_post=round (1/sps)*postheat

getsigs.rfunc=function(x,y,z,q) ({
x=x
y=y

getzero=sapply(x,function(s) (which(y$vid_fnum>=s)[1]))
getspikes=cbind (getzero-spshist,getzero+spshist)    

if (q=="heat") ({
getspikes=cbind (getzero-spsheat_pre,getzero+spsheat_post) 
})
                                                                                                                                                

tidx=seq(1,nrow(getspikes))

if(length (which (is.na(getspikes[,1])))>0) ({
naidx=which (is.na(getspikes[,1]))
getspikes=getspikes[-naidx,]
tobind=matrix (rep(NA,length(z)*length(naidx)),nrow=length(naidx))
colnames(tobind)=z
tobind=as.data.frame(tobind,stringsAsFactors=F)
tobind$rnum=naidx
tidx=tidx[-naidx]
})


sigsout=apply(getspikes,1,function(s) (y$z_base[s[1]:s[2]]))
sigsout=t(sigsout)

if (nrow(sigsout)>1) ({
sigsout=sigsout[,-ncol(sigsout)]
})


if (nrow(sigsout)==1) ({
sigsout=sigsout[,-length(sigsout)]
sigsout=as.data.frame(rbind (sigsout),stringsAsFactors=F)
})


colnames(sigsout)=z
sigsout=as.data.frame(sigsout,stringsAsFactors=F)
sigsout$rnum=tidx

if (nrow(sigsout)!=length(x)) ({
newsigs=as.data.frame(rbind(tobind,sigsout),stringsAsFactors=F)
newsigs=newsigs[order(newsigs$rnum),]

sigsout=newsigs
})

sigsout=sigsout[,-ncol(sigsout)]

return(sigsout)
})

allcue_onset=getsigs.rfunc(cuedf$f_start,sig,bins,"noheat")
write.csv(allcue_onset,paste (outpath,paste (paste("allcue_onset",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)


if (length (cuedf$f_start[which(cuedf$shock=="no")])>0) ({
allcue_onset_heat_noshock=getsigs.rfunc(cuedf$f_start[which(cuedf$shock=="no")],sig,binsheat,"heat")
}) else ({
allcue_onset_heat_noshock=rep (NA, length (binsheat))
})

if (length (cuedf$f_start[which(cuedf$shock=="yes")])>0) ({
allcue_onset_heat_shock=getsigs.rfunc(cuedf$f_start[which(cuedf$shock=="yes")],sig,binsheat,"heat")
}) else ({
allcue_onset_heat_shock=rep (NA, length (binsheat))
})


esc_mov_onset=getsigs.rfunc(cuedf$f_emov,sig,bins,"noheat")
write.csv(esc_mov_onset,paste (outpath,paste (paste("esc_mov_onset",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)
at_esc=getsigs.rfunc(cuedf$esc_frame,sig,bins,"noheat")
write.csv(at_esc,paste (outpath,paste (paste("at_esc",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)
at_maxvel=getsigs.rfunc(cuedf$maxframe,sig,bins,"noheat")
write.csv(at_maxvel,paste (outpath,paste (paste("at_maxvel",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)


if(fdat$start[1]<0) ({
fdat=fdat[-1,]
})

cue_freeze=getsigs.rfunc(fdat$start_f[which(fdat$cue=="yes")],sig,bins,"noheat")
write.csv(cue_freeze,paste (outpath,paste (paste("cue_freeze",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)

post_esc_freeze_idx=fdat$start_f[which(fdat$ftype=="post_e") [which (which(fdat$ftype=="post_e") %in% which(fdat$ttype=="esc"))]]
idx_eventout=fdat$tint[which(fdat$ftype=="post_e") [which (which(fdat$ftype=="post_e") %in% which(fdat$ttype=="esc"))]]

if (length(post_esc_freeze_idx)>0) ({
post_esc_freeze=getsigs.rfunc(post_esc_freeze_idx,sig,bins,"noheat")
}) else ({
post_esc_freeze=rep(NA,length(bins))
})

write.csv(post_esc_freeze,paste (outpath,paste (paste("post_esc_freeze",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)



splitmids=median(cuedf$esc_lat,na.rm=T)

fast=cuedf$tint [which(cuedf$esc_lat<splitmids)]
slow=cuedf$tint [which(cuedf$esc_lat>splitmids)]

fastT=fdat[which (fdat$tint %in% cuedf$tint[fast]),]
slowT=fdat[which (fdat$tint %in% cuedf$tint[slow]),]

fastT_cue=fastT[which(fastT$cue=="yes"),] 
slowT_cue=slowT[which(slowT$cue=="yes"),]

fast_esc_freeze=getsigs.rfunc(fastT$start_f,sig,bins,"noheat")
write.csv(fast_esc_freeze,paste (outpath,paste (paste("fast_esc_cue_freeze",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)
slow_esc_freeze=getsigs.rfunc(slowT$start_f,sig,bins,"noheat")
write.csv(slow_esc_freeze,paste (outpath,paste (paste("slow_esc_cue_freeze",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)

fast_esc_cue_freeze=getsigs.rfunc(fastT_cue$start_f,sig,bins,"noheat")
write.csv(fast_esc_cue_freeze,paste (outpath,paste (paste("fast_esc_cue_freeze",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)
slow_esc_cue_freeze=getsigs.rfunc(slowT_cue$start_f,sig,bins,"noheat")
write.csv(slow_esc_cue_freeze,paste (outpath,paste (paste("slow_esc_cue_freeze",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)

shock_num=length(which(cuedf$shock=="yes"))
if (shock_num!=0) ({
cuedf$dur [which(cuedf$shock=="yes")]=15
})

esc_num=length(which(cuedf$shock=="no"))
esc_lat=round (mean (cuedf$esc_lat [which(cuedf$shock=="no")],na.rm=T),digits=2)
e_lat=round (mean (cuedf$esc_lat [which(cuedf$shock=="no")],na.rm=T),digits=2)
e_mov_lat=round (mean (cuedf$emov_cuelat [which(cuedf$shock=="no")],na.rm=T),digits=2)
lat_tomax=round (mean (cuedf$lat_tomax [which(cuedf$shock=="no")],na.rm=T),digits=2)

cue_freeze_dur_total=round (sum (fdat$dur[which(fdat$cue=="yes")],na.rm=T),digits=2)
cue_freeze_num=nrow (fdat[which(fdat$cue=="yes"),])
cue_freeze_dur_avg=round (mean(fdat$dur[which(fdat$cue=="yes")],na.rm=T),digits=2)
total_freeze_time=round (sum (fdat$dur,na.rm=T),digits=2)
total_cue_time=sum(cuedf$dur)

write.csv(cuedf,paste (outpath,paste (paste("cues",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)
write.csv(fdat,paste (outpath,paste (paste("freeze",outname,sep="_"),".csv",sep=""),sep=""),row.names=F)

bout=data.frame (cbind(shock_num,esc_num,e_lat,e_mov_lat,lat_tomax,cue_freeze_dur_total,cue_freeze_num,cue_freeze_dur_avg,total_freeze_time,total_cue_time),stringsAsFactors=F)
bigb[[i]]=data.frame (cbind(finf$Video_Name,finf$stype,bout),stringsAsFactors=F)


getmean.rfunc=function(s) ({
#at_esc some NA
#post_esc_freeze all NA
#
if (is.null (ncol(s))==T) ({
return(s)
}) else ({
return (apply (s,2,function(s) (mean(s,na.rm=T))))
})
})

avgsig=rbind (getmean.rfunc(allcue_onset),getmean.rfunc(esc_mov_onset),getmean.rfunc(at_esc),getmean.rfunc(at_maxvel),getmean.rfunc(cue_freeze),
getmean.rfunc(post_esc_freeze),getmean.rfunc(fast_esc_cue_freeze),getmean.rfunc(slow_esc_cue_freeze))
avgsig=data.frame(avgsig,stringsAsFactors=F)
colnames(avgsig)=bins
avgsig$event=c("allcue_onset","esc_mov_onset","at_esc","at_maxvel","cue_freeze","post_esc_freeze","fast_esc_cue_freeze","slow_esc_cue_freeze")
avgsig=data.frame (cbind(id=finf$Video_Name,etype=finf$stype,avgsig),stringsAsFactors=F)
colnames(avgsig)=c("id","etype",bins,"event")

bigsig[[i]]=avgsig

forheatcue[[i]]=allcue_onset_heat_noshock
forheatshock[[i]]=allcue_onset_heat_shock

allcue_onset_out[[i]]=cbind (bindout,event="cue_onset",allcue_onset)
esc_mov_onset_out[[i]]=cbind (bindout,event="esc_mov_onset",esc_mov_onset)
at_esc_out[[i]]=cbind (bindout,event="at_esc",at_esc)
at_maxvel_out[[i]]=cbind (bindout,event="at_maxvel",at_maxvel)
cue_freeze_out[[i]]=cbind (bindout [fdat$tint [which(fdat$cue=="yes")],],event="cue_freeze",cue_freeze)

if (length(idx_eventout)>0) ({
post_esc_freeze_out[[i]]=cbind (bindout [idx_eventout,],event="post_esc_freeze",post_esc_freeze)
})

fast_esc_cue_freeze_out[[i]]=cbind (bindout [fastT_cue$tint,],event="fast_esc_cue_freeze",fast_esc_cue_freeze)
slow_esc_cue_freeze_out[[i]]=cbind (bindout [slowT_cue$tint,],event="slow_esc_cue_freeze",slow_esc_cue_freeze)
fast_esc_freeze_out[[i]]=cbind (bindout [fastT$tint,],event="fast_esc_freeze",fast_esc_freeze)
slow_esc_freeze_out[[i]]=cbind (bindout [slowT$tint,],event="slow_esc_freeze",slow_esc_freeze)

onlycuefreeze=fdat
freeze_split=split(onlycuefreeze,onlycuefreeze$tint)
totalcuedur=unlist (lapply(freeze_split,function(x) (sum(x$cuedur,na.rm=T))))
cue_freeze_dur=round(totalcuedur,digits=2)
forcuetint=as.vector(unlist (lapply(freeze_split,function(x) (x$tint[1]))))
if (max (forcuetint)>nrow(cuedf)) ({
forcuetint=forcuetint[-length(forcuetint)]
cue_freeze_dur=cue_freeze_dur[-length(cue_freeze_dur)]
})


bindout$cue_freeze_dur=0.00
bindout$cue_freeze_dur [match(forcuetint,as.numeric (bindout$tint))]=cue_freeze_dur
bindout$cue_dist=cuedf$cue_dist

allfreeze_out[[i]]=bindout




})


allb=data.frame (do.call("rbind",bigb),stringsAsFactors=F)
colnames(allb)[c(1,2)]=c("id","etype")

allb=data.frame (do.call("rbind",split(allb,allb$etype)),stringsAsFactors=F)
row.names(allb)=seq(1,nrow(allb))

write.csv(allb,paste (outpath,paste (paste("summary","behavior",sep="_"),".csv",sep=""),sep=""),row.names=F)

sigsout=data.frame (do.call("rbind",bigsig),stringsAsFactors=F)
colnames(sigsout)=c("id","etype",bins,"event")
sigsout$id=unlist(as.vector(sigsout$id))
sigsout$etype=unlist(as.vector(sigsout$etype))


sigsout$spliti=paste(sigsout$event,sigsout$etype,sep="_")

splitsigs=split(sigsout,sigsout$spliti)
sigsout=data.frame (do.call("rbind",splitsigs),stringsAsFactors=F)
colnames(sigsout)=c("id","etype",bins,"event","spliti")

sigsout$idx=unlist (lapply (strsplit(sigsout$spliti,"_"),function(x) (x[length(x)])))

splitsigs=split(sigsout,sigsout$idx)
sigsout=data.frame (do.call("rbind",splitsigs),stringsAsFactors=F)
colnames(sigsout)=c("id","etype",bins,"event","spliti","idx")
row.names(sigsout)=seq(1,nrow(sigsout))

sigsout=sigsout[,-grep("idx",colnames(sigsout))]
colnames(sigsout)=c("id","etype",bins,"event","atype")

outhist=paste (paste("summary","ZS","perievent_raw",sep="_"),".csv",sep="")
write.csv(sigsout,paste (outpath,outhist,sep=""),row.names=F)

newsplit=split(sigsout,sigsout$atype)

getmeans.rfunc=function(x) ({
formean=seq (grep("etype",colnames(x))+1,grep("event",colnames(x))-1)
meanies=apply (x[,formean],2,function(s) (mean(s,na.rm=T)))
sds=apply (x[,formean],2,function(s) (sd(s,na.rm=T)))
se= apply (x[,formean],2,function(s) (sd(s,na.rm=T)/sqrt(length(s))))

out=data.frame (rbind(round (meanies,digits=3),round (sds,digits=3),round (se,digits=3)),stringsAsFactors=F)
out$measure=c("mean","sd","se")
out$event=x$event[1]
out$atype=x$atype[1]
out=data.frame (cbind(rep(x$id[1],nrow(out)),rep(x$etype[1],nrow(out)),out),stringsAsFactors=F)
out[,1]=unlist(as.vector(out[,1]))
out[,2]=unlist(as.vector(out[,2]))

colnames(out)=colnames(x)
out=out[,-1]
return(out)
})


sumsig=lapply(newsplit,function(x) (getmeans.rfunc(x)))
sumhist=data.frame (do.call("rbind",sumsig),stringsAsFactors=F)
row.names(sumhist)=seq(1,nrow(sumhist))

colnames(sumhist)=c("etype",bins,"measure","event","atype")

outhist=paste (paste("summary","ZS","perievent_mean",sep="_"),".csv",sep="")
write.csv(sumhist,paste (outpath,outhist,sep=""),row.names=F)


cleannames.rfunc=function(x) ({
colnames (x)=gsub ("X.","",colnames (x))
return(x)
})

allcue_onset_out=data.frame (do.call("rbind",allcue_onset_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(allcue_onset_out),paste (outpath,"summary_ZS_allcue_onset.csv",sep=""),row.names=F)

esc_mov_onset_out=data.frame (do.call("rbind",esc_mov_onset_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(esc_mov_onset_out),paste (outpath,"summary_ZS_esc_mov_onset.csv",sep=""),row.names=F)

at_esc_out=data.frame (do.call("rbind",at_esc_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(at_esc_out),paste (outpath,"summary_ZS_at_esc.csv",sep=""),row.names=F)

at_maxvel_out=data.frame (do.call("rbind",at_maxvel_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(at_maxvel_out),paste (outpath,"summary_ZS_at_maxvel.csv",sep=""),row.names=F)

cue_freeze_out=data.frame (do.call("rbind",cue_freeze_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(cue_freeze_out),paste (outpath,"summary_ZS_cue_freeze.csv",sep=""),row.names=F)

post_esc_freeze_out=data.frame (do.call("rbind",post_esc_freeze_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(post_esc_freeze_out),paste (outpath,"summary_ZS_post_esc_freeze.csv",sep=""),row.names=F)

fast_esc_cue_freeze_out=data.frame (do.call("rbind",fast_esc_cue_freeze_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(fast_esc_cue_freeze_out),paste (outpath,"summary_ZS_fast_esc_cue_freeze.csv",sep=""),row.names=F)

slow_esc_cue_freeze_out=data.frame (do.call("rbind",slow_esc_cue_freeze_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(slow_esc_cue_freeze_out),paste (outpath,"summary_ZS_slow_esc_cue_freeze.csv",sep=""),row.names=F)

fast_esc_freeze_out=data.frame (do.call("rbind",fast_esc_freeze_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(fast_esc_freeze_out),paste (outpath,"summary_ZS_fast_esc_freeze.csv",sep=""),row.names=F)

slow_esc_freeze_out=data.frame (do.call("rbind",slow_esc_freeze_out),stringsAsFactors=F)
write.csv(cleannames.rfunc(slow_esc_freeze_out),paste (outpath,"summary_ZS_slow_esc_freeze.csv",sep=""),row.names=F)



all_cue_freeze_out=data.frame (do.call("rbind",allfreeze_out),stringsAsFactors=F)
write.csv(all_cue_freeze_out,paste (outpath,"summary_all_cue_freeze_out.csv",sep=""),row.names=F)


outheat_1=paste (paste("summary","ZS","heat_matrix_noshock",sep="_"),".csv",sep="")
bigheat_noshock=data.frame (do.call("rbind",forheatcue),stringsAsFactors=F)
colnames (bigheat_noshock)= gsub ("X.","",colnames (bigheat_noshock))
write.csv(bigheat_noshock,paste (outpath,outheat_1,sep=""),row.names=F)

outheat_2=paste (paste("summary","ZS","heat_matrix_shock",sep="_"),".csv",sep="")
bigheat_shock=data.frame (do.call("rbind",forheatshock),stringsAsFactors=F)
colnames (bigheat_shock)= gsub ("X.","",colnames (bigheat_shock))
write.csv(bigheat_shock,paste (outpath,outheat_2,sep=""),row.names=F)








