root="//eu.boehringer.com/users/bib/users1/duhoffma/Desktop/AA_pipeline/aa_test/"

load (paste (paste (root,"rout",sep=""),"cueout.rdat",sep="/"))
load (paste (paste (root,"rout",sep=""),"fout.rdat",sep="/"))

outpath=paste(root,"output_dump",sep="")

tdir=paste(root,"rtrack",sep="")
tfiles=list.files(tdir)

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

lidx=match (gsub(".rdat",".wmv",tfiles),pointer$Video_Name)

bigb=list()

for (i in (1:length(cueout))) ({
finf=pointer[lidx [i],]

inf=gsub(".rdat","",finf$Video_Name)
splitinf=strsplit(inf," ")[[1]]
didx=splitinf[1]
id=paste(splitinf[length(splitinf)],sep="_")
outname=paste(didx,id,finf$stype,sep="_")
outname=gsub(".wmv","",outname)


cuedf=cueout[[i]]
fdat=fout[[i]]

if (length(cuedf)!=1) ({

shock_num=length(which(cuedf$shock=="yes"))

if (shock_num!=0) ({
cuedf$dur [which(cuedf$shock=="yes")]=15
})

esc_num=length(which(cuedf$shock=="no"))
esc_lat=round (mean (cuedf$esc_lat [which(cuedf$shock=="no")],na.rm=T),digits=2)
e_lat=round (mean (cuedf$esc_lat [which(cuedf$shock=="no")],na.rm=T),digits=2)
e_mov_lat=round (mean (cuedf$emov_cuelat [which(cuedf$shock=="no")],na.rm=T),digits=2)
lat_tomax=round (mean (cuedf$lat_tomax [which(cuedf$shock=="no")],na.rm=T),digits=2)
total_cue_time=sum(cuedf$dur)

cue_freeze_dur_total=sum (round (fdat$cuedur[which(fdat$cue=="yes")],digits=2))

cue_freeze_num=nrow (fdat[which(fdat$cue=="yes"),])
cue_freeze_dur_avg=round (mean(fdat$cuedur[which(fdat$cue=="yes")],na.rm=T),digits=2)
total_freeze_time=round (sum (fdat$dur,na.rm=T),digits=2)
bout=data.frame (cbind(shock_num,esc_num,e_lat,e_mov_lat,lat_tomax,cue_freeze_dur_total,cue_freeze_num,cue_freeze_dur_avg,total_freeze_time,total_cue_time),stringsAsFactors=F)
}) else ({
bout=data.frame (rbind (rep(NA,10)),stringsAsFactors=F)
colnames(bout)=c ("shock_num","esc_num","e_lat","e_mov_lat","lat_tomax","cue_freeze_dur_total","cue_freeze_num","cue_freeze_dur_avg","total_freeze_time","total_cue_time")
})


write.csv(cuedf,paste (outpath,paste (paste("cues",outname,sep="_"),".csv",sep=""),sep="/"),row.names=F)
write.csv(fdat,paste (outpath,paste (paste("freeze",outname,sep="_"),".csv",sep=""),sep="/"),row.names=F)

bigb[[i]]=data.frame (cbind(finf$Video_Name,finf$stype,bout),stringsAsFactors=F)

})


allb=data.frame (do.call("rbind",bigb),stringsAsFactors=F)
colnames(allb)[c(1,2)]=c("id","etype")

allb$e_lat [which (allb$e_lat=="NaN")]=NA
allb$e_mov_lat [which (allb$e_mov_lat=="NaN")]=NA
allb$lat_tomax [which (allb$lat_tomax=="NaN")]=NA
allb$cue_freeze_dur_avg [which (allb$cue_freeze_dur_avg=="NaN")]=NA


newsplit=split(allb,allb$etype)

allb=data.frame (do.call("rbind",newsplit),stringsAsFactors=F)
row.names(allb)=seq(1,nrow(allb))

allb$id=gsub(".wmv","",allb$id)
write.csv(allb,paste (outpath,paste (paste("summary","behavior",sep="_"),".csv",sep=""),sep="/"),row.names=F)
















