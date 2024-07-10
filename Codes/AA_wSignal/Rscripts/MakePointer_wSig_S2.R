root="C:/Users/Jun/Documents/R/AA_wSignal/"  #this is the root folder

###########################you should only need to change the root directory to poinit towards your working directory
#following sub folders rawtrack, rawtrack_ttl, allsigs, sig_ttls 

tdir=paste(root,"rawtrack",sep="")   #put topscan files in rawtrack
tfiles=list.files(tdir)

tdirttl=paste(root,"rawtrack_ttl",sep="")  #put anymaze TTL files in  rawtrack_ttl
ttlfl=list.files(tdirttl)

formatch=gsub(".csv","",ttlfl)

sigdir=paste(root,"rawsig",sep="") #put all of your GC and AF files in allsigs
sigmatch= gsub (".csv","",list.files(sigdir))


sigpulse=paste(root,"rawsig_ttl",sep="")  #put all of your Synapse cue and shock files in sig_ttls
pfiles=list.files(sigpulse)
pfiles=gsub("CUE","cue",pfiles)
pfiles=gsub("SHOCK","shock",pfiles)
pulsematch=gsub(".csv","",pfiles)


allids=list()
for (i in (1:length(tfiles))) ({

tfname=paste(tdir,tfiles[i],sep="/")

getinf=read.table (tfname,nrows=22, sep = '\t', header = F,,strip.white=T,blank.lines.skip=T,fill=F,stringsAsFactors=F)
fps=as.numeric (strsplit (getinf[7,],":")[[1]][2])
starttime=strsplit (getinf[10,],":")[[1]][2]
stime=gsub("[[:space:]]", "",starttime)
stime=as.numeric (gsub("\\(s)", "",stime))

splitinf=strsplit (getinf[1,],"\\\\")[[1]]
vidname=splitinf[length(splitinf)]
inf=gsub (".BMP","",splitinf[length(splitinf)])
inf=gsub(".wmv","",inf)
inf=gsub(".avi","",inf)

dateraw=strsplit (gsub(" ","_",inf),"_")[[1]][1]
edate=as.Date(dateraw,format="%m%d%y")

am=formatch [match(inf,formatch)]
afmatch=sigmatch [match (paste(inf,"_AF",sep=""),sigmatch)]
gcmatch=sigmatch [match (paste(inf,"_GC",sep=""),sigmatch)]

cue_pulse=pulsematch [match (paste(inf,"_cue",sep=""),pulsematch)]
shock_pulse=pulsematch [match (paste(inf,"_shock",sep=""),pulsematch)]

out=c(as.character (edate),vidname,rep(NA,6),tfiles[i],am,afmatch,gcmatch,cue_pulse,shock_pulse)

allids[[i]]=out
})

bigout=data.frame (do.call("rbind",allids),stringsAsFactors=F)
colnames(bigout)=c("expt_date","Video_Name","Behavior_Type","Behavior_Subtype","Experiment_Group","Experiment_Group_Code","Experiment_Subgroup","Light_Status",
"tfile","am_file","af_file","gc_file","cue_pulse","shock_pulse")


bigout[order (as.Date (bigout$expt_date)),]

fnameout=paste(root,"pointer_raw",sep="")

write.csv(bigout,file=paste(fnameout,"finalpointer.csv",sep="/"),row.names=F)











