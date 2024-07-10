



cleantrack.rfunc=function(s,a) ({

if (a=="nose") ({
colmatchx=4
colmatchy=5
})


if (a=="center") ({
colmatchx=2
colmatchy=3
})


tscan=s
#fr=a
#tscan=tscan[1:nrow(tscan),c(c(1,3,4,13),seq(15,19))]
#tscan=data.frame (rbind (tscan),stringsAsFactors=F)
#colnames(tscan)=c("ttime","posx","posy","zone","in","rec","rewA","startA","centerA")
#tscan$fnum=tscan$ttime
#tscan$ttime=round (tscan$ttime* (1/fr),digits=3)
missingdat=which(tscan[,colmatchx]=="-1")


if (length(missingdat)!=0)  ({
if (missingdat[1]==1 & length(missingdat)!=0 ) ({
getrid=rle(tscan[,colmatchx])
tscan=tscan[-getrid$lengths[1],]
row.names(tscan)=seq(1,nrow(tscan))
missingdat=missingdat[-seq(1,getrid$lengths[1])]
})
if (length(missingdat)!=0 ) ({
if (missingdat[length(missingdat)]==nrow(tscan) ) ({
getrid=rle(tscan[,colmatchx])
getridL=getrid$lengths[length(getrid$lengths)]
endidx=seq ((nrow(tscan)-(getridL-1)),nrow(tscan))
tscan=tscan[-endidx,]
row.names(tscan)=seq(1,nrow(tscan))
missingdat=missingdat[-match(endidx,missingdat)]
})
})
tscan[,colmatchx][missingdat]=0
tscan[,colmatchy][missingdat]=0
tscan$mpts=NA
tscan[,colmatchx]=as.numeric(tscan[,colmatchx])
tscan[,colmatchy]=as.numeric(tscan[,colmatchy])
missmax=0
if (length(missingdat)>0) ({
tscan$mpts[missingdat]="bad"
rlepts=rle(tscan$mpts)
forep=seq(1,length(rlepts$lengths))
tscan$mint=rep (forep,rlepts$lengths)
last=tscan$mint [which(tscan$mpts=="bad")]-1
lastint=tscan[which (tscan$mint%in% last),]
splitlast=split(lastint,lastint$mint)
lastk=do.call ("rbind",lapply(splitlast,function(x) (x[nrow(x) ,c(colmatchx,colmatchy)])))
nexts=tscan$mint [which(tscan$mpts=="bad")]+1
nextint=tscan[which (tscan$mint%in% nexts),]
splitnext=split(nextint,nextint$mint)
nextk=do.call ("rbind",lapply(splitnext,function(x) (x[1 ,c(colmatchx,colmatchy)])))
bpts=tscan[which(tscan$mpts=="bad"),]
bsplit=split(bpts,bpts$mint)
missmax=max (as.vector (unlist (lapply(bsplit,FUN=nrow))))
int.rfunc=function(x,y,z,q) ({
if (x==0) ({
out=y
})
if (q>1) ({
out=round (seq ((y+x),(z-x),length.out=q))
})
if (q==1) ({
out=round (y+x)
})
return(out)
})

colnames(nextk)=c("posx","posy")
colnames(lastk)=c("posx","posy")

newx=list()
newy=list()
for (p in (1:length(bsplit))) ({
totalpts=nrow(bsplit[[p]])
direx=(nextk$posx[p]-lastk$posx[p])/totalpts
direy=(nextk$posy[p]-lastk$posy[p])/totalpts
if (totalpts==1) ({
direy=direy/2
})
if (totalpts==1) ({
direx=direx/2
})
newx[[p]]=int.rfunc(direx,lastk$posx[p],nextk$posx[p],totalpts)
newy[[p]]=int.rfunc(direy,lastk$posy[p],nextk$posy[p],totalpts)
})
tscan[,colmatchx][missingdat]=unlist(newx)
tscan[,colmatchy][missingdat]=unlist(newy)
})
#tscan$posx=tscan$posx*0.1
#tscan$posy=tscan$posy*0.1
#dists=unlist (lapply (as.list (seq(1,nrow(tscan)-1)),function(x) (sqrt ( ((tscan$posx[x+1]-tscan$posx[x])^2)+((tscan$posy[x+1]-tscan$posy[x])^2)))))
#tscan$dists=c(NA,dists)
#tscan$dists=round (tscan$dists,digits=3)

})
tdat=tscan
return(tdat)



})
#})


save(cleantrack.rfunc,file=paste(funcdir,"cleantrack.rfunc",sep=""))
