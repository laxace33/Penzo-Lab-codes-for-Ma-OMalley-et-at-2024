old_root="C:/Users/Jun/Documents/R/AA_behavior/" #you will need to change this 
new_root="C:/Users/Jun/Documents/R/" #you will need to change this 

next_experimentName= "PL_pPVT_opto"  #you will need to change this 

###############################################################
newdir=paste(new_root,next_experimentName,sep="")

dir_list=list.dirs(old_root)

if (dir.exists(newdir)==F) ({

fnames=gsub(old_root,newdir,dir_list)

for (q in (1:length(fnames))) ({
dir.create(fnames[q],recursive=T)
})


tfunc=grep ( "track_func",dir_list)
mov_lsquare=fnames[tfunc]
file.copy(paste(dir_list[tfunc],list.files (dir_list[tfunc]),sep="/"), mov_lsquare)


script_idx=grep("Rscripts",dir_list)
mov_scripts=fnames[script_idx]
file.copy(paste(dir_list[script_idx],list.files (dir_list[script_idx]),sep="/"),mov_scripts)

})
