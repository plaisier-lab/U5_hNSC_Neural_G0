runSsGSEAwithPermutation<-function(profile_data_file, number_perms){

    # NEW ADD 1/3
    profile_data_path = dirname(file.path(profile_data_file))
    profile_data_name = basename(profile_data_file)
    MOD_file = paste(system.file("data", package = "ssgsea.GBM.classification"),"/",sep='')


    #please keep the first and second column as the same with NAME and Description
    data<-read.table(paste(profile_data_path,profile_data_name,sep='/'),header=T,row.names=1,skip=2); # read data from a gct file and ignore first 2 lines and 1 column
    size<-dim(data[-1])            #size of data
    num_rows<-size[1]          #number of columns
    num_cols<-size[2]          # number of rows
    num_perm<-number_perms     # number of permutations
    colnames<-colnames(data)   #g et the column names

    print("Data loading finished!")
     
    line1<-'#1.2'
    line2<-paste(num_rows,num_perm,sep='\t')
    line3<-paste('X',1:num_perm,sep='',collapse='\t')
    line3<-paste('NAME','Description',line3, sep='\t')
    write.table(line1, paste(profile_data_path,'/random_profile_',profile_data_name,'.gct',sep=""),quote=F,col.name=F,row.name=F) #first line for the gct file
    write.table(line2, paste(profile_data_path,'/random_profile_',profile_data_name,'.gct',sep=""),quote=F,col.name=F,row.name=F,append=T)  # 2nd line for the gct file
    write.table(line3, paste(profile_data_path,'/random_profile_',profile_data_name,'.gct',sep=""),quote=F,col.name=F,row.name=F,append=T)  # 2nd line for the gct file


    random_profile<-data[1]    
    for (i in 1:num_perm){
        a<-data[-1][cbind(seq(1:num_rows),sample(1:num_cols,num_rows,replace=T))]  # the first column is the decription column, ignore it.
        random_profile<-cbind(random_profile,a)
        if(i%%100==0){
           print(i)
        }
    }
    write.table(random_profile, paste(profile_data_path,'/random_profile_',profile_data_name,'.gct',sep="") ,sep='\t',quote=F,col.name=F, append=T)


    print("Random profiles was genereated!")

    # source(msig_library_file)
    selected.models<-c("Proneural","Classical","Mesenchymal")
    OPAM.apply.model.2(  
        input.ds           = paste(profile_data_path,'/random_profile_',profile_data_name,'.gct',sep=''),
        models.dir         = MOD_file,
        models             = selected.models,
        raw.score.outfile  = paste(profile_data_path,'/random_raw.score_',profile_data_name,'.gct',sep=""),
        norm.score.outfile = "", model.score.outfile= "", prob.outfile= "",graphics.off= T)

    print("SsGSEA was performed on random profiles!");

    OPAM.apply.model.2(  
        input.ds           = paste(profile_data_path, profile_data_name,sep='/'),
        models.dir         = MOD_file,
        models             = selected.models,
        raw.score.outfile  = paste(profile_data_path,'/raw.score_',profile_data_name,'.gct',sep=""),
        norm.score.outfile = "", model.score.outfile= "", prob.outfile= "",graphics.off= T)

    print("SsGSEA was performed on the original profiles!");

    random_result<-read.table(paste(profile_data_path,'/random_raw.score_',profile_data_name,'.gct',sep=""),header=T,row.names=1,skip=2,sep="\t");
    random_result<-random_result[-1];
    random_result<-t(random_result);
    original_result<-read.table(paste(profile_data_path,'/raw.score_',profile_data_name,'.gct',sep=""),header=T,row.names=1,skip=2,sep="\t");
    original_result<-original_result[-1];
    original_result<-t(original_result)

    
    p_result<-original_result;
    for (i in 1:dim(original_result)[1]){
         p_result[i,]<-colSums(sweep(random_result,2,original_result[i,])>=0);
    }


    # NEW ADD 2/3
    file.remove(paste(profile_data_path,'/random_profile_',profile_data_name,'.gct',sep=""));
    file.remove(paste(profile_data_path,'/random_raw.score_',profile_data_name,'.gct',sep=""));
    file.remove(paste(profile_data_path,'/raw.score_',profile_data_name,'.gct',sep=""));

    
    # NEW ADD 3/3
    colnames(p_result)=paste(colnames(p_result),"pval",sep="_");
    p_result=t(apply(p_result,1,function(x){(x+1)/(number_perms+1)}))


    write.table(cbind(original_result,p_result), paste(profile_data_path,'/p_result_',profile_data_name,'.txt',sep="") ,sep='\t',quote=F)

    print("P_values for each subtype was calculated!")
    
    print ("Finished!")
    return(cbind(original_result,p_result))

}
