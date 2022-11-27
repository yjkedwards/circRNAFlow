#!/usr/bin/env Rscript


getPValsInOrderOfAdj=function(padj_sorted,padj_unsorted,pvals_unsorted) {
	new_p_vals=c()
	padj_sorted_uniq=unique(padj_sorted)
	for(i in 1:length(padj_sorted_uniq)) {
		the_adj=padj_sorted_uniq[i]
		index_in_adj_unsorted=which(padj_unsorted==the_adj)
		the_pvals=pvals_unsorted[index_in_adj_unsorted]
		new_p_vals=c(new_p_vals,the_pvals)
		}
	if(length(new_p_vals)>length(padj_unsorted)) {
		#too many?  then chop extras and mark last with -99
		new_p_vals=new_p_vals[1:(length(new_p_vals)-1)]
		new_p_vals=c(new_p_vals,-99)
		}
	else if(length(new_p_vals)<length(padj_unsorted)) {
		#too few?  then add -99 to make up for the difference
		while(length(new_p_vals)<length(padj_unsorted)) {
			new_p_vals=c(new_p_vals,-99)
			}				
		}
	return(new_p_vals)
	}


getTestResultsAsDF=function(test) {
	#summary
	summ_df=test$summary_table
	#get pvals
	pvals_proper_order=getPValsInOrderOfAdj(summ_df[,"sig_p"],test$p.adj,test$p.val)
	pvals_proper_order_df=data.frame(p.val=pvals_proper_order)	
	#get sigd
	sigd_item=test$sig.dat
	sigd_df=sigd_item[,4:(dim(test$sig.dat)[2])]
	#########################################3
	#merge summ and sigd
	summ_sigd=cbind(summ_df,sigd_df)
	#merge that with the pvals
	summ_sigd_pv_df=cbind(summ_sigd,pvals_proper_order_df)
	return(summ_sigd_pv_df)
	}



args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
count_file=""
linear_file=""
coord_file=""
inc_grp_1_idx=0
if (length(args)!=6) {
  stop("Three arguments expected <COUNT> <LINEAR> <COORD> <INC_GRP1_IDX> <GRP1> <GRP2>", call.=FALSE)
} else  {
count_file=args[1]
linear_file=args[2]
coord_file=args[3]
inc_grp_1_idx=as.numeric(args[4])
grp_1=args[5]
grp_2=args[6]
}

print(paste("Using count file ",count_file,sep=""))
print(paste("Using linear file ",linear_file,sep=""))
print(paste("Using coordinate file ",coord_file,sep=""))
print(paste("Using  inc_grp_1_idx ",inc_grp_1_idx,sep=""))

#############################################
#read the files
CircRNACount <- read.delim(count_file,header=T)
LinearCount <- read.delim(linear_file,header=T)
CircCoordinates <- read.delim(coord_file,header=T)


###############################################
#perform group assignment
group_indicators=c()
group_indicator_names=c()
for(c_idx in 4:length(colnames(CircRNACount))) {
	if(c_idx<=inc_grp_1_idx) {
		#for indices at and below the given cut-off, put them in group 1
		group_indicators=c(group_indicators,1)
		group_indicator_names=c(group_indicator_names,grp_1)
		} else {
		#otherwise, put them in group 2
		group_indicators=c(group_indicators,2)
		group_indicator_names=c(group_indicator_names,grp_2)
		}
	}
if(length(group_indicators)<=0) {
	stop(paste("Unable to define group indicators?!  Invalid choice possibly : ",inc_grp_1_idx,sep=""))
	}
print("Using the following group assignments for indices 4 and higher :")
print(group_indicators)
Nreplicates=(inc_grp_1_idx+1)-4
print(paste("Nreplicates computed as ",Nreplicates,sep=""))
tot_num_samps=length(colnames(CircRNACount))-3
#verify that that each group has the same number of replicates
#if((Nreplicates*2)==tot_num_samps) {
#	print(paste("The two groups share the same number of replicates : ",Nreplicates))
#	} else {
#	#bad
#	msg=paste("Nreplicates computed as ",Nreplicates,".  tot_num_samps computed as ",tot_num_samps,".  It is expected that 2*Nreplicates=tot_num_samps so that each replicate group has the same size!",sep="")
#	stop(msg)
#	}




###############################################
#perform filtering
library(CircTest)
library(aod)
library("ggplot2")


#perform filtering with default parameters
CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount, Nreplicates = Nreplicates)
CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]
LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]
num_rows_after_filter=dim(CircRNACount_filtered)[1]
circTestWillHaveData=1
if(num_rows_after_filter==0) {
	print("No data to test for circtest plain!!  It all got filtered out!")
	circTestWillHaveData=0
}
################
# these defaults copied from the circtest filter script/repo  https://github.com/dieterich-lab/CircTest/blob/master/R/Circ.filter.R
filter.sample=4
filter.count=5
percentage=1
alpha=0.05
out_file_root=paste("filtercount",filter.count,"filtersample",filter.sample,"percentage",percentage,"alpha",alpha,sep="_")
outfile=paste(out_file_root,".tsv",sep="")
if(circTestWillHaveData==0) {
	system(paste("touch ",outfile,sep=""))
} else {
	# Run the test
	test=Circ.test(CircRNACount_filtered,LinearCount_filtered,CircCoordinates_filtered,group=group_indicators,alpha=alpha,circle_description = c(1:3))
	test_whole_df=getTestResultsAsDF(test)

	print(paste("Writing output to ",outfile))
	write.table(test_whole_df,outfile,sep="\t",row.names = FALSE)
	system(paste("touch ",outfile,sep=""))
	}

