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
###############################################
# some documentation from the README https://github.com/dieterich-lab/CircTest
#Nreplicates specifies the number of replicates in each condition.
#filter.sample specifies the number of samples the circle has to have enough circular reads in to be considered.
#filter.count specifies the circular read count threshold.
#percentage specifies the minimum circle to host-gene ratio.
#circle_description tells the function which columns are NOT filled with read counts but the circle's annotation.
###############################################
# some documentation from the source from the following commit (most recent as of today) : https://github.com/dieterich-lab/CircTest/blob/a22301476e6f41bc1ecec3caae98971d643eccdd/R/Circ.filter.R
#' @param circ CircRNACount file. A file of circRNA read count table. First three columns are circRNA coordinates, and followed by columns for circRNA read counts, each sample per column.
#' @param linear LinearCount file. A file of circRNA host gene expression count table. Same configuration as CircRNACount file.
#' @param Nreplicates Number of replicates in your data. Expect each group have the same number of replicates.
#' @param filter.count The minimum read count used for filtering.
#' @param filter.sample The minimum number of samples need to have above filter.count number of circRNA supporting reads.
#' @param percentage The minimum percentage of circRNAs account for the total transcripts in at least one group.
#' @param circle_description Column indices which do not carry circle/linear read counts.
#' @export Circ.filter
############################################### 
# DEFAULT PARAMETERS
###Circ.filter <- function(circ=circ,linear=linear,Nreplicates=3,filter.sample=4,filter.count=5,percentage=1, circle_description=c(1:3)

#parameter space for filtering
percentage_bins=c(0.01,0.05,0.1,0.25,0.5,0.75,1.0)
filter.count_bins=c(1,5,10)
filter.sample_bins=1:(length(group_indicators))

print("percentage bins :")
print(percentage_bins)
print("filter.count_bins")
print(filter.count_bins)
print("filter.sample_bins")
if(length(filter.sample_bins)>=10) {
        new_filter.sample_bins=c()
	partial_bins=c()
	for(t in 1:5) {
		the_frac=t/5.0
		partial_bins=c(partial_bins,as.integer(the_frac*length(filter.sample_bins)))
		}
        for(i in 1:length(filter.sample_bins)) {
		temp_fsval=filter.sample_bins[i]
		if(temp_fsval==1 || temp_fsval==length(group_indicators) || temp_fsval %in% partial_bins) {
			new_filter.sample_bins=c(new_filter.sample_bins,temp_fsval)
			}
		}
	filter.sample_bins=new_filter.sample_bins
	}
print(filter.sample_bins)
num_combos=length(percentage_bins)*length(filter.count_bins)*length(filter.sample_bins)
print(paste("Num combos : ",num_combos))
run_filter_counter=0

for(fci  in 1:length(filter.count_bins)) {
	filter.count=filter.count_bins[fci]
	#if 1 is min, 5 is default, then consider loop 1....10
	cat(paste("filter.count\t",filter.count,"\n",sep=""))
	for(fsi in 1:length(filter.sample_bins)) {
		filter.sample=filter.sample_bins[fsi]
		cat(paste("filter.sample\t\t",filter.sample,"\n",sep=""))
		#allow filter.sample to range from 1 to the total number of samples
		for(p_i in 1:length(percentage_bins)) {
			percentage=percentage_bins[p_i]
			cat(paste("percentage\t\t\t",percentage,"\n",sep=""))
			print("BEFORE FILTER")
			print(dim(CircRNACount))
			print(paste("Running a filter ; filter.count=",filter.count," ; filter.sample=",filter.sample," ; percentage=",percentage," ; run_filter_counter=",run_filter_counter,sep=""))
			run_filter_counter=run_filter_counter+1
			CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount, Nreplicates = Nreplicates, filter.count = filter.count,filter.sample=filter.sample,percentage=percentage)
			print("AFTER FILTER")
			print(dim(CircRNACount_filtered))
			num_rows_after_filter=dim(CircRNACount_filtered)[1]
			if(num_rows_after_filter==0) {
				print("No data to test!")
				next
				}
			CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]
			LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]
			print("Proceed to test")

			alpha_bins=c(1.0)
			for(ai in 1:length(alpha_bins)) {
				alpha=alpha_bins[ai]
				print(paste("Running test with alpha ",alpha))
				test=Circ.test(CircRNACount_filtered,LinearCount_filtered,CircCoordinates_filtered,group=group_indicators,alpha=alpha,circle_description = c(1:3))
				print("circ.test.done")
				#print(test$summary_table)
				num_candidates=dim(test$summary_table)[1]
				print(paste("num_candidates is ",num_candidates," for setting filter.count",filter.count,", filter.sample",filter.sample,", and percentage",percentage))
				if(num_candidates==0) {
					next
					}
				out_file_root=paste("filtercount",filter.count,"filtersample",filter.sample,"percentage",percentage,"alpha",alpha,sep="_")
				#if(num_candidates>=1) {
				#	#need to plot
				#	pdf_outfile=paste(out_file_root,".plots.pdf",sep="")
				#	pdf(pdf_outfile)
				#	for (i in rownames(test$summary_table)) {
				#		Circ.ratioplot(CircRNACount_filtered,LinearCount_filtered,CircCoordinates_filtered,
				#			plotrow=i,groupindicator1=group_indicator_names,lab_legend='CF Status',gene_column= 4)
				#		}
				#	dev.off()
				#	}
				test_whole_df=getTestResultsAsDF(test)
				outfile=paste(out_file_root,".tsv",sep="")
				print(paste("Writing output to ",outfile))
				write.table(test_whole_df,outfile,sep="\t",row.names = FALSE)
				}
			
			}		
		}

	}
