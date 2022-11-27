#!/usr/bin/env Rscript

inFile=""
fcCol="log2FoldChange"
fcCutPoint=0.0
SymCol="GENESYMBOL"
ThreshCol="pvalue"
threshVal=0.01


print("Arguments (required) <inFILE> : optional with Defaults : <FC_COL> <FC_CUT> <SYMCOL> <THRESHCOL>")
args <- commandArgs(trailingOnly = TRUE)

if (length(args)<1) {
  stop("inFile is required!")
}
inFile=args[1]
print(paste("Using inFile ",inFile))
#############################################
if(length(args)>=2) {
fcCol=as.character(args[2])
}
print(paste("Using FC_COL ",fcCol))
#############################################
if(length(args)>=3) {
fcCutPoint=as.numeric(args[3])
}
print(paste("Using FC_CUT ",fcCutPoint))


#############################################
if(length(args)>=4) {
SymCol=as.character(args[4])
}
print(paste("Using SYMCOL ",SymCol))


#############################################
if(length(args)>=5) {
ThreshCol=as.character(args[5])
}
print(paste("Using THRESHCOL ",ThreshCol))



#GET FILE SEPARATOR
cpsSep=""
if("CPSSEP" %in% names(Sys.getenv())) {
	cpsSep=Sys.getenv()["CPSSEP"]
        print(paste("Using separator ***",cpsSep,"***",sep=""))
} else {
stop("Error, require environment variable 'CPSSEP' to be set and hold the separator ; try e.g. CPSSEP=$'\\t'")
}


# Read input file
df = read.csv(inFile, header=TRUE,sep=cpsSep)
print(paste("Loaded ",dim(df)[1]," rows from the file."))
print("Found column names ")
print(colnames(df))

#verify names are present in columns
if(!(fcCol %in% colnames(df))) {
    stop(paste("Did not find FC_COL column ",fcCol," in the column names!"))
}
if(!(SymCol %in% colnames(df))) {
    stop(paste("Did not find SYMCOL column ",SymCol," in the column names!"))
}
if(!(ThreshCol %in% colnames(df))) {
    stop(paste("Did not find THRESHCOL column ",ThreshCol," in the column names!"))
}

#load KEGG DF
sym_to_kegg_gid_df=0
if("KEGG_CP_FILE" %in% names(Sys.getenv())) {
	sym_to_kegg_gid_file=Sys.getenv()["KEGG_CP_FILE"]
	sym_to_kegg_gid_df=read.csv(sym_to_kegg_gid_file)
        print(paste("Using KEGG file ",sym_to_kegg_gid_file,sep=""))
} else {
stop("Error, require environment variable 'KEGG_CP_FILE' to be set for KEGG file lookup !")
}


#load the libraries
suppressMessages(library(biomaRt))
suppressMessages(library(DOSE))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(ggplot2))
suppressMessages(library(enrichplot))
suppressMessages(library(clusterProfiler))
KEGG_ORGANISM=""
if("KEGG_ORGANISM" %in% names(Sys.getenv())) {
	KEGG_ORGANISM=as.character(Sys.getenv()["KEGG_ORGANISM"])
        print(paste("Using KEGG_ORGANISM ",KEGG_ORGANISM,sep=""))
} else {
stop("Error, require environment variable 'KEGG_ORGANISM' (e.g. mmu for mouse or hsa for human) to be set for KEGG data lookup !")
}
GO_DB=""
if(KEGG_ORGANISM=="mmu") {
	GO_DB=org.Mm.eg.db
	print("Using GO DB org.Mm.eg.db")
} else if(KEGG_ORGANISM=="hsa") {
	GO_DB=org.Hs.eg.db
	print("Using GO DB org.Hs.eg.db")
} else {
	stop(paste("UNKNOWN KEGG organism ",KEGG_ORGANISM," ! expected 'mmu' or 'hsa'"))
}

print("Now scanning/dropping for ensemble names...")
invwhich <- function(indices, totlength) is.element(seq_len(totlength), indices)
indices_of_ensemble=grep("^ENS[A-Z]+\\d{11}",df[,SymCol],perl=TRUE)
if(length(indices_of_ensemble)>0) {
	ens_ids_found=df[indices_of_ensemble,SymCol]
	print("Found the following ensembl names")
	print(ens_ids_found)
	print("Dim before drop")
	print(dim(df))
	df=df[!invwhich(indices_of_ensemble,dim(df)[1]),]
	print("Dim after drop")
	print(dim(df))	
	}


threshVals=c(0.001,0.005,0.01,0.05,0.1)
print("To run ClusterProfiler on these p-value thresholds :")
print(threshVals)
for(t_idx in 1:length(threshVals)) {

	#cycle the threshal among the pre-set values above
	threshVal=threshVals[t_idx]
	print("###################################################")
	print(paste("Now using significance threshold value ",threshVal))

	# find significant rows (based on threshold and its cutoff)
	naThreshIndics=is.na(df[,ThreshCol])
	df_this_thresh=df[!naThreshIndics,]
	sig_rows=df_this_thresh[df_this_thresh[,ThreshCol]<=threshVal,]
	print(paste("Based on threshold ",threshVal," found ",dim(sig_rows)[1]," significant rows!"))
	head(sig_rows)

	#from a given DF, have a function to extract the gene symbols using the symbol column name
	getGeneNames=function(the_df,sym_col) {
	    genes=the_df[,sym_col]
	    genes=na.omit(genes)
	    genes=unique(genes)
	    gene_list = sort(genes, decreasing = TRUE)
	    return(gene_list)
	}


	# find the up/down regulated data
	downRegSig=sig_rows[sig_rows[,fcCol]<(-abs(fcCutPoint)),]
	UpRegSig=sig_rows[sig_rows[,fcCol]>abs(fcCutPoint),]
	print(paste("By additionally filtering on fold change ; less than the (",fcCutPoint,") cutoff : ",dim(downRegSig)[1]))
	print(paste("By additionally filtering on fold change ; greater than the (", fcCutPoint,") cutoff : ",dim(UpRegSig)[1]))


	singleSymToVecOfKeggs=function(sym,sym_to_kegg_gid_df) {
	    lone=sym_to_kegg_gid_df[sym_to_kegg_gid_df$SYMBOL==sym,2]
	    temp_vec=c()
	    for(l in 1:length(lone)) {
		test_skip_short=nchar(as.character(lone[l]))<=4
		if(length(test_skip_short)==0 || is.na(test_skip_short) || test_skip_short) {
		    #if NA or XXX: or less characters go to next
		    next
		}
                kegg_db_pattern="^[a-z]+:"
		if(length(temp_vec)>=1) {
		temp_vec=c(sub(kegg_db_pattern,"",lone[l]),temp_vec)
		    } else {
		    temp_vec=c(sub(kegg_db_pattern,"",lone[l]))
		}
	    }
	    return(temp_vec)
	}
	symsToVecOfKeggs=function(syms,sym_to_kegg_gid_df) {
	    allKeggs=c()
	    for(s in 1:length(syms)) {
		sym=syms[s]
		the_keggs=singleSymToVecOfKeggs(sym,sym_to_kegg_gid_df)
		allKeggs=c(allKeggs,the_keggs)        
	    }
	    return(allKeggs)
	}



	# for each of the 3 modes , have a function to plot the bubles for a given gene_list
	plotBubbles=function(genes,name_root) {
	    #print(paste("Now plotting for "))
	    #print(gene_list)
            write(genes,paste(name_root,".input_genes_list.txt",sep=""))
	    classes=c("BP","CC","MF")
	    #classes=c("BP")
	    for(c_i in 1:length(classes)) {
		this_cat=classes[c_i]
		#print(paste("now processing category ",this_cat," with input size ",length(genes)))
		tryCatch({
				go_enrich <- enrichGO(gene = genes, OrgDb = GO_DB, 
						  keyType = "SYMBOL", readable = FALSE, ont = this_cat, 
						  pAdjustMethod="BH", pvalueCutoff = 1.0, 
						  minGSSize = 5,
						  maxGSSize = 2000,
						  qvalueCutoff = 1.0)        
				enrich_table_file=paste(name_root,this_cat,"table.txt",sep=".")
				tiff_file_name=paste(name_root,this_cat,"tiff",sep=".")
				print(paste("Attempt to write to files ",tiff_file_name," and ",enrich_table_file,sep=""))
				write.table(go_enrich,file = enrich_table_file,sep="\t",row.names = FALSE)
				ggsave(tiff_file_name, plot = (dotplot(go_enrich, showCategory=10,
								label_format = function(x) stringr::str_wrap(x, width=60))), 
					   device = "tiff", path = NULL, scale = 1, width = 15, 
					   height = 8, units=("in"), dpi = 100, limitsize = TRUE)
			},error=function(cond) {
				err_msg=paste("Error in running/plotting for GO enrich with class=",this_cat," name root = ",name_root," table=",enrich_table_file," tiff=",tiff_file_name,sep="")
				print(err_msg)
				print(cond)
				message(cond)
				message(err_msg)
			},warning=function(cond) {
				warn_msg=paste("Warning in running/plotting for GO enrich with class=",this_cat," name_root = ",name_root," table=",enrich_table_file," tiff=",tiff_file_name,sep="")
				print(warn_msg)
				message(cond)
				message(warn_msg)
			})
	    }


	tryCatch({
		#do the work for KEGG
		kegg_ids=symsToVecOfKeggs(genes,sym_to_kegg_gid_df)
		kegg_enrich <- enrichKEGG(gene         = kegg_ids,
			organism     = KEGG_ORGANISM,
			pvalueCutoff = 0.05,
			keyType = "kegg")
		enrich_table_file=paste(name_root,"KEGG","table.txt",sep=".")
		tiff_file_name=paste(name_root,"KEGG","tiff",sep=".")
		print(paste("To write KEGG data to files ",tiff_file_name," and ",enrich_table_file,sep=""))
		write.table(kegg_enrich,file = enrich_table_file,sep="\t",row.names = FALSE)
		ggsave(tiff_file_name, plot = (dotplot(kegg_enrich, showCategory=10,
						label_format = function(x) stringr::str_wrap(x, width=60))), 
			   device = "tiff", path = NULL, scale = 1, width = 15, 
			   height = 8, units=("in"), dpi = 100, limitsize = TRUE)
		},
		error=function(cond) {
			err_msg=paste("Error in running/plotting for KEGG enrich with class=",this_cat," name root = ",name_root," table=",enrich_table_file," tiff=",tiff_file_name,sep="")
			print(err_msg)
			print(cond)
			message(cond)
			message(err_msg)
		},
		warning=function(cond) {
			warn_msg=paste("Warning in running/plotting for KEGG enrich with class=",this_cat," name_root = ",name_root," table=",enrich_table_file," tiff=",tiff_file_name,sep="")
			print(warn_msg)
			print(cond)
			message(cond)
			message(warn_msg)
		})
		



	}
	#define the tiff file name root
	name_root=paste(inFile,"sig","thresh",threshVal,sep=".")
	name_root_up=paste(inFile,"sig.up","thresh",threshVal,sep=".")
	name_root_down=paste(inFile,"sig.down","thresh",threshVal,sep=".")

	###############################
	#acquire the gene name lists
	#any
	pvalThreshGenes=getGeneNames(sig_rows,SymCol)
	print(paste("Found ",length(pvalThreshGenes)," unique gene names sig. diff expr."))
	print(pvalThreshGenes)
	#down
	downRegGenes=getGeneNames(downRegSig,SymCol)
	print(paste("Found ",length(downRegGenes)," unique gene names sig. diff expr DOWN."))
	print(downRegGenes)
	#up
	upRegGenes=getGeneNames(UpRegSig,SymCol)
	print(paste("Found ",length(upRegGenes)," unique gene names sig. diff expr UP."))
	print(upRegGenes)

	#perform the analyses and plot bubbles!
	#print("Analyzing for significantly differentially expressed genes")
	if(length(pvalThreshGenes)>0) {
		plotBubbles(pvalThreshGenes,name_root)
		} else { 
		print("Skipping for differentially expressed genes because none found!") 
	}
	#print("Analyzing for significantly differentially expressed genes (down)")
	if(length(downRegGenes)>0) {
		plotBubbles(downRegGenes,name_root_down)
		} else {
		print("Skipping for differentially expressed genes (down) because none found!")
		}
	#print("Analyzing for significantly differentially expressed genes (up)")
	if(length(upRegGenes)>0) {
		plotBubbles(upRegGenes,name_root_up)
		} else {
		print("Skipping differentially expressed genes (up) because none found!")	
		}



}


