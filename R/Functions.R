
########################################################################
#
# Check, install, and load multiple packages from cran and bioconductor
# 
# ######################################################################



cran.libraries <- function(...){
        # For cran libraries
        list.of.cran.packages<- c("easypackages",...)
        new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)) install.packages(new.packages)
        library(easypackages)
        libraries(list.of.cran.packages)
        }



bio.libraries <- function(...){
        # For bioconductor libraries
        if(!("easypackages" %in% installed.packages()[,"Package"]))
                install.packages("easypackages")
        source("https://bioconductor.org/biocLite.R")
        list.of.bio.packages <- c(...)
        new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)) biocLite(new.packages)
        library(easypackages)
        libraries(list.of.bio.packages)
        }

# #####################################################################
# 
#  detect OS and set enviroment
#
# ####################################################################

initialize.project <- function(x){
        if (Sys.info()[['sysname']]=="Darwin"){
                WD <- paste0("/Users/yah2014/Dropbox/Public/Olivier/R/",x,"/datasets")
                setwd(WD);print(getwd());print(list.files())}
        if (Sys.info()[['sysname']]=="Windows"){
                WD <- paste0("C:/Users/User/Dropbox/Public/Olivier/R/",x,"/datasets")
                setwd(WD);print(getwd());print(list.files())}
        if (Sys.info()[['sysname']]=="Linux"){
                WD <- paste0("/home/yah2014/Dropbox/Public/Olivier/R/",x,"/datasets")
                setwd(WD);print(getwd());print(list.files())}
}

#===convert srt file to txt===
convert_srt_to_txt <- function(x){
        rldm <- read.delim(x,header = FALSE,stringsAsFactors =F)
        n <- seq(3,nrow(rldm),3)
        script <- rldm[n,]
        text <- paste(unlist(script), collapse =" ")
        text
}
setwd("/Users/yah2014/Downloads/Outro Subtitles")
list_files = list.files()
list_files = grep("lang_en_",list_files,value = T)
(list_files_0 <- list_files)

# add 0 to files
less_then_10 <- as.numeric(sub(' - .*', '', list_files)) <10
(list_files_0[less_then_10] = paste0("0",list_files[less_then_10]))
file.rename(from = list_files, to = list_files_0)

list_files = sort(list_files_0)

Text = convert_srt_to_txt(list_files[1])
for(i in 2:length(list_files)){
        Text = c(Text,"\n\n\n\n",convert_srt_to_txt(list_files[i]))
}
write.table(Text,"subtitles.txt")
#=====Clean memory======================
#WGCNA::collectGarbage()

#=== Merge all files in one folder ===
MergeAllTxt1 <- function(folder = "",select = NULL,header = FALSE,
                         method = "full_join"){
        # folder = "Danwei", for example. 
        # folder = "" will merge all txt files in current folder
        # select = 1:2 columns for example
        list_files <- list.files(paste0("./",folder),".txt")
        if(length(list_files) < 1) {return("No file in the folder")}
        print(head(list_files))
        counts <- read.delim2(paste0("./",folder,"/",list_files[1]),
                              header = header, fill = TRUE, skipNul = TRUE)
        if(!is.null(select)) counts <- counts[,select]
        colnames(counts)[1] <- "Gene"
        
        if(method == "full_join"){Mergefun = dplyr::full_join}
        if(method == "inner_join"){Mergefun = dplyr::inner_join}
        
        for(i in 2:length(list_files)){
                temp <- read.delim2(paste0("./",folder,"/",list_files[i]),
                                    header,fill = TRUE, skipNul = TRUE)
                if(!is.null(select)) temp <- temp[,select]
                colnames(temp)[1] <- "Gene"
                counts <- Mergefun(counts,temp,by="Gene")
                print(paste0("Complished file ",i," in ",length(list_files)))
        }
        print(dim(counts))
        return(counts)
}

#ALL.counts <- MergeAllTxt1(folder = "",header = TRUE,method = "inner_join")


#===================================================
# rename all extentions from pileup to txt file in sub folder
Rename <- function(folder, from = ".pileup", to = ".txt"){
        #===change file extention======
        # folder = "Theunissen_2016", for example
        list_files <- list.files(paste0("./",folder),from)
        if(length(list_files) < 1) {return("No file in the folder")}
        print(head(list_files))
        
        for(i in 1:length(list_files)){
                File.Name <- tools::file_path_sans_ext(list_files[i])
                file.rename(from = paste0("./",folder,"/",list_files[i]),
                            to = paste0("./",folder,"/",File.Name,to))
        }
}
#Rename("Takashima_Y_2014", from = ".bam", to = ".txt")
