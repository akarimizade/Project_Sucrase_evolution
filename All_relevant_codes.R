#This file consists of three parts: 
#The first one is the codes that I used to produce the tree
#The second part is the code that I used to download sequences from NCBI and makes it much easier 
#The third part are some other relevant and useful function  

################################################### Part one phylogeny 
####################
####### fishtree tool
####### https://fishtreeoflife.org/about/
####### https://github.com/jonchang/fishtree
####################


##### 1. install and load libraries 

packages = c("fishtree", "ape", "tidytree", "readxl", "grid", "gridExtra", 
             "ggplot2", "ggtree", "tidyverse", "treeio","taxize","usethis","myTAI","tidyverse")
package.check = lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

############### 2.Run before starting 
FontSize = 2.5
AxisTxFontSizeSize_s = 2.5
AxisTxFontSizeSize = 2
AxisTitleFontSizeSize = 7
Factor_mmtoin = 0.0393701
Width_HalfCol = 85*Factor_mmtoin 

openfile <- function(filepath){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  type = tolower(os);
  switch(os, Linux = system(paste0("xdg-open ", filepath)),
         Windows = system(paste0("open \"", filepath, "\"")),
         osx = system(paste0("open \"", filepath, "\""))
  )}




###################### This code is here to help find intersection of species that are fund in Robosky and NCBI
#### Before starting download NCBI list of Actinopterygian species
phy <- fishtree_phylogeny()
treespname <- phy$tip.label # extract only species names

csv_file<-read.delim(choose.files())
vector_organism_name<-csv_file$X.Organism.Name
vector_phy<-phy$tip.label
vector_organism_name_edited<-gsub(" ","_",vector_organism_name)

new<-final_species_order[2]
colnames(new)<-c()
new<-unlist(new)
#list of the species that are found in both the CSV file and in the phylogeny
intersect1<-intersect(vector_phy,vector_organism_name_edited)

phy_filtered<-fishtree_phylogeny(species=intersect1)
write.tree(phy_filtered,"phylogent.tre")

#open the tree that has been made
tree_f<-read.tree(choose.files())

#Finding the order of the species function
#This function takes the species names as a vector and will give you whatever taxonomic info that you want to put on the tree 

##find order and add it to the list
library("taxize")
library("usethis")
library("myTAI")
library("tidyverse")

############################################################################################################### Finding species ranks 
#This is a function to find the the taxonomic ranks of a vector of species from NCBI or itis 
#In the species list add the list of species either in "Mus_musculus" or "Mus musculus" format 
#in the rank_name add the rank you want. For instance class, order, or family 
#in the source_name add either itis or ncbi 
#All tree arguments in the function should be strings 
#If the function cannot find the asked information it will give you either -1 or 0 

order_of_the_species<-function(species_list,rank_name,source_name){
  sp<-gsub("_"," ",species_list)
  vector_species<-c("Arabisopsis")
  vector_order<-c("first_one")
  data_frame_order_name<-data.frame(species_name=vector_species,species_order=vector_order)
  for(i in sp){
    test<-myTAI::taxonomy( organism = i
                           ,db       = source_name,
                           output   = "classification" )
    if(dim(test)[1]<=1&dim(test)[2]<=1){
      test<-data.frame(spe=-1,rank=rank_name)
    }
    check<-filter(test,rank==rank_name)
    check<-unlist(check)
    if(length(check)==0){
      test<-data.frame(spe=0,rank=rank_name)
      
    }
    
    test1<-dim(test)
    
    
    
    
    v1<-i
    v2<-filter(test, rank==rank_name)[1]
    colnames(v2)<-c()
    species_order2<-data.frame(species_name=v1,species_order=v2)
    
    data_frame_order_name<-rbind(data_frame_order_name,species_order2)
    
  }
  data_frame_order_name<-data_frame_order_name[-1,]
  return(data_frame_order_name)
}



########################################################################################################

#example
data_frame_order2<-order_of_the_species(species_list=c("Homo_sapiens","Mus_musculus","Gallus_gallus"),rank_name = "order",source_name = "itis"); data_frame_order2 

#Writing the result 
write_xlsx(data_frame_order2,"species_order.xlsx")


#Plot out and add the label
#First read the species_order_excel_file
data_frame_order3<-read.delim(choose.files())
graph_path = paste0(getwd(), "\\","species_tree12", ".pdf")       # all data
pdf(graph_path, width = Width_HalfCol*4, height = Width_HalfCol*20, pointsize = AxisTxFontSizeSize, onefile = TRUE)
plot(tree_f, show.tip.label = TRUE) 
#to add the species order label
tiplabels(data_frame_order3$species_order,adj = 1.5, cex = 0.8, font = 10, bg="green")

dev.off()

openfile(graph_path)



########################################################## Part 2: NCBI predicted sequences 

library(tidyverse)
library("writexl")
library(dplyr)
#change format function
#Some times each sequence in a fasta file have multiple lines after the accession 
#It is annoying because it makes the code hard to do analysis on 
#This function is used to format any kind of fasta or text data file with accessions to become in the format of: First line accession and the following line the sequence(in one line)
change_format<-function(PATH){
  file_database<-PATH
  file<-readLines(file_database)
  file_result2<-c()
  num<-which(startsWith(file,">"))
  number<-length(num)-1
  if(number==0){
    firstline<-file[num]
    num2<-num+1
    to_collapse<-file[num2:length(file)]
    lastline<-str_c(to_collapse,collapse = "")
    file_result2<-paste(sep="\n",firstline,lastline)
    writeLines(file_result2,"file_result2.txt")
    file_result2<-readLines("file_result2.txt")
    unlink("file_result2.txt")
    return(file_result2)
  }else{
    for (i in 1:number){
      n<-i+1  
      number0<-num[i]
      number1<-num[i]+1
      number2<-num[n]-1
      firstline<-file[number0]
      file_collapsed<-file[number1:number2]
      line_to_line<-str_c(file_collapsed,collapse="")
      file_result<-paste(sep="\n",firstline,line_to_line)
      file_result2<-append(file_result2,file_result)
    }
    find_it<-num[length(num)]
    first_last_line<-file[find_it]
    find_it_2<-find_it+1
    last<-length(file)
    last_one<-file[find_it_2:last]
    last_one<-str_c(last_one,collapse = "")
    last_lines<-paste(sep="\n",first_last_line,last_one)
    file_result2<-append(file_result2,last_lines)
    writeLines(file_result2,"file_result2.txt")
    file_result2<-readLines("file_result2.txt")
    unlink("file_result2.txt")
    return(file_result2)
  }
}

#To run the code on a fasta or text file 
file_for_me<-change_format(choose.files())
#Number of sequences in the file 
number_of_acessions<-length(file_for_me)/2


#Give me the accessions of a downloaded fasta file from NCBI (using the option coding sequences accession files) ##required to run the last function 
give_the_accession<-function(file_for_me_function){
  accession2<-c()
  number_of_acessions<-length(file_for_me_function)/2
  for (i in length(number_of_acessions)) {
    n<-i-1
    n<-2*n+1  
    split1<-unlist(str_split(file_for_me_function[n],"cds"))[1]
    splitted<-unlist(str_split(split1,""))
    m<-length(splitted)-1
    accession<-splitted[6:m]
    
    accesssion<-str_c(accession,collapse = "")
    accession2<-append(accesssion,accession2)
    
  }
  return(accession2)
  
}

#example
give_the_accession(readLines(choose.files()))




#A function for downloading the predicted sequences and formatting them: the last two functions must have been executed
setwd(choose.dir())
name_gene<-c()
accession3<-c()
name_species<-c()
PATH1<-choose.files() #Address of the sequence file downloaded from NCBI 

save_the_file<-function(species_name,gene_name,PATH=PATH1){
  species_name<-gsub(" ","_",species_name)
  gene_name<-gsub(" ","_",gene_name)
  file_for_me2<-change_format(PATH)
  accession1<-give_the_accession(file_for_me2)
  species_name1<-paste(">",species_name,sep="")
  first_line<-paste(paste(species_name1,gene_name,sep="_"),accession1,sep = ".")
  second_line<-file_for_me2[2]
  
  text1<-c(first_line,second_line)
  
  name_of_file<-paste(species_name,gene_name,"fasta",sep = ".")
  if(file.exists(name_of_file)){
    return("File exists")    
    
  }
  else{
    new_file<-file(name_of_file)
    writeLines(text1,new_file)  
  }
  
  unlink(PATH)
  name_species<<-append(name_species,species_name,)
  name_gene<<-append(name_gene,gene_name)
  accession3<<-append(accession3,accession1)
  data_frame_gene<<-data.frame(name_species,name_gene,accession3)
}



#run this
file.exists("sequence.txt")


#copy the name of the species in to your clipboard or right the species name down as the first argument instead of read_clipboard 
readClipboard()


save_the_file(readClipboard(),"Sucrase_isomaltase")

write_xlsx(data_frame_gene,"excel_sequences.xlsx")


##################### Part three: relevant functions 


#puts together a list of separated files containing sequences into one file 
put_them_altogether<-function(address_files=choose.dir(),output=choose.dir()){
  append_vector<-c()
  names_files<-list.files(address_files)
  for(i in names_files){
    x<-readLines(i)
    append_vector<-append(x,append_vector)
    
  }
  to_write<-str_c(collapse = "\n",append_vector)
  file_output<-file(paste(output,"\\","output.txt",sep=""))
  
  writeLines(to_write,file_output)
}

put_them_altogether()


######### Running Muscle, trimal, and IQtree 

#write the code to run in the command prompt in to my clipboard
##################Muscle
#Define the exe path

muscle_exe <- "D:\\Sucrase_evolution\\Muscle\\muscle5.1.win64.exe"


# Define the input and output files
input_file_1 <- "D:\\Sucrase_evolution\\predicted_sequences\\final_sequences\\for_align_sequences_translation.fasta"

# Define the Muscle command
muscle_cmd <- paste("-align", input_file_1, "-output", "for_align_sequences_alignment.fasta",sep=" ")

# Run the Muscle command 
make_ready_muscle<-function(){
  writeLines(paste(muscle_exe,muscle_cmd, sep = " "),"temp.txt")
  temp<-readLines("temp.txt")
  writeClipboard(temp)
  unlink("temp.txt")
  
}

make_ready_muscle()


################## Trimal
# Define the path to the Trimal executable
trimal_exe <- "C:\\Users\\ORIGINAL COMPUTER\\Downloads\\Compressed\\trimal.v1.2rev59\\trimAl\\bin\\trimal.exe"

# Define the input and output files
input_file <- "D:\\Sucrase_evolution\\predicted_sequences\\final_sequences\\Small_tree\\for_align_sequences_alignment.fasta"
output_file <- "D:\\Sucrase_evolution\\predicted_sequences\\final_sequences\\Small_tree\\for_align_sequences_alignment_trimal.fasta"

# Define the Trimal command
trimal_cmd <- paste(" -in ", input_file, " -out ", output_file, "-noallgaps",sep=" ")

# Run the Trimal command using system2()
system2(trimal_exe, trimal_cmd, stdout = TRUE, stderr = TRUE)



########################## IQ tree
#running IQ tree
iqtree_exe <- '"C:\\Users\\ORIGINAL COMPUTER\\Downloads\\Compressed\\iqtree-2.2.0-Windows\\iqtree-2.2.0-Windows\\bin\\iqtree2.exe"'


# Define the input and output files
input_file <- "D:\\Sucrase_evolution\\predicted_sequences\\final_sequences\\Small_tree\\for_align_sequences_alignment_trimal.fasta"

# Define the Iqtree command
iqtree_cmd <- paste("-s", input_file, "-m", "TEST", "-alrt","1000","-bb", "5000","-bcor","0.9")

# Run the Iqtree command 
make_ready_IQ_tree<-function(){
  writeLines(paste(iqtree_exe,iqtree_cmd, sep = " "),"temp.txt")
  temp<-readLines("temp.txt")
  writeClipboard(temp)
  unlink("temp.txt")
  
}
#Copy the code in to my clipboard
make_ready_IQ_tree()


####################Function for p values to work with Maggie's code
#give this function the name of the two ancestors that you want to compare: The result is a data frame that the third column "result" is the p.value
test_to_perform<-function(a0,b0){
  a1<-subset(dose_l_ori_SV2,(dose_l_ori_SV2$Enzyme==a0|dose_l_ori_SV2$Enzyme==b0)  & dose_l_ori_SV2$Ligand=="Maltose");a1<-a1[order(a1$Enzyme),] 
  a2<-subset(dose_l_ori_SV2,(dose_l_ori_SV2$Enzyme==a0|dose_l_ori_SV2$Enzyme==b0)  & dose_l_ori_SV2$Ligand=="Sucrose");a2<-a2[order(a2$Enzyme),] 
  num1<-as.numeric(table(a1$Enzyme))[1]
  num2<-as.numeric(table(a1$Enzyme))[2]
  num3<-num1+1
  num4<-num1+num2
  a1<-t.test(a1$AUC[1:num1],a1$AUC[num3:num4])$p.value
  a2<-t.test(a2$AUC[1:num1],a2$AUC[num3:num4])$p.value
  s<-paste(a0," vs ",b0)
  name_of_subset<-c(s,s)
  kind_of_sugar<-c("Maltose","Sucrose")
  result<-c(a1,a2)
  return(data.frame(name_of_subset,kind_of_sugar,result))
}
#a column can be added using this function that is TRUE if the p value is singificant
add_column <- function(column) {
  new_column <- column < 0.05
  new_column <- ifelse(new_column, "TRUE", "FALSE")
  return(new_column)
}

#comparing all combinations
vec<-c("N1","N11","N31","N32","N33","N49del","N79","N83","N84","N113")

for (n in 1:45){
  vec_name<-(combn(vec,2))[,n]
  c2<-test_to_perform(vec_name[1],vec_name[2])
  c1<-rbind(c1,c2)
}
c1<-c1[-1,]
c1<-c1[-2,]
c3<-c1

c3$significant <- add_column(c3$result)
#calculating the adjusted p value based on Holm correction
adjusted_p_value<-data.frame(p.adjust(c1$result, method = "holm", n = length(c1$result)))

c3<-cbind(c1,adjusted_p_value)
colnames(c3)[5]<-"adjusted.p.value"
c3$significant_adjusted<-add_column(c3$adjusted.p.value)
colnames(c3)[6]<-"significant_adjusted"
#change the order based on significant
c3<-c3[order(c3$significant),] 



