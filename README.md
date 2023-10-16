# LD Score regression analysis
# Tool background

This tool was designed to study the heritability and the genetic correlation from GWAS summary statistics to other traits and genetic functional annotations. All the credit of this bioinformatic tool is in the following github link: <https://github.com/bulik/ldsc>. The important point of this tool is that it takes into account the linkage disequilibrium (LD) blocks of the genome, which are inherited together. In this sense, is not useful to find an enrichment of variants located in a massive LD block when not taken into account the size of this genomic region, as it will lead to misinterpretation of the data. Instead, including this as a covariate helps us to distinguish properly the effects of this GWAS genetic associations.

# First steps: Installing the package

Firstly, you will need to download the github repository:

```{bash, eval = F}
git clone https://github.com/bulik/ldsc.git
cd ldsc
```

Is important that this package works with python, so you will need the Anaconda Python distribution to be installed in your server. Then you have to create a conda enviroment within the ldsc folder.

```{bash, eval = F}
conda env create --file environment.yml
source activate ldsc
```

**\*Note:** if you exit that environment, to start it again you will need to do the following and log into ldsc environment:

```{bash, eval = F}
conda env update --file environment.yml --prune
```

# Step 1: Prepare the summary statistics of the desired disease/trait

The first step is to make sure the program is able to read your summary statistics. To do this, you will use the munge_sumstats.py. For this, you will require a list of SNPs to be included in the study. At the end, you would like this SNPs to merge well with the SNPs that are included in the posterior regression analysis, so make sure they match as much as possible. Same if you are going to compare the heritability between two diseases/traits, you would probably want the summary statistics files to share many of these SNPs. Each summary statistics would have its particular issues, but this would be an example code line for this part:

```{bash, eval = F}
python munge_sumstats.py --out DiseaseA --merge-alleles SNPsToBeIncluded.txt --N 26679 --sumstats DiseaseA_summstats.txt
```

For more information of this analysis, please refer to: <https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics>

The information munge_sumstats.py needs:

| Variant ID (rs)        | Allele 1 (a1)     | Allele 2 (a2)         | Sample size (N)                                                                                   | P-value (p) | OR (or)                                                                                          |
|------------|------------|------------|------------|------------|------------|
| Usually, the RS number | The effect allele | The non-effect allele | It may vary from SNP to SNP (If you pass it through the previous line you don't need this column) | The pvalue  | The effect of each variant. Can be either OR or Beta. It will be used to calculate the Z scores. |

As you may notice, it is possible to prepare it by yourself, however, if you run the previous line you will make sure it does not give you any unspecified problem in following analysis.

**Little trick** 
If you're struggling with your summary statistics data because it takes too long, try this flag: --chunksize 500000

# Step 2: Prepare the annotations

If you want to assess if there is an enrichment of the associated variants within specific cell functional annotations, you need to calculate the LD scores of this data. In this sense, you will need a .bed file of each annotation. Here is an example how I processed the annotations DNase I hypersensitive sites from this paper: <https://www.nature.com/articles/s41586-020-2559-3>. If you want to use this files, they are available here: <https://zenodo.org/record/3838751#.X69tfEJKg6U>.

This data is distributed as following:

**DHS_Index_and_Vocabulary_hg19_WM20190703.txt -** This is the .bed format file for each studied region and some additional information, such as the number of cell types studied in that particular region and the mean signal values, which does not interest us much. It has over 3 million rows, one per each region.

**DHS_Index_and_Vocabulary_hg19_metadata.gz -** This includes all the information for each cell type included in the study. It has over 700 hundred cells.

**dat_bin_FDR01_hg19.txt.gz -** This is the resulting matrix. Each row is a genomic region and each column is each cell annotation. If there is an annotation of this cell type (column) in a certain region (rows), it will be indicated with a 0. Instead, it will be a 0.

## However, what do we want from these files?

We want a .bed file that contains with the exact regions that had a signal en each cell type. The script is only going to learn the chromosome, the start, the end, and the identifier columns, so we only need this those rows. To do this, I prepared the following script.

```{r, eval = F}
library(data.table) #To import the datasets faster

setwd("~/MPRA/Enrichmentanalyses/ldsc/DHS")
DHSMatrix <- fread(input = "dat_bin_FDR01_hg19.txt.gz", header = F)
DHSIndex <- fread(input = "DHS_Index_and_Vocabulary_hg19_WM20190703.txt.gz", header = T)
DHSMetadata <- fread(input = "DHS_Index_and_Vocabulary_metadata.tsv", header = T, nrows = 733) #The final argument is because it includes an additional empty row

colnames(DHSMatrix) <- paste(DHSMetadata$Organ, '.', DHSMetadata$`Biosample name`, sep = '') #Just to make the matrix more accesible
row.names(DHSIndex$identifier) #Just to make the matrix more accesible

colnames(DHSIndex) <- c('#seqnames', colnames(DHSIndex)[2:ncol(DHSIndex)]) 
#Bedtools to work needs the first line to start with #, so I change it here to not cause any trouble when preparing the annotations.

#Select the cell types you are interested in
Dis <- which(DHSMetadata$`Biological state`!='Cancer') #First filter I used
Age <- which(DHSMetadata$`Growth stage`=='Adult') #Second filter I used

TOKEEP <- Dis[Dis%in%Age] #This vector contains the exact rows you are interested in.

#Prepare a loop to take the information you want
for(i in TOKEEP){
  
  PICOS <- which(DHSMatrix[,..i]==1) #Takes out the rows [genetic regions] of which is cell type [columns] has an annotation from the big matrix
  
  write.table(DHSIndex[PICOS,], 
              file=paste('DHSDatasetsToCompute/', DHSMetadata$Organ[i], '.', DHSMetadata$`Biosample name`[i],'.', DHSMetadata$`Altius Library ID`[i],'.bed', sep = ''), 
              col.names=T, row.names=F, 
              quote=F, sep = '\t') #Rewrites the bed file including only the desired regions and names it properly
}

```

**NOTE:** It is also important to save a file that contains the name of the files you created in the previous loop in each row. This will be used in the following steps to do the regressions of all the annotations automatically.

```{r, eval = F}
filenames <- paste(DHSMetadata$Organ[TOKEEP], '.', DHSMetadata$`Biosample name`[TOKEEP],'.', DHSMetadata$`Altius Library ID`[TOKEEP], sep = '') #Make sure that the name matches with the name you save the file in the loop!
write.table(filenames, "filenames.txt", col.names = F, row.names = F, quote = F, sep = "\n")
```

This will generate a file as following:

```         
Blood.GM06990.LN1203
Blood.hTH1.LN1222
Mammary.HMEC.LN1376
Blood.GM12878.LN1614
Lung.SAEC.LN1924
```

# **Step 3: Run the LD regression for each annotation through each chromosome**

Now, it is time to calculate the LD score for each cell type annotation. Additionally, this is done for each chromosome, so if you want to examine 100 cell type annotations, you will need to run this command 2,200 times. Instead, I prepared a script that does this in an automatic way. However, it will take long, as each LD score for each chromosome lasts around 40-50 minutes.

## Previous things you need to know about the script:

**\--bed-file DHS/DHSDatasetsToCompute/\$file_name.bed** - This is the .bed file previously generated in Step 2.

**\--bimfile 1000G_plinkfiles/1000G.mac5eur.\${i}.bim** - This is the .bim file from the 1000 genomes project of european population

**\--annot-file DHS_AnnotationFiles/\$file_name.\${i}.annot** - This is the output for the annotations. It takes in consideration both the cell type and the chromosme

**\--bfile 1000G_plinkfiles/1000G.mac5eur.\${i}** - Refers to the plink bfiles of the 1000 genomes project from European populations that will be used to

**\--ld-wind-cm** 1 - This is the genomic window to estimate the LD scores. It is measured in cM.

**\--thin-annot** - This is to create your own annotation files from your bed, is mandatory.

**\--print-snps hapmap3_snps/hm.\${j}.snp** -This is to use the hapmap3 SNPs annotation type, which will be determinant to match properly with your previous summary statistics file.

## MakeAnnotations.sh

```{bash, eval = F}

#!/bin/bash

read -p "Please enter the name and the path to the annotations names file: " input_file #This file contains all the datasets you want to prepare for your enrichment analysis

if [ ! -f "$input_file" ]; then #Check if the input file exists, if not, it exists the script.
  echo "Input file '$input_file' not found"
  exit 1
fi
while IFS= read -r file_name; do #This will go through each .bed file name
  if [ -f "DHS/DHSDatasetsToCompute/$file_name.bed" ]; then #Checks if it exists
  echo "Processing file $file_name" 
      for i in {1..22} #This will be going through each autosomic chromosme
      do
      (
        python make_annot.py --bed-file DHS/DHSDatasetsToCompute/$file_name.bed --bimfile 1000G_plinkfiles/1000G.mac5eur.${i}.bim --annot-file DHS_AnnotationFiles/$file_name.${i}.annot
        gzip -c DHS_AnnotationFiles/$file_name.${i}.annot >> DHS_AnnotationFiles/$file_name.${i}.annot.gz  #Gunzip the file to be the output of the next step
        echo "Finished annotation for: $file_name and chromosome ${i}"
      ) & #It is between brackets to perform each iteration in parallel to save time
    done
    wait #This is to make sure it does not start the next step until it finishes this one
  else
    echo "File not found: $file_name"
  fi
  rm DHS_AnnotationFiles/*.annot #This will remove all the intermediate .annot files
  
  for j in {1..22}
  do
  (
    if [ -f DHS_AnnotationFiles/"$file_name.${j}.annot.gz" ]; then
      python ldsc.py --l2 --bfile 1000G_plinkfiles/1000G.mac5eur.${j} --ld-wind-cm 1 --annot DHS_AnnotationFiles/$file_name.${j}.annot.gz --thin-annot --out DHS_AnnotationFiles/$file_name.${j} --print-snps hapmap3_snps/hm.${j}.snp
    else
      echo "File not found: $file_name.${j}.annot.gz"
    fi
  ) & #Same here, it parallelizes the process. This particular step is very time consuming.
  done
  
done < "$input_file"
```

**\*Reminder:** This paths are particular for my analysis, you will have to modify them for yours.

### Disclosure

The script assumes that the SNPs are in the same order. In this case, as we are going to use the same SNPs datasets, we assume that they are. However, this should be taken into account when splitting this step into two different phases: the **make_annot.py** and the **ldsc.py**

# Step 4: Cell type enrichment analysis

In this step, you will use the LD scores generated before to see the correlation bewteen your disease and the annotations you are interested in.

## Previous things you need to know about the script:

**\--ref-ld-chr** - This is the reference ld score you want your summary statistics to be compared.

**\--frqfile-chr** - This is the file of each SNPs frequency. In this case, is the 1000 genomes project data.

**\--overlap-annot** - This is to tell that the categories you used to create the annotation file overlap with eachother.

**\--w-ld-chr** - This is the file that contains the weight of each SNP in the regression.

## LDScore_CellTypeAnnotations.sh

```{bash, eval = F}
#!/bin/bash

#Asks for the summstats file
read -p "Please enter the name and the path to the summary statistic file: " summstat

read -p "Please enter the name and the path to the annotations names file: " input_file #This file contains all the datasets you want to prepare for your enrichment analysis


if [ ! -f "$input_file" ]; then #Check if the input file exists, if not, it exists the script.
  echo "Input file '$input_file' not found"
  exit 1
fi

if [ ! -f "$summstat" ]; then #Check if the input file exists, if not, it exists the script.
  echo "Summary statistics '$summstat' not found"
  exit 1
fi

annotation_values="" #This will contain all the annotations to evaluate

#Loop through all the names of the input file containing the name of the Annotations
whle IFS= read -r name; do
  #Adds the name to the previous variable
  annotation_values="$annotation_values DHS_AnnotationFiles/$name.," < "$input_file" #The point is to include all the chromosomes
done

annotation_values=${annotation_values# } #Removes the initial space

echo "You are going to evaluate the following annotations: $annotation_values"
wait

read -p "Give me an output name for your results: " output

#Construct the command
ldsc_command="python ldsc.py --h2 $summstat --ref-ld-chr $annotation_values --out $output --overlap-annot  --frqfile-chr 1000G_frq/1000G.mac5eur. --w-ld-chr weights_hm3_no_hla/weights. --print-coefficients"

#Execute the analysis
eval "$ldsc_command"

```

**NOTE:** This will give you a single table with each annotation in a row, and the statistics coefficients in each column

# Extra step: Heritability and genetic correlation between diseases and traits

This refers to the analysis performed in the Fig 1 of <https://www.nature.com/articles/s41467-023-36306-5>. This analysis gives you an estimated shared heritability across different diseases and traits. It is easy to calculate and to implement once you have the summary statistics.

This analysis is detailed here: <https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#estimating-heritability-genetic-correlation-and-the-ld-score-regression-intercept>

```{bash, eval = F}
python ldsc.py --rg disease1.summstat,disease2.summstat --ref-ld-chr eur_w_ld_chr --w-ld-chr weights_hm3_no_hla/weights. --out disease1.disease2 
```
**\*NOTE:** Be very careful with the ref-ld-chr, as it is quite frequent that an error: list out of range appear, and it has to be with this parameter.
