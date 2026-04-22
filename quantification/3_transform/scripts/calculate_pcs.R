library(tidyverse)
library(data.table)
library(PCAtools)

########## PARSE COMMAND LINE ARGUMENTS ##########
option_list <- list(
    optparse::make_option(c("--bed_file"), type="character", default=NULL,
                        help="BED file from transform.py"),
    optparse::make_option(c("--output_prefix"), type="character", default=NULL,
                        help="Output prefix for PC file")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

bed_file <- opt$bed_file
prefix <- opt$output_prefix
phenotype_pcs_out <- paste0(prefix,'_phenotype_PCs.tsv')

message(paste0('Writing to: ', phenotype_pcs_out))

######## FUNCTIONS ########
compute_pcs <- function(expression_df){
    subsetted_expression_dat <- expression_df %>% select(-c(1,2,3,4))
    pca_standardized <- PCAtools::pca(subsetted_expression_dat)
    n_pcs <- chooseGavishDonoho(subsetted_expression_dat, 
                                var.explained = pca_standardized$sdev^2, 
                                noise = 1)
    message(paste0('Using ', n_pcs,' PCs'))
    
    pca_out <- pca_standardized$rotated %>% 
       data.frame() %>%
       select(1:n_pcs) %>% 
       rownames_to_column('ID') %>% 
       mutate(ID = str_remove(ID,'X'))
    
    pca_out
}

####### ANALYSIS BEGIN ########
bed_df <- readr::read_tsv(bed_file)
message('Bed file loaded')

message('Computing PCs')
PCA_data <- compute_pcs(bed_df)

message('Writing PCs')
PCA_data %>% write_tsv(phenotype_pcs_out)