library(dplyr)
library(tidyr)
library(optparse)

option_list <- list(
  make_option("--dragen_dir", type = "character",
              help = "path to dragen results"),
  make_option("--output", type = "character", default = ".",
              help = "path to output file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dragen_dir <- opt$dragen_dir
out_fn <- opt$output

tumor_dir <- paste0(dragen_dir,"/tumor_normal/")
normal_dir <- paste0(dragen_dir,"/normal/")

# -------------------------------------------------------------------------
# Tumor coverage
cov_fn <- list.files(tumor_dir, pattern = "*wgs_overall_mean_cov_tumor.csv",
                     recursive = T, full.names = T)
if(length(cov_fn) !=0){
  cov <- lapply(cov_fn, read.csv, header =F)
  names(cov) <- cov_fn
  
  cov_df <- bind_rows(cov, .id = "path")
  
  cov_df1 <- cov_df %>%
    mutate(sample = gsub(".wgs_overall_mean_cov_tumor.csv","",basename(path)))  %>%
    select(-V1) %>%
    rename(coverage = V2) 
} else{
  cov_df1 <- data.frame()
}

# Normal coverage
cov_fn <- list.files(normal_dir, pattern = "*wgs_overall_mean_cov.csv",
                     recursive = T, full.names = T)
cov <- lapply(cov_fn, read.csv, header =F)
names(cov) <- cov_fn

cov_df <- bind_rows(cov, .id = "path")

cov_df2 <- cov_df %>%
  mutate(sample = gsub(".wgs_overall_mean_cov.csv","",basename(path)))  %>%
  select(-V1) %>%
  rename(coverage = V2) 

# Purity and ploidy
purity_fn <- list.files(tumor_dir, pattern = "*.cnv_metrics.csv",
                        recursive = T, full.names = T)
purity1 <- lapply(purity_fn, read.csv, header =T, row.names = NULL)
purity <- lapply(purity1, function(tmp) tmp[grepl("Estimated tumor purity",tmp[,3]),4])
ploidy <- lapply(purity1, function(tmp) tmp[grepl("Overall ploidy",tmp[,3]),4])
ids <- gsub(".cnv_metrics.csv","",basename(purity_fn))

purity_df <- data.frame(sample = ids, purity = unlist(purity), ploidy = unlist(ploidy))

# Peak fragment size
frag_len_fn <- list.files(dragen_dir, pattern = ".fragment_length_hist.csv", recursive = T, full.names = T)
frag_len <- lapply(frag_len_fn, read.csv, skip = 1)
names(frag_len) <- gsub(".fragment_length_hist.csv","",basename(frag_len_fn))

frag_length <- bind_rows(frag_len, .id = "sample")
peak_frag_length <- frag_length %>%
  group_by(sample) %>%
  summarise(insert_size = FragmentLength[which.max(Count)])

QC_df <- bind_rows(cov_df1,cov_df2) %>%
  left_join(purity_df) %>%
  left_join(peak_frag_length) %>%
  distinct()

write.csv(QC_df, paste0(out_fn,"QC_metrics.csv"), row.names =F)
