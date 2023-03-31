##########################################################
################# CELL BARCODE MATCHER ###################
##########################################################

run <- "samplename"
reads <- "25000000"

# input file paths
R1_file_path <- "/path/to/sub_samplename_25000000_R1.fastq"
R2_file_path <- "/path/to/sub_samplename_25000000_R2.fastq"
CB_file_path <- "/path/to/cellbender_cell_barcodes.csv"

# output file paths
df_whitelist_output_path <- "/path/to/df_whitelist_samplename.txt"
reads_per_SB_output_path <- "/path/to/reads_per_SB_samplename.csv"
matcher_summary_output_path <- "/path/to/matcher_summary_samplename.txt"



############### DATA AND PACKAGE SETUP ###############

# load required R packages.
library(ShortRead)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggthemes)  
library(plyr)
library(fuzzyjoin)
library(tidyr)

# read in sequencing and cell barcode data. perform some pre-processing to get desired data format
R2 <- readFastq(R1_file_path) # contains spatial barcode information
R1 <- readFastq(R2_file_path) # contains cell barcode information
CB <- read.csv(CB_file_path) # list of 10x cell barcodes
#CB <- read.table(CB_file_path) # read txt file if cell ranger output not cellbender

colnames(CB) <- "CB_region"
CB$CB_region = substr(CB$CB_region,1,nchar(CB$CB_region)-2) 

R2_seq <- as.data.frame(sread(R2))
R1_seq <- as.data.frame(sread(R1))
colnames(R2_seq) <- "R2_seq"
colnames(R1_seq) <- "R1_seq"
R2_seq$row_position <- row.names(R2_seq)
R1_seq$row_position <- row.names(R1_seq)

# check dimensions of R1 and R2 are the same, plus check dimensions of CB
dim(R2_seq)
dim(R1_seq)
dim(CB)

# read in GEX CB / FB CB matches (first transfer from google bucket)
CB_dictionary <- read.table("/path/to/3M-february-2018.txt") # list of 10x cell barcodes --> https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz

# rename columns 
colnames(CB_dictionary) <- c("GEX_CB", "CB_region") # variant 1 is CB from GEX capture seq, variant 2 is from FB capture seq 



############### PREPARE R2 FOR CB-R2 JOIN ############### 

# Create new column in R2_seq dataframe with the expected 10x CB region that we will use for joining (CB is positions 1 to 16 of read 2). 
R2_seq <- R2_seq %>%
  mutate(CB_region = substr(R2_seq,1,16))

# Assign group ID to each CB_region sequence. Remove duplicates 
R2_seq <- R2_seq %>%                         
  group_by(CB_region) %>%
  dplyr::mutate(ID = cur_group_id())
R2_CB_region <- R2_seq %>%
  select(row_position, CB_region, ID)
CB$cell_number <- row.names(CB)
R2_CB_region_minus_dups <- R2_CB_region[!duplicated(R2_CB_region$ID),]

# Check number of unique CB regions in read 2
dim(R2_CB_region_minus_dups)



############### DICTIONAIRY RECODE ############### 

print("start dictionary recoding:")
Sys.time()

# left join
R2_CB_region_minus_dups_corrected <- left_join(x = R2_CB_region_minus_dups, y = CB_dictionary, by = "CB_region")

# NAs in join results
paste("NAs in left joined GEX_CB:")
sum(is.na(R2_CB_region_minus_dups_corrected$GEX_CB))
sum(is.na(R2_CB_region_minus_dups_corrected$GEX_CB))/length(R2_CB_region_minus_dups_corrected$GEX_CB)

# drop rows with NA 
R2_CB_region_minus_dups_corrected_NA_dropped <- R2_CB_region_minus_dups_corrected %>% drop_na(GEX_CB)

print("finish dictionary recoding:")
Sys.time()



############### GRAB R2_SEQ READS ############### 
#  join R2_seq reads to 10x cell barcodes list by only looking at one representative CB per group

# some column renaming for inner_join
colnames(CB) <- c("CB_region_GEX", "cell_number")
colnames(R2_CB_region_minus_dups_corrected_NA_dropped) <- c("row_position", "CB_region", "ID", "CB_region_GEX")

# join
joined <- inner_join(R2_CB_region_minus_dups_corrected_NA_dropped, CB, by = "CB_region_GEX")

# Count how many of the 10x CBs were matched with reads.
number_recovered_CB <- as.data.frame(joined$cell_number)
number_recovered_CB <- number_recovered_CB[!duplicated(number_recovered_CB$`joined$cell_number`),]
paste("how many 10x CBs matched with reads:")
print(length(number_recovered_CB))



############### GRAB READS CORRESPONDING WITH MATCHED CB FROM R2 AND R1 ############### 

# Join R2-CB grouped_join with original R2_seq file by ID so we recover all reads with CB_region matching 10x CB. 
R2_cells_h2 <- inner_join(joined, R2_seq, by = "ID")

# Print number of reads in R2_seq and number of reads with a matching 10x CB 
dim(R2_seq)[1]
dim(R2_cells_h2)[1]
dim(R2_cells_h2)[1]/dim(R2_seq)[1]*100

# Select columns of interest
R2_cells_h2 <- R2_cells_h2 %>% 
  select(CB_region.x,
         ID,
         CB_region.y,
         cell_number,
         R2_seq,
         row_position.y,
         CB_region_GEX
  )

colnames(R2_cells_h2) <- c("R2_CB_region_seq",
                           "R2_CB_region_ID",
                           "tenx_CB_seq",
                           "tenx_CB_number",
                           "R2_seq",
                           "row_position",
                           "CB_region_GEX")

# Join with R1 to get all R1 and R2 reads with matched CB.
tenx_CB_reads <- inner_join(R1_seq, R2_cells_h2, by = "row_position")

# Break up the read sequence into component parts and create full length SB. 
tenx_CB_reads <- tenx_CB_reads %>%
  mutate(R2_CB_UMI_seq = substr(R2_seq,17,28)) %>%
  mutate(R1_SBa_seq = substr(R1_seq,1,8)) %>%
  mutate(R1_SBb_seq = substr(R1_seq,27,32))
tenx_CB_reads <- tenx_CB_reads %>%
  mutate(R1_SB_concat = paste(R1_SBa_seq, R1_SBb_seq, sep = ""))



############### PREPARE OUTPUT FILE WITH WHITELISTED CELL AND SPATIAL BARCODES  ############### 

# Write a file with nUMI (per CB_SB combination), CB, and SB.
final_df <- tenx_CB_reads %>%
  select(tenx_CB_seq,
         CB_region_GEX,
         R1_SB_concat,
         R2_CB_UMI_seq)

# remove duplicates
final_df_nodups <- final_df %>% distinct() 

# create column with CB_SB concat
final_df_nodups <- final_df_nodups %>% 
  mutate(CB_SB = paste(tenx_CB_seq, R1_SB_concat, sep = "")) 

# count number of unique CB_SB, which should equal number of UMIs per CB_SB because we removed duplicate reads
UMIs_per_CB_SB <- data.frame(table(final_df_nodups$CB_SB)) 
colnames(UMIs_per_CB_SB) <- c("CB_SB", "nUMI") 

# split reads up again
UMIs_per_CB_SB <- UMIs_per_CB_SB %>%
  mutate(cell_bc_10x = substr(CB_SB,1,16)) %>%
  mutate(bead_bc = substr(CB_SB,17,30))

# write table
write.table(UMIs_per_CB_SB, df_whitelist_output_path, quote = FALSE, sep = "\t")



############### PREPARE OUTPUT FILE FOR PUCK MATCHER ############### 

# Write a file with number of UMIs per SB for puck matcher.
reads_per_SB_t <- table(final_df$R1_SB_concat)
reads_per_SB_t <- t(as.matrix(reads_per_SB_t))
write.csv(reads_per_SB_t, reads_per_SB_output_path, quote = FALSE)



############### PREPARE OUTPUT FILE WITH SUMMARY STATISTICS ############### 

# Write a file with summary statistics for the 10x matcher

# hamming (from legacy fuzzy matching)
hamming_dist <- 0

# UP site matches
UP_site_matches <- dim(R2_seq)[1]

# Number cell ranger CBs
cell_bender_cb <- dim(CB)[1]

# Number of CBs matched
cb_matched <- length(number_recovered_CB)

# Number of reads with CB
cb_matched_reads <- dim(R2_cells_h2)[1]

# Number of unique SB_CB pairings
CB_SB_unique_pairings <- length(unique(UMIs_per_CB_SB$CB_SB))

# Proportion reads not in CB dictionary (thrown out, appear as NAs after left join)
proportion_cb_region_not_in_dictionary <- sum(is.na(R2_CB_region_minus_dups_corrected$GEX_CB))/length(R2_CB_region_minus_dups_corrected$GEX_CB)

# Assemble metrics
`slide-tags-10x-matcher-summary` <- data.frame(hamming_dist, UP_site_matches, cell_bender_cb, cb_matched, cb_matched_reads, CB_SB_unique_pairings, proportion_cb_region_not_in_dictionary)
print(`slide-tags-10x-matcher-summary`)

# Write table
write.table(`slide-tags-10x-matcher-summary`, matcher_summary_output_path, quote = FALSE, sep = "\t")