library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

#annotation (bed) file, if present, will need 4 columns, tab delimited and no header
#4 fields: chrom_name, start, end, feature_name
#this is used to check for overlap of the fragments with specific features/positions of interest.
#the name of this file is given as an optional 3rd argument in the bash wrapper script.

min_frag_size = 180
max_frag_size = 480

#run this for each possible combination of selectives
total.frags.in.range <- function(min.frag.size,
                             max.frag.size,
                             selectives,
                             rev = FALSE,
			     side1_re_name,
			     side2_re_name
                             ) {
  
  #filter the fragments for those with at least one side "2" adaptor and with both overhangs matching an adaptor
  if (rev == FALSE) {
    frags_filt <- frags_all %>% 
      filter(side1 == 2 | side2 == 2) %>%
      filter(seq1 %in% selectives & seq2 %in% selectives) %>%
      filter(length >= min.frag.size & length <= max.frag.size)
  
    frag_2k <- frags_filt %>%
      filter(length <= 2000)

    if (side1_re_name != side2_re_name) {
     pdf(file = paste0(args[1], "_out/", side1_re_name, "_", side2_re_name, "/", side1_re_name, "_", side2_re_name, "_", selectives[1], ":", selectives[2], "-frag2k.pdf"))
     }
    else {
     pdf(file = paste0(args[1], "_out/", side1_re_name, "/", side1_re_name, "_", side2_re_name, "_", selectives[1], "-frag2k.pdf"))
    }		

    hist(frag_2k$length, main = paste0(side1_re_name, "_", side2_re_name, " under 2k"))
    dev.off()
  }
  
  if (rev == TRUE) {
    frags_filt <- frags_all %>% 
      filter(side1 == 1 | side2 == 1) %>%
      filter(seq1 %in% selectives & seq2 %in% selectives) %>%
      filter(length >= min.frag.size & length <= max.frag.size)

    frag_2k <- frags_filt %>%
      filter(length <= 2000)
    
    if (side1_re_name != side2_re_name) {
     pdf(file = paste0(args[1], "_out/", side1_re_name, "_", side2_re_name, "/", side2_re_name, "_", side1_re_name, "_", selectives[2], ":", selectives[1], "-frag2k.pdf"))
     }
    else {
     pdf(file = paste0(args[1], "_out/", side1_re_name, "/", side1_re_name, "_", side1_re_name, "_", selectives[1], "-frag2k.pdf"))
    }
    hist(frag_2k$length, main = paste0(side2_re_name, "_", side1_re_name, " under 2k"))
    dev.off()
  }
  
  num_frags <- nrow(frags_filt)
  
  
  #return a data frame for the current set of selective adaptors that can be rbound to a master df
  if(rev == FALSE) {
    return(data.frame("genome" = args[1],
                "RE1" = side1_re,
                "RE2" = side2_re,
                "adaptor1" = selectives[1],
                "adaptor2" = selectives[2],
                "num_in_range" = num_frags))
  }
  
  if(rev == TRUE) {
    return(data.frame("genome" = args[1],
               "RE1" = side2_re,
               "RE2" = side1_re,
               "adaptor1" = selectives[2],
               "adaptor2" = selectives[1],
               "num_in_range" = num_frags))
  }
}


#bed_df has "chr", "start", "end", "feature_ID"
feature.overlaps <- function(selectives, 
                         features = args[4], 
                         rev = FALSE,
                         min.frag.size,
                         max.frag.size) {
  
  summary_by_feature <- data.frame("feature_ID" = character(),
                                   "feature_length" = integer(),
                                   "sites_covered" = integer(),
                                   "pct_covered" = double())
    
  #unnest bed positions
  bed <- read_tsv(features, col_names = FALSE, col_types=c("ciic"))
  names(bed) <- c("chr", "start", "end", "feature_ID")
  
  
  for (x in unique(bed$chr)) {
    #get vector of positions ('chr_frag_positions') in fragments recovered 
    #for the current chr
    chr_frags <- frags_all %>%
      filter(chr == x & length >= min.frag.size & length <= max.frag.size)
    
    if (rev == FALSE){
      chr_frags <- chr_frags %>% 
        filter(side1 == 2 | side2 == 2) %>%
        filter(seq1 %in% selectives & seq2 %in% selectives)
    }
    
    if (rev == TRUE){
      chr_frags <- chr_frags %>% 
        filter(side1 == 1 | side2 == 1) %>%
        filter(seq1 %in% selectives & seq2 %in% selectives)
    }
      
    chr_frag_positions <- mapply(FUN = function(start, end) {
      start:end },
      start = chr_frags$start_pos, 
      end = chr_frags$end_pos
    ) %>% unlist() %>% as.vector()
    
    #doing this to try to keep the RAM requirement down for bed files with lot of positions.
    
    unique_features <- unique(bed$feature_ID)
    nfeatures <- length(unique_features)
    nfeatures_per_loop <- ceiling(nfeatures/50)
    
    for (i in 0:48) {
      start_idx <- i*nfeatures_per_loop+1
      end_idx <- i*nfeatures_per_loop+nfeatures_per_loop
      target_features <- unique_features[start_idx:end_idx]
        
      bed_sub <- bed %>% filter(feature_ID %in% target_features)
      #bed_sub <- bed %>% dplyr::slice((i*nrows_per_loop)+1:(i*nrows_per_loop)+nrows_per_loop)
      feat_pos_df <- bed_sub %>% 
        filter(chr == x)
        
       if (nrow(feat_pos_df) > 0) {
        feat_pos_df <- feat_pos_df %>%
        group_by(r = row_number()) %>% 
        mutate("feat_pos" = list(start:end)) %>% 
        ungroup %>% 
        unnest(feat_pos) %>%
        select(feature_ID, feat_pos) %>%
        mutate("covered" = ifelse(feat_pos %in% chr_frag_positions, 1, 0))
      
      tmp_summary <- feat_pos_df %>%
        group_by(feature_ID) %>%
        summarise("feature_length" = n(),
                  "sites_covered" = sum(covered),
                  "pct_covered" = sites_covered/feature_length)
      
      summary_by_feature <- rbind(summary_by_feature, tmp_summary)
     }  
    }
    
    start_idx <- 49*nfeatures_per_loop+1
    #end_idx <- 49*nfeatures_per_loop+nfeatures_per_loop
    target_features <- unique_features[start_idx:length(unique_features)]
    
    bed_sub <- bed %>% filter(feature_ID %in% target_features)
    
    feat_pos_df <- bed_sub %>% 
      filter(chr == x)
    
    if (nrow(feat_pos_df) > 0) {
      feat_pos_df <- feat_pos_df %>%
      group_by(r = row_number()) %>% 
      mutate("feat_pos" = list(start:end)) %>% 
      ungroup %>% 
      unnest(feat_pos) %>%
      select(feature_ID, feat_pos) %>%
      mutate("covered" = ifelse(feat_pos %in% chr_frag_positions, 1, 0))
    
    tmp_summary <- feat_pos_df %>%
      group_by(feature_ID) %>%
      summarise("feature_length" = n(),
                "sites_covered" = sum(covered),
                "pct_covered" = sites_covered/feature_length)
    
    summary_by_feature <- rbind(summary_by_feature, tmp_summary)
   } 
  }
  
  summary_by_feature
}


###############################################################

#assume there is no features file
#if there is, 'bed' is assigned the name of the file
bed <- NULL
if (!is.na(args[4])) {
  bed <- args[4]
}

#the full set of cut sites identified from the bash script are stored as 'cuts'
cuts <- read_tsv(paste0(args[1], "_out/", args[2], "/cutsites.tsv"), col_names=FALSE, col_types=c("icc"))
names(cuts) <- c("pos", "seq", "chr")

meta <- read_tsv(args[3], col_names = FALSE)
names(meta) <- c("seq", 
                 "RE_name", 
                 "side")

#add in RE_name and RE_side (arbitrary) to cuts data
merged_all <- left_join(cuts, meta, by = "seq")

get_pairs1 <- function(x) {
  c(x[1:(length(x) - 1)], NA)
}

get_pairs2 <- function(x) {
  c(x[2:length(x)], NA)
}

#check to make sure all chromosomes have at least 2 cut sites. If not, script will break.
dup_chrs <- merged_all$chr[duplicated(merged_all$chr)] %>% unique()
merged_all <- merged_all %>% filter(chr %in% dup_chrs)

frags_all <- merged_all %>% 
  group_by(chr) %>% 
  mutate("start_pos" = get_pairs1(x = pos)) %>%
  mutate("end_pos" = get_pairs2(x = pos)) %>%
  mutate("length" = end_pos - start_pos) %>%
  mutate("RE1" = get_pairs1(x = RE_name)) %>%
  mutate("RE2" = get_pairs2(x = RE_name)) %>%
  mutate("side1" = get_pairs1(x = side)) %>%
  mutate("side2" = get_pairs2(x = side)) %>%
  mutate("seq1" = get_pairs1(x = seq)) %>%
  mutate("seq2" = get_pairs2(x = seq)) %>%
  drop_na() %>% ungroup()

#get the names of the RE's, and order them based on the side they represent
#this isn't used after the next ordering step
re_names <- frags_all %>% 
  select(RE1,side1) %>% 
  distinct() %>% 
  unite(col = "RE_side", sep = "_") %>%
  unlist() %>% as.character()

#Order the re_names - this stores just the names of the RE's
side1_re <- re_names[grep("_1$", re_names)]
side1_re <- gsub("_1", "", side1_re)
side2_re <- re_names[grep("_2$", re_names)]
side2_re <- gsub("_2", "", side2_re)

#now start filtering
#for fragments with at least one side 2 (think that's required for amplification)
#for fragments for which seq1 and seq2 will both match the adaptor overhang used 
#evaluate for each pair of possible adaptor overhangs (only a single pair if no ambguities in RE recognition seq)

#get set of potential overhangs for RE1
RE1_selectives <- meta %>% 
  filter(RE_name == side1_re) %>%
  #filter(RE_name == gsub("(.+)(_\\d)", "\\1", re_names[1])) %>%
  select(seq) %>%
  distinct() %>% 
  unlist() %>% 
  as.character()

#get set of potential overhangs for RE2
RE2_selectives <- meta %>%
  filter(RE_name == side2_re) %>%
  select(seq) %>%
  distinct() %>% 
  unlist() %>% 
  as.character()

master_df <- data.frame("genome" = character(), "RE1"= character(), "RE2" = character(), "adaptor1" = character(), "adaptor2" = character(), "num_in_range" = integer())

#iterate over all possible combinations of adaptor overhangs and consider
#each RE on each side
#for each, bind the results to master_df

for (i in 1:length(RE1_selectives)) {
  for (j in 1:length(RE2_selectives)) {
    
    #get the current pair of overhangs
    selectives <- c(RE1_selectives[i], RE2_selectives[j])
    
    current_selectives_df <- total.frags.in.range(selectives = selectives,
                                                  min.frag.size = min_frag_size,
                                                  max.frag.size = max_frag_size, 
                                                  side1_re_name = side1_re, 
                                                  side2_re_name = side2_re)

    current_selectives_df_rev <- total.frags.in.range(selectives = selectives, 
                                                      min.frag.size = min_frag_size,
                                                      max.frag.size = max_frag_size,
                                                      rev = TRUE, 
                                                      side1_re_name = side1_re, 
                                                      side2_re_name = side2_re)
    
    master_df <- rbind(master_df, current_selectives_df, current_selectives_df_rev)
  }
}  


#if there's a bed file, check feature overlap

if (!is.null(bed)) {
  
  #iterate over all possible combinations of adaptor overhangs and consider
  #each RE on each side, just like above for master_df - will give same # of rows for cbind
  
  overall_summary_master <- data.frame()
  summary_by_feature <- data.frame("feature_ID" = character(),
                                   "feature_length" = integer(),
                                   "sites_covered" = integer(),
                                   "pct_covered" = double(),
                                   "RE1" = character(),
                                   "RE2" = character(),
                                   "overhang1" = character(),
                                   "overhang2" = character())
  
  for (i in 1:length(RE1_selectives)) {
    for (j in 1:length(RE2_selectives)) {
  
      #get the current pair of overhangs
      selectives <- c(RE1_selectives[i], RE2_selectives[j])
      
      
      #first get info on the individual features
      
      coverage <- feature.overlaps(selectives = selectives,
                                  min.frag.size = min_frag_size,
                                  max.frag.size = max_frag_size)
      
      coverage_rev <- feature.overlaps(selectives = selectives,
                                      min.frag.size = min_frag_size,
                                      max.frag.size = max_frag_size,
                                      rev = TRUE)
      
      summary <- as.data.frame(coverage)
      summary_rev <- as.data.frame(coverage_rev)
      
      summary$RE1 <- side1_re
      summary$RE2 <- side2_re
      summary$overhang1 <- RE1_selectives[i]
      summary$overhang2 <- RE2_selectives[j]
      
      summary_rev$RE1 <- side2_re
      summary_rev$RE2 <- side1_re
      summary_rev$overhang1 <- RE2_selectives[i]
      summary_rev$overhang2 <- RE1_selectives[j]
      
      summary_by_feature <- rbind(summary_by_feature, summary, summary_rev)
      
      
      #next get the overall summary to col bind to master_df
      #need to get # features evaluated, # partially covered, # fully covered, total # covered, % features covered
  
      feats_eval <- nrow(coverage)
      partial_feat_cov <- coverage %>% filter(pct_covered > 0 & pct_covered < 100) %>% nrow()
      whole_feat_cov <- coverage %>% filter(pct_covered == 100) %>% nrow()
      tot_feat_cov <- coverage %>% filter(pct_covered > 0) %>% nrow()
      pct_cov <- tot_feat_cov/feats_eval
      
      overall_summary <- data.frame("features_evaluated" = feats_eval, 
                                    "features_covered_partial" = partial_feat_cov, 
                                    "features_covered_whole" = whole_feat_cov,
                                    "total_features_covered" = tot_feat_cov,
                                    "pct_features_covered" = pct_cov)
      
      feats_eval <- nrow(coverage_rev)
      partial_feat_cov <- coverage_rev %>% filter(pct_covered > 0 & pct_covered < 100) %>% nrow()
      whole_feat_cov <- coverage_rev %>% filter(pct_covered == 100) %>% nrow()
      tot_feat_cov <- coverage_rev %>% filter(pct_covered > 0) %>% nrow()
      pct_cov <- tot_feat_cov/feats_eval
      
      overall_summary_rev <- data.frame("features_evaluated" = feats_eval, 
                                        "features_covered_partial" = partial_feat_cov, 
                                        "features_covered_whole" = whole_feat_cov,
                                        "total_features_covered" = tot_feat_cov,
                                        "pct_features_covered" = pct_cov)
      
      overall_summary_master <- rbind(overall_summary_master, overall_summary, overall_summary_rev)

    }
  }
  
  write.table(summary_by_feature,
              file = paste0(args[1], "_out/", side1_re, "_", side2_re, "_features_detail.txt"),
              sep = "\t", 
              quote = FALSE, 
              col.names = TRUE, 
              row.names = FALSE)
  
  master_df <- cbind(master_df, overall_summary_master)
  write.table(master_df,
              file = paste0(args[1], "_out/", side1_re, "_", side2_re, "_summary.txt"),
              sep = "\t", 
              quote = FALSE, 
              col.names = TRUE, 
              row.names = FALSE)
  
}    
      
  
if (is.null(bed)) {  
 write.table(master_df,
             file = paste0(args[1], "_out/", side1_re, "_", side2_re, "_summary.txt"),
             sep = "\t", 
             quote = FALSE, 
             col.names = TRUE, 
             row.names = FALSE)
}  
