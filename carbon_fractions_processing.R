### R script to calculate carbon fractions leaf-level CABO data from Fulcrum
# and to compare calculated C fractions in Fulcrum to those here
# Etienne Laliberté
# April 7, 2020


### Step 1: Export all C fractions data from Fulcrum ----
# See "carbon_fractions" folder for an example exported on April 7
# I exported the app "Carbon fractions" as a .csv and selected all projects
# But you can export for a single project, or just a few. Script will also work.


### Step 2: Load libraries ----
library(tidyverse)


### Step 3: Read the data ----
batches <- read_csv('carbon_fractions/carbon_fractions.csv') # parent records: each record is a batch (i.e. a run)
bags <- read_csv('carbon_fractions/carbon_fractions_bags.csv') # child records: each record is a bag in a batch


#### Step 4: Remove CABO-test batches ----
# Find the fulcrum_id of good batches (not CABO-test ones, not deleted ones)
good_batches <- batches %>% 
  filter(project != 'CABO-test',  status != 'deleted') %>% 
  select(fulcrum_id) %>% 
  distinct() %>% 
  pull()

# Filter the batches and bags
batches <- filter(batches, fulcrum_id %in% good_batches)
bags <- filter(bags, fulcrum_parent_id %in% good_batches) # in bags, fulcrum_parent_id = id of batch


### Step 5: Check the data ----
# How many good batches?
nrow(batches) # 90 in this export dated 2020-04-07

# Find number of blanks per batch
batches_n_blanks <- bags %>% 
  group_by(fulcrum_parent_id, sample_type) %>% 
  count() %>% 
  filter(sample_type == 'blank') # All batches have exactly ONE blank

# Same but for samples
batches_n_samples <- bags %>% 
  group_by(fulcrum_parent_id, sample_type) %>% 
  count() %>% 
  filter(sample_type == 'sample') # Max is 23 (+1 blank = 24 bags)

# How many batches have ONE blank?
batches_n_blanks %>% 
  filter(n == 1) %>% 
  nrow() # 90, all of them

# How many batches have at least one sample
batches_n_samples %>% 
  filter(n >= 1) %>% 
  nrow() # 90, all of them


### Step 6: Write a function to calculate C fractions for one batch ----
calc_c_frac <- function(x) {
  
  # Get all blanks correction factors
  corr_facts <- x %>% 
    filter(sample_type == 'blank') %>% # get blanks
    mutate(ndf_correction_factor_R =  ndf_weight_g / empty_bag_weight_g,
           adf_correction_factor_R = adf_weight_g / empty_bag_weight_g,
           adl_correction_factor_R = adl_weight_g / empty_bag_weight_g,
           ash_weight_g_R = crucible_ash_weight_g - empty_crucible_weight_g,
           ashing_correction_factor_R = ash_weight_g_R / empty_bag_weight_g)
    
    # Get average correction factors per batch
    corr_facts_avg <- corr_facts %>% 
      select(ndf_correction_factor_R,
             adf_correction_factor_R,
             adl_correction_factor_R,
             ashing_correction_factor_R) %>% 
      summarise_all(mean, na.rm = T) # for cases where there would be multiple blanks per batch
       
  # Calculate NDF, ADF, ADL, ash,
  # and soluble C, hemicellulose, cellulose, lignin, recalcitrants, sum of fractions
  calc_fracs <- x %>% 
    filter(sample_type == 'sample') %>% 
    mutate(ndf_perc_R = (ndf_weight_g - (empty_bag_weight_g * corr_facts_avg$ndf_correction_factor_R)) * 100 / sample_weight_g,
              adf_perc_R = (adf_weight_g - (empty_bag_weight_g * corr_facts_avg$adf_correction_factor_R)) * 100 / sample_weight_g,
              adl_perc_R = (adl_weight_g - (empty_bag_weight_g * corr_facts_avg$adl_correction_factor_R)) * 100 / sample_weight_g,
              ash_weight_g_R = crucible_ash_weight_g - empty_crucible_weight_g,
              soluble_perc_R = 100 - ndf_perc_R,
              hemicellulose_perc_R = ndf_perc_R - adf_perc_R,
              cellulose_perc_R = adf_perc_R - adl_perc_R,
              recalcitrants_perc_R = (ash_weight_g_R - (empty_bag_weight_g * corr_facts_avg$ashing_correction_factor_R)) * 100 / sample_weight_g,
              recalcitrants_perc_R = replace(recalcitrants_perc_R, recalcitrants_perc_R < 0, 0), # replace negative recalcitrants values by 0
              lignin_perc_R = adl_perc_R - recalcitrants_perc_R,
              sum_fractions_perc_R = soluble_perc_R + hemicellulose_perc_R + cellulose_perc_R + lignin_perc_R + recalcitrants_perc_R)
  
  
  # Bind blanks and samples together
  fracs_bind <- bind_rows(corr_facts, calc_fracs) %>% 
    
    # Select only relevant variables and order
    select(bag_id = fulcrum_id,
           fulcrum_parent_id,
           bag_number,
           sample_type,
           leaf_chemistry_sample,
           bottle_id,
           sample_remarks,
           empty_bag_weight_g,
           sample_weight_g,
           total_initial_weight_g,
           ndf_weight_g,
           ndf_correction_factor,
           ndf_correction_factor_R,
           adf_weight_g,
           adf_correction_factor,
           adf_correction_factor_R,
           adl_weight_g,
           adl_correction_factor,
           adl_correction_factor_R,
           empty_crucible_weight_g,
           crucible_adl_weight_g,
           crucible_ash_weight_g,
           ash_weight_g,
           ash_weight_g_R,
           ashing_correction_factor,
           ashing_correction_factor_R,
           soluble_perc,
           soluble_perc_R,
           hemicellulose_perc,
           hemicellulose_perc_R,
           cellulose_perc,
           cellulose_perc_R,
           lignin_perc,
           lignin_perc_R,
           recalcitrants_perc,
           recalcitrants_perc_R,
           sum_fractions_perc,
           sum_fractions_perc_R) %>% 
    arrange(bag_number) # order by bag number
}


### Step 7: Calculate C fractions per sample for each batch ----
c_frac_all <- bags %>% 
  group_by(fulcrum_parent_id) %>% # Group bags by batch
  do(calc_c_frac(.)) # Runs the function above on each batch
  

#### Step 8: Merge with batch data
# Select only relevant batch variables
batches_relev <- batches %>% 
  select(batch_id = fulcrum_id,
         status,
         project,
         analysis_id,
         measured_by,
         date_started,
         analysis_remarks,
         ndf_correction_factors:ashing_correction_factors,
         verified_by,
         date_verified
         )

# Merge the two data sets
bags_all <- batches_relev %>% 
  left_join(c_frac_all, by = c('batch_id' = 'fulcrum_parent_id')) %>% 
  arrange(date_started,
          bag_number)


### Step 9: Save as .csv file ----
dir.create('output')
write_csv(bags_all, 'output/carbon_fractions_calculated.csv')
