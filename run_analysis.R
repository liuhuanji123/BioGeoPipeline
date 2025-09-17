# run_analysis.R

# Load all the functions from your functions file
source("tree_biogeography_version_1.R")

# --- 1. Define related groups ---
sister_group_family <- c("Silphidae_135", "Leiodidae_185", "Agyrtidae_185")
distantly_group <- c("Scarabaeidae", "Tenebrionidae")

# --- 2. Define input file paths (must be relative to the project root) ---
allsequences_path <- "data/test_staphylinidae_barcodes.fasta"
mitogenomes_path <- "data/test_staphylinidae_mitogenomes.fasta"
timeperiods_filepath <- "data/my_timeperiods_20mya.txt"
dispersal_multipliers_filepath <- "data/my_dispersal_multipliers_20mya.txt"

# --- 3. Set computational resources ---
cpu_threads <- 8

# Pass the configuration variables as arguments to the main function
results <- tree_biogeography_pipeline(allsequences_path=allsequences_path,
                                      mitogenomes_path=mitogenomes_path,
                                      mito_length=11000,
                                      sister_group_family=sister_group_family,
                                      distantly_group=distantly_group,
                                      threads = cpu_threads,
                                      timeperiods_filepath = timeperiods_filepath,
                                      dispersal_multipliers_filepath = dispersal_multipliers_filepath)
