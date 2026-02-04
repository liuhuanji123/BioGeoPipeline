#20260204
# The pipeline consists of the following steps:
# construction_tree: Constructs the phylogenetic tree.
# rooting_tree: Roots the tree generated in the previous step.
# dating_tree: Dates the backbone tree. While it's possible to date the reconstruction tree directly,
# for the reconstruction tree, we opt to use the dating results from the backbone via the fix_backbone_dating_to_reconstruction function.
# biogeography_analysis_pastml: Performs biogeographical analysis using PastML.
# biogeography_analysis_biogeobears: Performs biogeographical analysis using BioGeoBEARS.

# =============================================================================
#  ENVIRONMENT PRE-CHECK (runs automatically on source)
# =============================================================================
#
#  This function checks all dependencies before running the pipeline.
#  It will report all missing components at once, so you can fix them all
#  before starting a long analysis.
#
# =============================================================================

check_pipeline_environment <- function(verbose = TRUE) {

  cat("\n")
  cat("###########################################################################\n")
  cat("#              PIPELINE ENVIRONMENT PRE-CHECK                            #\n")
  cat("###########################################################################\n\n")

  errors <- c()
  warnings <- c()

  # --- 1. Check R packages ---
  cat("--- 1. Checking R packages ---\n")

  required_packages <- c(
    "ape", "phytools", "stringr", "dplyr", "tidyr", "ggplot2",
    "BioGeoBEARS", "cladoRcpp", "pheatmap", "reshape2", "gridExtra", "grid"
  )

  optional_packages <- c(
    "ggtree", "treeio", "tidyverse", "circlize", "MultinomialCI"
  )

  missing_required <- c()
  missing_optional <- c()

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_required <- c(missing_required, pkg)
    }
  }

  for (pkg in optional_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_optional <- c(missing_optional, pkg)
    }
  }

  if (length(missing_required) > 0) {
    errors <- c(errors, paste0("Missing REQUIRED R packages: ", paste(missing_required, collapse = ", ")))
    cat("  [FAIL] Missing required packages:", paste(missing_required, collapse = ", "), "\n")
  } else {
    cat("  [OK] All required R packages installed\n")
  }

  if (length(missing_optional) > 0) {
    warnings <- c(warnings, paste0("Missing optional R packages: ", paste(missing_optional, collapse = ", ")))
    cat("  [WARN] Missing optional packages:", paste(missing_optional, collapse = ", "), "\n")
  }

  # --- 2. Check WSL ---
  cat("\n--- 2. Checking WSL ---\n")

  wsl_check <- suppressWarnings(system("wsl echo OK", intern = TRUE, ignore.stderr = TRUE))
  if (length(wsl_check) == 0 || !any(grepl("OK", wsl_check))) {
    errors <- c(errors, "WSL is not available or not configured properly")
    cat("  [FAIL] WSL not available\n")
  } else {
    cat("  [OK] WSL is available\n")
  }

  # --- 3. Check PastML in WSL ---
  cat("\n--- 3. Checking PastML (in WSL) ---\n")

  pastml_path <- ""
  tryCatch({
    con <- pipe('wsl bash -ic "which pastml"', "r")
    lines <- readLines(con, warn = FALSE)
    close(con)
    for (line in lines) {
      line <- trimws(line)
      if (grepl("^/.*pastml$", line)) {
        pastml_path <- line
        break
      }
    }
  }, error = function(e) {})

  if (nchar(pastml_path) == 0) {
    errors <- c(errors, "PastML not found in WSL. Install with: wsl pip install pastml")
    cat("  [FAIL] PastML not found\n")
  } else {
    cat("  [OK] PastML found:", pastml_path, "\n")
  }

  # --- 4. Check RAxML-NG in WSL ---
  cat("\n--- 4. Checking RAxML-NG (in WSL) ---\n")

  raxml_path <- ""
  tryCatch({
    con <- pipe('wsl bash -ic "which raxml-ng"', "r")
    lines <- readLines(con, warn = FALSE)
    close(con)
    for (line in lines) {
      line <- trimws(line)
      if (grepl("raxml-ng", line)) {
        raxml_path <- line
        break
      }
    }
  }, error = function(e) {})

  if (nchar(raxml_path) == 0) {
    warnings <- c(warnings, "RAxML-NG not found in WSL (only needed for tree construction)")
    cat("  [WARN] RAxML-NG not found (only needed for tree construction)\n")
  } else {
    cat("  [OK] RAxML-NG found:", raxml_path, "\n")
  }

  # --- 5. Check treePL in WSL ---
  cat("\n--- 5. Checking treePL (in WSL) ---\n")

  treepl_path <- ""
  tryCatch({
    con <- pipe('wsl bash -ic "which treePL"', "r")
    lines <- readLines(con, warn = FALSE)
    close(con)
    for (line in lines) {
      line <- trimws(line)
      if (grepl("treePL", line)) {
        treepl_path <- line
        break
      }
    }
  }, error = function(e) {})

  if (nchar(treepl_path) == 0) {
    warnings <- c(warnings, "treePL not found in WSL (only needed for tree dating)")
    cat("  [WARN] treePL not found (only needed for tree dating)\n")
  } else {
    cat("  [OK] treePL found:", treepl_path, "\n")
  }

  # --- Summary ---
  cat("\n")
  cat("###########################################################################\n")
  cat("#                         PRE-CHECK SUMMARY                              #\n")
  cat("###########################################################################\n\n")

  if (length(errors) == 0 && length(warnings) == 0) {
    cat("  [ALL OK] Environment is ready for full pipeline execution!\n\n")
    return(invisible(TRUE))
  }

  if (length(warnings) > 0) {
    cat("  WARNINGS (pipeline may still work for some functions):\n")
    for (w in warnings) {
      cat("    -", w, "\n")
    }
    cat("\n")
  }

  if (length(errors) > 0) {
    cat("  ERRORS (must fix before running pipeline):\n")
    for (e in errors) {
      cat("    -", e, "\n")
    }
    cat("\n")
    cat("  To install missing R packages, run:\n")
    cat("    install.packages(c(\"ape\", \"phytools\", \"stringr\", \"dplyr\", \"tidyr\", \"ggplot2\", \"pheatmap\", \"reshape2\", \"gridExtra\"))\n")
    cat("    # For BioGeoBEARS (from GitHub):\n")
    cat("    devtools::install_github(\"nmatzke/BioGeoBEARS\")\n")
    cat("\n")
    cat("  To install PastML in WSL, run:\n")
    cat("    wsl pip install pastml\n")
    cat("\n")
    stop("Environment check failed. Please fix the errors above before running the pipeline.")
  }

  return(invisible(TRUE))
}

# Run environment check automatically on source
# Set CHECK_ON_SOURCE <- FALSE before source() to skip
if (!exists("CHECK_ON_SOURCE") || CHECK_ON_SOURCE != FALSE) {
  check_pipeline_environment()
}

# =============================================================================
#  PIPELINE FUNCTIONS START HERE
# =============================================================================

# Step 1: Construct a backbone tree (mitogenomes only) and a reconstruction tree (mitogenomes + barcode).
# This step uses multiple distant groups (for rooting), one sister group (for dating), and sequences of the target family.
# The format for sequences and tree tip labels is: genetype_id_taxon_biogeographic_realm, e.g., MMG_GBDL00960_Cleridae_NA, MMG_BIOD06517_Cleridae_NT, GMT_GMNGF130-16_Cleridae_OC.
# MMG represents mitogenomes, and GMT represents barcode sequences.

# The first step is to build the backbone and reconstruction trees.
# First, obtain an initial tree from a preliminary search. This initial tree is sourced from a lineage on a large tree database,
# where the lineage predominantly belongs to a single family.
# nwk_file<-"C:\\Users\\16575\\Documents\\Rdocuments\\test_20250703\\Cleridae_test_20250725/Cleridae_1.nwk"
convert_path_to_wsl <- function(win_path) {
  wsl_path <- gsub("\\\\", "/", win_path)
  wsl_path <- sub("([A-Za-z]):", "/mnt/\\L\\1", wsl_path, perl = TRUE)
  return(wsl_path)
}

extract_family <- function(tip) {
  parts <- str_split_fixed(tip, "_", 4)
  parts[,3]  # The second-to-last element is the taxon name.
}

#' @title Sanitize a string to be a valid filename
#' @description Replaces potentially problematic characters in a string with underscores.
#' This is useful for creating safe filenames from variable names or descriptions.
#'
#' @param filename A character string to be sanitized.
#' @param replacement_char The character to use for replacing special characters.
#' Defaults to "_".
#'
#' @return A sanitized character string, safe for use as a filename.
#'
sanitize_filename <- function(filename, replacement_char = "_") {
  # Define a regular expression that matches common problematic characters.
  # This includes: space, +, /, \, :, *, ?, ", <, >, |
  # The \\ is an escaped backslash.
  bad_chars_regex <- "[+ /\\\\:*?\"<>|]"
  
  # Use gsub() to find all occurrences of the bad characters and replace them.
  sanitized_name <- gsub(pattern = bad_chars_regex, 
                         replacement = replacement_char, 
                         x = filename)
  
  return(sanitized_name)
}

# =============================================================================
#
#  BIOGEOGRAPHY PRE-PROCESSING FUNCTIONS
#
# =============================================================================
#
# These functions prepare inputs for PastML and BioGeoBEARS analyses.
# They are designed to be called independently, allowing step-by-step execution.
#
# Workflow:
#   Step A: prune_tree_to_target_family()  -> Remove non-target family tips
#   Step B: prepare_geography_for_pastml() -> Create PastML-compatible CSV
#   Step C: prepare_geography_for_biogeobears() -> Create BioGeoBEARS-compatible CSV + tree
#   Step D: pastml_process_tree()          -> Run PastML analysis
#   Step E: run_biogeobears_pipeline()     -> Run BioGeoBEARS analysis
#
# =============================================================================


# =============================================================================
# Step A: prune_tree_to_target_family
# =============================================================================
#
# Purpose:
#   Remove tips that do not belong to the target family.
#   This is typically used to remove outgroups/sister groups used for dating/rooting.
#
# Input:
#   - tree_path: Path to the input tree file (.nwk or .tre)
#   - target_family: Name of the target family (e.g., "Cleridae")
#                    If NULL, automatically selects the most frequent family.
#   - output_path: Path for the output tree. If NULL, adds "_pruned" suffix.
#
# Output:
#   - Writes the pruned tree to output_path
#   - Returns the output path (invisibly)
#
# Note:
#   This function assumes tip labels follow the format: genetype_id_taxon_realm
#   where taxon (parts[3]) is the family name.
#
# =============================================================================

prune_tree_to_target_family <- function(tree_path,
                                         target_family = NULL,
                                         output_path = NULL) {
  library(ape)
  library(stringr)

  cat("\n")
  cat("###########################################################################\n")
  cat("#  Step A: Prune Tree to Target Family                                   #\n")
  cat("###########################################################################\n\n")

  # Validate input

if (!file.exists(tree_path)) {
    stop("ERROR: Tree file not found: ", tree_path)
  }

  # Read tree
  tree <- ape::read.tree(tree_path)
  original_n_tips <- length(tree$tip.label)
  cat("Input tree:", tree_path, "\n")
  cat("Original number of tips:", original_n_tips, "\n")

  # Extract family from tip labels
  all_tips <- tree$tip.label
  families <- sapply(all_tips, function(tip) {
    parts <- str_split_fixed(tip, "_", 4)
    parts[, 3]
  }, USE.NAMES = FALSE)
  names(families) <- all_tips

  # Determine target family
  if (is.null(target_family)) {
    family_table <- table(families[families != "" & !is.na(families)])
    if (length(family_table) == 0) {
      stop("ERROR: Could not extract any family names from tip labels.")
    }
    target_family <- names(which.max(family_table))
    cat("Target family (auto-detected):", target_family, "\n")
  } else {
    cat("Target family (user-specified):", target_family, "\n")
  }

  # Identify tips to remove
  tips_to_drop <- names(families)[families != target_family | is.na(families) | families == ""]
  tips_to_drop <- intersect(tips_to_drop, all_tips)

  if (length(tips_to_drop) == 0) {
    cat("No tips to remove. Tree unchanged.\n")
  } else {
    cat("Removing", length(tips_to_drop), "tips from other families...\n")
    tree <- ape::drop.tip(tree, tips_to_drop)
    cat("Remaining tips:", length(tree$tip.label), "\n")
  }

  # Determine output path
  if (is.null(output_path)) {
    output_path <- sub("\\.(nwk|tre|tree)$", "_pruned.\\1", tree_path)
    if (output_path == tree_path) {
      output_path <- paste0(tools::file_path_sans_ext(tree_path), "_pruned.nwk")
    }
  }

  # Write output
  ape::write.tree(tree, file = output_path)
  cat("Output tree saved to:", output_path, "\n")
  cat("\n")

  return(invisible(output_path))
}


# =============================================================================
# Step B: prepare_geography_for_pastml
# =============================================================================
#
# Purpose:
#   Prepare a geography CSV file compatible with PastML requirements.
#
# PastML requirements:
#   - Each tip must have exactly ONE state (realm)
#   - Tips without realm data are allowed (will have empty state)
#   - All tips in the tree must be present in the CSV
#
# Input:
#   - tree_path: Path to the tree file
#   - geography_csv: Path to the input geography CSV (columns: tip, realm)
#                    Multiple realms can be separated by ";" or ","
#   - output_csv: Path for output CSV. If NULL, auto-generated.
#
# Processing:
#   - For tips with multiple realms: keep only the FIRST one
#   - For tips not in CSV: add them with empty realm
#   - All tree tips are preserved
#
# Output:
#   - Writes PastML-compatible CSV (columns: ID, realm)
#   - Returns the output path (invisibly)
#
# =============================================================================

prepare_geography_for_pastml <- function(tree_path,
                                          geography_csv,
                                          output_csv = NULL) {
  library(ape)

  cat("\n")
  cat("###########################################################################\n")
  cat("#  Step B: Prepare Geography for PastML                                  #\n")
  cat("###########################################################################\n\n")

  # Validate inputs
  if (!file.exists(tree_path)) {
    stop("ERROR: Tree file not found: ", tree_path)
  }
  if (!file.exists(geography_csv)) {
    stop("ERROR: Geography CSV not found: ", geography_csv)
  }

  # Read tree
  tree <- ape::read.tree(tree_path)
  tree_tips <- tree$tip.label
  cat("Tree:", tree_path, "\n")
  cat("Number of tips in tree:", length(tree_tips), "\n")

  # Read geography CSV
  geog_raw <- utils::read.csv(geography_csv, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(geog_raw) <- tolower(trimws(colnames(geog_raw)))

  # Identify columns
  if ("tip" %in% colnames(geog_raw)) {
    id_col <- "tip"
  } else if ("id" %in% colnames(geog_raw)) {
    id_col <- "id"
  } else {
    id_col <- colnames(geog_raw)[1]
  }

  if ("realm" %in% colnames(geog_raw)) {
    realm_col <- "realm"
  } else if ("state" %in% colnames(geog_raw)) {
    realm_col <- "state"
  } else {
    realm_col <- colnames(geog_raw)[2]
  }

  cat("Geography CSV:", geography_csv, "\n")
  cat("ID column:", id_col, ", Realm column:", realm_col, "\n")
  cat("Rows in CSV:", nrow(geog_raw), "\n")

  # Process: for each tip, get the FIRST realm only
  geog_processed <- data.frame(ID = character(), realm = character(), stringsAsFactors = FALSE)

  for (tip in tree_tips) {
    matching_rows <- geog_raw[[id_col]] == tip
    if (any(matching_rows)) {
      realm_str <- geog_raw[[realm_col]][which(matching_rows)[1]]
      # Split by ; or , and take only the first
      realms <- trimws(unlist(strsplit(as.character(realm_str), "[;,]")))
      realms <- realms[realms != "" & !is.na(realms) & tolower(realms) != "unk"]
      first_realm <- if (length(realms) > 0) realms[1] else ""
    } else {
      first_realm <- ""
    }
    geog_processed <- rbind(geog_processed, data.frame(ID = tip, realm = first_realm, stringsAsFactors = FALSE))
  }

  # Statistics
  n_with_realm <- sum(geog_processed$realm != "")
  n_without_realm <- sum(geog_processed$realm == "")
  cat("\nProcessing complete:\n")
  cat("  Tips with realm:", n_with_realm, "\n")
  cat("  Tips without realm:", n_without_realm, "(will have empty state in PastML)\n")

  # Determine output path
  if (is.null(output_csv)) {
    output_csv <- file.path(dirname(geography_csv),
                            paste0(tools::file_path_sans_ext(basename(geography_csv)), "_pastml.csv"))
  }

  # Write output
  utils::write.csv(geog_processed, file = output_csv, row.names = FALSE, quote = FALSE)
  cat("Output CSV saved to:", output_csv, "\n\n")

  return(invisible(output_csv))
}


# =============================================================================
# Step C: prepare_geography_for_biogeobears
# =============================================================================
#
# Purpose:
#   Prepare geography CSV and tree compatible with BioGeoBEARS requirements.
#
# BioGeoBEARS requirements:
#   - Each tip MUST have at least one realm (tips without realm are removed)
#   - A tip CAN have multiple realms (each realm on a separate row in internal format)
#   - Tree must only contain tips that have realm data
#
# Input:
#   - tree_path: Path to the tree file
#   - geography_csv: Path to input geography CSV (columns: tip, realm)
#                    Multiple realms can be separated by ";" or ","
#   - output_csv: Path for output CSV. If NULL, auto-generated.
#   - output_tree: Path for output tree. If NULL, auto-generated.
#
# Processing:
#   - Remove tips from tree that have no realm data
#   - Expand multi-realm entries (each realm on separate row)
#
# Output:
#   - Writes BioGeoBEARS-compatible CSV (columns: ID, realm - one realm per row)
#   - Writes pruned tree (only tips with realm data)
#   - Returns list(csv = output_csv, tree = output_tree)
#
# =============================================================================

prepare_geography_for_biogeobears <- function(tree_path,
                                               geography_csv,
                                               output_csv = NULL,
                                               output_tree = NULL) {
  library(ape)

  cat("\n")
  cat("###########################################################################\n")
  cat("#  Step C: Prepare Geography for BioGeoBEARS                             #\n")
  cat("###########################################################################\n\n")

  # Validate inputs
  if (!file.exists(tree_path)) {
    stop("ERROR: Tree file not found: ", tree_path)
  }
  if (!file.exists(geography_csv)) {
    stop("ERROR: Geography CSV not found: ", geography_csv)
  }

  # Read tree
  tree <- ape::read.tree(tree_path)
  tree_tips <- tree$tip.label
  original_n_tips <- length(tree_tips)
  cat("Tree:", tree_path, "\n")
  cat("Number of tips in tree:", original_n_tips, "\n")

  # Read geography CSV
  geog_raw <- utils::read.csv(geography_csv, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(geog_raw) <- tolower(trimws(colnames(geog_raw)))

  # Identify columns
  if ("tip" %in% colnames(geog_raw)) {
    id_col <- "tip"
  } else if ("id" %in% colnames(geog_raw)) {
    id_col <- "id"
  } else {
    id_col <- colnames(geog_raw)[1]
  }

  if ("realm" %in% colnames(geog_raw)) {
    realm_col <- "realm"
  } else if ("state" %in% colnames(geog_raw)) {
    realm_col <- "state"
  } else {
    realm_col <- colnames(geog_raw)[2]
  }

  cat("Geography CSV:", geography_csv, "\n")
  cat("ID column:", id_col, ", Realm column:", realm_col, "\n")

  # Process: expand multi-realm entries, filter to tree tips
  geog_expanded <- data.frame(ID = character(), realm = character(), stringsAsFactors = FALSE)
  tips_with_realm <- character()

  for (i in seq_len(nrow(geog_raw))) {
    tip_id <- geog_raw[[id_col]][i]

    # Skip if tip not in tree
    if (!(tip_id %in% tree_tips)) next

    realm_str <- geog_raw[[realm_col]][i]
    # Split by ; or ,
    realms <- trimws(unlist(strsplit(as.character(realm_str), "[;,]")))
    realms <- realms[realms != "" & !is.na(realms) & tolower(realms) != "unk"]

    if (length(realms) > 0) {
      tips_with_realm <- c(tips_with_realm, tip_id)
      for (r in realms) {
        geog_expanded <- rbind(geog_expanded, data.frame(ID = tip_id, realm = r, stringsAsFactors = FALSE))
      }
    }
  }

  tips_with_realm <- unique(tips_with_realm)

  # Identify tips to remove (no realm data)
  tips_without_realm <- setdiff(tree_tips, tips_with_realm)

  cat("\nProcessing complete:\n")
  cat("  Tips with realm data:", length(tips_with_realm), "\n")
  cat("  Tips without realm data:", length(tips_without_realm), "(will be removed from tree)\n")

  if (length(tips_with_realm) == 0) {
    stop("ERROR: No tips have valid realm data. Cannot proceed.")
  }

  # Prune tree
  if (length(tips_without_realm) > 0) {
    tree <- ape::drop.tip(tree, tips_without_realm)
    cat("  Pruned tree now has:", length(tree$tip.label), "tips\n")
  }

  # Filter CSV to only include tips in pruned tree
  geog_expanded <- geog_expanded[geog_expanded$ID %in% tree$tip.label, ]

  # Statistics
  n_unique_tips <- length(unique(geog_expanded$ID))
  n_multi_realm <- sum(table(geog_expanded$ID) > 1)
  cat("  Tips with multiple realms:", n_multi_realm, "\n")
  cat("  Total rows in output CSV:", nrow(geog_expanded), "\n")

  # Determine output paths
  if (is.null(output_csv)) {
    output_csv <- file.path(dirname(geography_csv),
                            paste0(tools::file_path_sans_ext(basename(geography_csv)), "_biogeobears.csv"))
  }
  if (is.null(output_tree)) {
    output_tree <- sub("\\.(nwk|tre|tree)$", "_biogeobears.\\1", tree_path)
    if (output_tree == tree_path) {
      output_tree <- paste0(tools::file_path_sans_ext(tree_path), "_biogeobears.nwk")
    }
  }

  # Write outputs
  utils::write.csv(geog_expanded, file = output_csv, row.names = FALSE, quote = FALSE)
  ape::write.tree(tree, file = output_tree)

  cat("\nOutput files:\n")
  cat("  CSV:", output_csv, "\n")
  cat("  Tree:", output_tree, "\n\n")

  return(invisible(list(csv = output_csv, tree = output_tree)))
}

# Path to the FASTA file containing all barcode sequences. This includes barcodes 
# for the target family, distant groups, and sister families. For taxa with full 
# mitogenomes, their corresponding barcode regions are also included in this alignment.
# allsequences_path <- "C:\\Users\\16575\\Documents\\Rdocuments\\aligned_fin_k13k14_barcod_renamed_20250815.fasta"

# Path to the FASTA file containing the complete mitogenome sequences (supermatrix). 
# This file includes all mitogenomes used in the study, encompassing the distant and sister groups.
# mitogenomes_path <- "C:\\Users\\16575\\Documents\\Rdocuments\\5_Supermatrix_13PCG_NT_cleaned_rename_20250822.fasta"

# Reference to the list of sister groups for each family (defined previously).
# sister_group_family  

# Reference to the list of distant groups. These will serve as outgroups for rooting, 
# while sister groups will be part of the ingroup for dating calibrations.
# distantly_group


construction_tree <- function(nwk_file=NULL, sister_group_family, distantly_group, allsequences_path, mitogenomes_path,threads=14,construction_model="GTR+G") {
  library(ape)
  library(stringr)
  library(seqinr)
  
  # nwk_file: Path to a tree from a preliminary search, containing all tips. 
  # This step can be skipped if starting from sequences only.
  # sister_group_family: List containing sister groups and their divergence times, used for dating.
  # distantly_group: List of distant taxa, used for rooting.
  # allsequences_path: Path to the FASTA file with barcode region sequences for all relevant taxa.
  # mitogenomes_path: Path to the FASTA file with complete mitogenome sequences.
  # raxmlng_executable_path_wsl: Absolute path to the raxml-ng executable in WSL is recommended.
  # It can be found using 'which raxml-ng' in WSL, but often 'raxml-ng' alone works.
  
  if (is.null(nwk_file) || nwk_file == "") {
    # If no initial tree is provided, use the sequence file to define the initial tip set.
    nwk_file <- allsequences_path
    original_tips <- names(read.fasta(allsequences_path))
    # --- Setup paths and names for the current iteration ---
    output_dir <- dirname(nwk_file)
    base_name <- sub("\\.fasta$", "", basename(nwk_file))
    
  } else {
    # Read the provided tree to get the initial set of tip labels.
    current_tree <- read.tree(nwk_file)
    original_tips <- current_tree$tip.label
    # --- Setup paths and names for the current iteration ---
    output_dir <- dirname(nwk_file)
    base_name <- sub("\\.nwk$", "", basename(nwk_file))
    
  }
  
  cat("\nProcessing file:", basename(nwk_file), "\n")
  
  # Dynamically find the path to the raxml-ng executable within WSL.
  raxml_ng_path_raw <- tryCatch({
    system("wsl bash -ic 'which raxml-ng'", intern = TRUE)
  }, error = function(e) {
    stop("Could not find raxml-ng. Please ensure it is installed in WSL.")
  })
  if(length(raxml_ng_path_raw) == 0) {
    stop("raxml-ng not found in WSL path. Please check the installation.")
  }
  raxmlng_executable_path_wsl <- tail(raxml_ng_path_raw, 1)
  # Print the located path for verification.
  print(raxmlng_executable_path_wsl)
  
  # Define output file names.
  backbone_tree_path <- file.path(output_dir, paste0(base_name, "_backbone.nwk"))
  final_tree_path <- file.path(output_dir, paste0(base_name, "_reconstruction.nwk"))
  
  # Check if the final reconstruction tree already exists; if so, skip this process.
  if (file.exists(final_tree_path)) {
    cat("  --> Skipping:", final_tree_path, "(already exists)\n")
    return(invisible(NULL))
  }
  
  # Define paths for temporary FASTA files for this iteration.
  backbone_fasta_path <- file.path(output_dir, paste0(base_name, "_mmg_for_backbone.fasta"))
  final_fasta_path <- file.path(output_dir, paste0(base_name, "_all_for_final.fasta"))
  
  # Load all barcode-region sequences.
  all_seqs_fasta <- read.fasta(allsequences_path, seqtype = "DNA", as.string = TRUE)
  all_seqs_list <- setNames(all_seqs_fasta, getName(all_seqs_fasta))
  
  # Load all mitogenome sequences.
  mito_seqs_fasta <- read.fasta(mitogenomes_path, seqtype = "DNA", as.string = TRUE)
  mito_seqs_list <- setNames(mito_seqs_fasta, getName(mito_seqs_fasta))
  
  # Get all tip labels that correspond to mitogenomes (MMG).
  all_mito_tips <- names(mito_seqs_list)
  # Intersect with barcode tips, as not all mitogenomes may have a corresponding barcode sequence.
  all_mito_tips <- intersect(all_mito_tips, names(all_seqs_list))
  
  # Create a lookup table mapping each MMG tip to its family.
  # This is done by parsing tip labels directly for robustness.
  all_mito_tips_families <- sapply(all_mito_tips, extract_family, USE.NAMES = FALSE)
  names(all_mito_tips_families) <- all_mito_tips
  
  # Determine the target family based on the most frequent family in the original tip set.
  families <- sapply(original_tips, extract_family)
  target_family <- names(which.max(table(families)))
  
  if (is.na(target_family)) {
    warning(paste("Could not determine the target family for", nwk_file, ". Skipping."))
    return(invisible(NULL))
  }
  
  # Find sister families for the target family.
  sister_group_family1 <- lapply(sister_group_family, function(sister_vector) {
    sub("_\\d+$", "", sister_vector)
  })
  sister_families <- sister_group_family1[[target_family]]
  distantly_families <- distantly_group[[target_family]] # These will be added to the tip set later.
  if (is.null(sister_families)) {
    warning(paste("No sister families found for", target_family, "in the sister_group_family list. Skipping."))
    return(invisible(NULL))
  }
  
  # Select one representative MMG tip for each sister family.
  sister_reps_tips <- c()
  for (sf in sister_families) {
    # Find all MMG tips belonging to this sister family.
    possible_tips <- names(all_mito_tips_families[all_mito_tips_families == sf & !is.na(all_mito_tips_families)])
    if (length(possible_tips) > 0) {
      # Randomly select one representative.
      selected_tip <- sample(possible_tips, 1)
      sister_reps_tips <- c(sister_reps_tips, selected_tip)
      cat("  - Selected", selected_tip, "as representative for sister family", sf, "\n")
    } else {
      cat("  - No MMG found for sister family:", sf, "\n")
    }
  }
  
  if (length(sister_reps_tips) == 0) {
    warning(paste("No representative mitogenomes found for any sister family of", target_family, ". Cannot build a robust backbone. Skipping."))
    return(invisible(NULL))
  }
  
  # Select only one sister group representative to avoid potential dating issues,
  # especially if sister groups are closely related and their sequences are intermingled with the target family.
  sister_reps_tips <- sister_reps_tips[1]
  
  # Select representatives from distant families to serve as outgroups.
  possible_tips <- names(all_mito_tips_families[all_mito_tips_families %in% distantly_families & !is.na(all_mito_tips_families)])
  if (length(possible_tips) > 5) {
    # Randomly select five representatives.
    distantly_tips <- sample(possible_tips, 5)
    cat("  - Selected", length(distantly_tips), "tips as representatives for distant families.\n")
  } else if (length(possible_tips) > 0) {
    distantly_tips <- possible_tips
    cat("  - Selected all available", length(distantly_tips), "tips as representatives for distant families.\n")
  } else {
    warning(paste("No representative mitogenomes found for any distant family of", target_family, ". Cannot build a robust backbone. Skipping."))
    return(invisible(NULL))
  }
  
  # Combine original tips with the selected sister and distant group representatives.
  all_relevant_tips <- unique(c(original_tips, sister_reps_tips, distantly_tips))
  
  # --- B. Build Backbone Tree with Mitogenomes ---
  
  # Get all mitogenome tips from our combined list for the backbone tree.
  backbone_tips <- intersect(all_relevant_tips, all_mito_tips)
  
  if (length(backbone_tips) < 4) { # RAxML requires at least 4 taxa.
    warning(paste("Fewer than 4 mitogenomes available for backbone construction for", target_family, ". Skipping."))
    return(invisible(NULL))
  }
  
  # Extract their full mitogenome sequences and write to a FASTA file.
  backbone_sequences <- mito_seqs_list[backbone_tips]
  write.fasta(sequences = backbone_sequences, names = names(backbone_sequences), file.out = backbone_fasta_path)
  cat("  - Created backbone alignment with", length(backbone_tips), "taxa.\n")
  
  # Construct and run the RAxML-NG command for the backbone tree.
  # A standard GTR+G model is used.
  raxml_backbone_cmd <- sprintf(
    "wsl %s --all --msa %s --model %s --prefix %s --seed 2 --threads %d --force",
    raxmlng_executable_path_wsl,
    convert_path_to_wsl(backbone_fasta_path),
    construction_model,
    convert_path_to_wsl(file.path(output_dir, paste0(base_name, "_backbone"))), # RAxML adds its own suffix.
    threads
  )
  # The final tree with support values will be used.
  
  # The --all command includes bootstrap analysis, which can be slow but provides a complete result.
  # RAxML will generate a ".raxml.support" file which contains branch support values.
  expected_backbone_output <- file.path(output_dir, paste0(base_name, "_backbone.raxml.support"))
  
  # If a previous output file exists, use it instead of re-running the analysis.
  if (file.exists(expected_backbone_output)) {
    cat("  - Found existing RAxML output for backbone. Skipping run.\n")
  } else {
    cat("  - Running RAxML for backbone tree...\n")
    system(raxml_backbone_cmd)
  }
  
  # Check for the output file and copy it to the desired final path.
  if (file.exists(expected_backbone_output)) {
    file.copy(expected_backbone_output, backbone_tree_path, overwrite = TRUE)
  } else {
    warning("RAxML backbone reconstruction failed for", basename(nwk_file))
    return(invisible(NULL))
  }
  
  # --- C. Build Final Tree with Barcodes under Constraint ---
  
  # Extract barcode sequences for all relevant tips.
  final_sequences <- all_seqs_list[all_relevant_tips]
  # Remove any tips that might not be present in the barcode file.
  final_sequences <- final_sequences[!sapply(final_sequences, is.null)]
  write.fasta(sequences = final_sequences, names = names(final_sequences), file.out = final_fasta_path)
  cat("  - Created final alignment with", length(final_sequences), "taxa for constrained search.\n")
  
  # Construct and run the RAxML-NG command with the backbone tree as a topological constraint.
  raxml_final_cmd <- sprintf(
    "wsl %s --all --msa %s --tree-constraint %s --model %s --prefix %s --seed 2 --threads %d --force --redo",
    raxmlng_executable_path_wsl,
    convert_path_to_wsl(final_fasta_path),
    convert_path_to_wsl(backbone_tree_path),
    construction_model,
    convert_path_to_wsl(file.path(output_dir, paste0(base_name, "_reconstruction"))), # Prefix for the final run.
    threads
  )
  
  expected_final_output <- file.path(output_dir, paste0(base_name, "_reconstruction.raxml.support"))
  
  # If a previous output file exists, use it. Otherwise, run the command.
  if (file.exists(expected_final_output)) {
    cat("  - Found existing RAxML output for final tree. Skipping run.\n")
  } else {
    cat("  - Running RAxML for final constrained tree...\n")
    system(raxml_final_cmd)
  }
  
  # Check for the output file and copy it to the final path.
  if (file.exists(expected_final_output)) {
    file.copy(expected_final_output, final_tree_path, overwrite = TRUE)
    cat("  - Successfully created final reconstructed tree:", basename(final_tree_path), "\n")
  } else {
    warning("RAxML final reconstruction failed for", basename(nwk_file))
    return(invisible(NULL))
  }
  
  # Optional: Clean up temporary FASTA files.
  # file.remove(backbone_fasta_path)
  # file.remove(final_fasta_path)
}

# nwk_file<-"C:\\Users\\16575\\Documents\\Rdocuments\\test_20250703\\Cleridae_test_20250725/Cleridae_1.nwk"
# Execute the tree construction function.
# construction_tree(nwk_file = nwk_file,
#                   sister_group_family = sister_group_family,
#                   distantly_group = distantly_group,
#                   allsequences_path ="C:\\Users\\16575\\Documents\\Rdocuments\\aligned_fin_k13k14_barcod_renamed_20250815.fasta",
#                   mitogenomes_path = "C:\\Users\\16575\\Documents\\Rdocuments\\5_Supermatrix_13PCG_NT_cleaned_rename_20250822.fasta"
# )

# --- Step 2: Root the backbone and reconstruction trees ---
# The trees generated in the previous step (backbone and reconstruction) are currently unrooted.

# This function initiates the rooting process, which requires the 'distantly_group' list.
diagnose_and_prune_outgroup <- function(phy, outgroup_tips, action = "prune") {
  
  # Validate inputs.
  if (!inherits(phy, "phylo")) {
    stop("Error: 'phy' must be a 'phylo' object.")
  }
  if (!all(outgroup_tips %in% phy$tip.label)) {
    stop("Error: Some or all specified outgroup_tips are not present in the tree's tip labels.")
  }
  
  # Check if the outgroup is already monophyletic.
  if (is.monophyletic(phy, outgroup_tips)) {
    message("Congratulations! The specified outgroup is already monophyletic. No action needed.")
    return(phy)
  }
  
  message("Non-monophyletic outgroup detected. Searching for intruders...")
  
  # 1. Find the most recent common ancestor (MRCA) of all specified outgroup tips.
  mrca_node <- getMRCA(phy, outgroup_tips)
  
  # 2. Extract the subtree descending from this MRCA node to get all its tips.
  clade_subtree <- extract.clade(phy, mrca_node)
  clade_tips <- clade_subtree$tip.label
  
  # 3. Identify intruders by finding the set difference between all tips in the clade and the defined outgroup tips.
  intruders <- setdiff(clade_tips, outgroup_tips)
  
  # Print diagnostic information.
  if (length(intruders) > 0) {
    message(paste("Found", length(intruders), "intruder(s):"))
    print(intruders)
  } else {
    # This case is theoretically unlikely unless there's a subtle difference between is.monophyletic and getMRCA logic.
    message("No intruders directly found, but the group is non-monophyletic, suggesting a more complex topology.")
    return(phy)
  }
  
  # Execute the specified action based on the 'action' parameter.
  if (tolower(action) == "diagnose") {
    message("\nAction mode is 'diagnose'. Reporting intruders without modifying the tree.")
    return(intruders)
  } else if (tolower(action) == "prune") {
    message("\nAction mode is 'prune'. Removing intruders from the tree...")
    pruned_tree <- drop.tip(phy, tip = intruders)
    message("Pruning complete.")
    return(pruned_tree)
  } else {
    stop("Error: The 'action' parameter must be either 'diagnose' or 'prune'.")
  }
}

rooting_tree <- function(undated_tree, distantly_group,sister_group_family) {
  library(ape)
  library(stringr)
  
  # undated_tree: The path to an unrooted tree file.
  # distantly_group: A list defining the distant taxa to be used for rooting.
  
  use_nwk_tree <- undated_tree

  print(use_nwk_tree)
  use_tree <- read.tree(use_nwk_tree)
  original_tips <- use_tree$tip.label
  
  # Determine the target family from the tree's tip labels.
  families <- sapply(original_tips, extract_family, USE.NAMES = FALSE)
  target_family_name <- names(which.max(table(families)))
  
  distantly_family <- distantly_group[[target_family_name]]
  sister_group_family1 <- lapply(sister_group_family, function(sister_vector) {
    sub("_\\d+$", "", sister_vector)
  })
  sister_family<-sister_group_family1[[target_family_name]]
  
  # Select tips belonging to the designated distant (outgroup) families.
  selected_tips <- use_tree$tip.label[families %in% distantly_family]
  # Enforce that only mitogenome-derived taxa (MMG) are used for rooting.
  selected_tips <- grep("MMG_", selected_tips, value = TRUE)
  
  sister_tips<-use_tree$tip.label[families %in% sister_family]
  
  my_outgroups <- selected_tips
  print(my_outgroups)
  
  if (length(my_outgroups) == 0) {
    message(paste0(target_family_name, " does not have a suitable outgroup."))
    return(invisible(NULL))
  }
  
  # use_tree<-root(use_tree, outgroup = my_outgroups[1], resolve.root = TRUE)
  
  pure_tips<-my_outgroups
  # Diagnose and prune the tree to ensure the outgroup is monophyletic.
  pruned_tree <- diagnose_and_prune_outgroup(use_tree, pure_tips, action = "prune")
  
  if (length(sister_tips)>1){
    pruned_tree <- diagnose_and_prune_outgroup(pruned_tree, sister_tips, action = "prune")
    
  }
  
  # Root the tree using the cleaned outgroup. resolve.root = TRUE ensures a bifurcating root, which is standard practice.
  rooted_tree <- root(pruned_tree, outgroup = my_outgroups, resolve.root = TRUE)

  sister_mrca_node <- NULL # Initialize the node variable.
  
  if (length(sister_tips) > 1) {
    # If the sister group has multiple members, use getMRCA.
    sister_mrca_node <- ape::getMRCA(rooted_tree, sister_tips)
    
  } else if (length(sister_tips) == 1) {
    # If the sister group is a single tip, find its corresponding node index directly.
    # The which() function returns the position in tree$tip.label, which is the tip's node number.
    sister_mrca_node <- which(rooted_tree$tip.label == sister_tips)
    
  } else {
    stop("Error: The 'sister_tips' vector is empty.")
  }
  
  if (is.null(sister_mrca_node) || length(sister_mrca_node) == 0) {
    stop("Could not find the sister group node in the tree.")
  }
  
  # Find the parent node of the sister clade's MRCA.
  #    In tree$edge, column 2 lists all daughter nodes and column 1 lists their corresponding parent nodes.
  parent_of_sister_node <- rooted_tree$edge[rooted_tree$edge[, 2] == sister_mrca_node, 1]
  
  if (length(parent_of_sister_node) == 0) {
    # If the sister clade's MRCA is a direct child of the root, its parent is the root node itself.
    parent_of_sister_node <- length(rooted_tree$tip.label) + 1
  }
  
  # Define all tips that should be retained in the final tree.
  #    First, get all descendants from the parent node (this includes all correct target + sister tips).
  tips_from_correct_ingroup_clade <- phangorn::Descendants(rooted_tree, parent_of_sister_node, type="tips")
  tips_from_correct_ingroup_clade <- rooted_tree$tip.label[unlist(tips_from_correct_ingroup_clade)]
  
  #    Then, combine them with the clean, distant outgroup tips.
  tips_to_keep <- c(tips_from_correct_ingroup_clade, my_outgroups)
  
  # Identify "rogue" tips that are not in the 'tips_to_keep' list.
  rogue_tips_to_prune <- setdiff(rooted_tree$tip.label, tips_to_keep)
  
  # If any rogue tips are found, prune them from the tree.
  if (length(rogue_tips_to_prune) > 0) {
    print("Found and pruning the following misplaced 'rogue' target tips:")
    print(rogue_tips_to_prune)
    
    final_cleaned_tree <- ape::drop.tip(rooted_tree, rogue_tips_to_prune)
    
    print(paste("Successfully pruned", length(rogue_tips_to_prune), "rogue tips."))
    
  } else {
    print("No rogue tips found. The tree structure is valid.")
    final_cleaned_tree <- rooted_tree
  }
  
  # Write the final rooted tree to a new file.
  file_dir  <- dirname(use_nwk_tree)
  file_base <- basename(use_nwk_tree)
  file_base_clean <- sanitize_filename(file_base)
  use_nwk_tree <- file.path(file_dir, file_base_clean)
  rooted_file <-ifelse(grepl("\\.nwk$", use_nwk_tree),
                       sub("\\.nwk$", "_rooted.nwk", use_nwk_tree),
                       paste0(use_nwk_tree, "_rooted.nwk"))
  
  write.tree(final_cleaned_tree, file = rooted_file)
}


# --- Input Data for Rooting ---
# The unrooted tree specified below is an output from the 'construction_tree' function 
# (in this case, the backbone tree).
# undated_tree <- "C:/Users/16575/Documents/Rdocuments/test_20250703/Cleridae_test_20250725/Cleridae_1_atcg_GTR+G_backbone.nwk"
# distantly_group
# 
# # --- Execute the Rooting Function ---
# rooting_tree(undated_tree = undated_tree,
#              distantly_group = distantly_group,
#              sister_group_family = sister_group_family)


# --- Step 3: Perform Tree Dating ---
# The strategy is to first date the robust backbone tree using a sister group calibration.
# Then, the node ages from this dated backbone are used as secondary calibration points 
# to date the larger reconstruction tree, leveraging the high confidence of the backbone's topology and dating.
# rooted_tree_path <- "C:/Users/16575/Documents/Rdocuments/test_20250703/Cleridae_test_20250725/Cleridae_1_atcg_GTR+G_backbone_rooted.nwk"

dating_tree <- function(rooted_tree_path, sister_group_family, num_sites = 11171) {
  
  # treepl_executable_path_wsl: Defines the WSL path to the treePL executable. An absolute path, 
  # obtainable via 'which treePL' in WSL, is recommended for stability.
  # rooted_tree_path: The Windows path to the tree file to be dated (must be rooted).
  # sister_group_family: A list containing divergence time information for calibrations.
  # num_sites: The number of sites (base pairs) in the alignment used to build the tree.
  # e.g., ~11,000 for mitogenomes, ~658 for a standard barcode.
  
  # Dynamically find the path to the treePL executable within WSL.
  treepl_path_raw <- tryCatch({
    system("wsl bash -ic 'which treePL'", intern = TRUE)
  }, error = function(e) {
    stop("Could not find treePL. Please ensure it is installed in WSL.")
  })
  if(length(treepl_path_raw) == 0) {
    stop("treePL not found in WSL path. Please check the installation.")
  }
  treepl_executable_path_wsl <- tail(treepl_path_raw, 1)
  # Print the located path for verification.
  print(treepl_executable_path_wsl)
  
  sister_group_family1 <- lapply(sister_group_family, function(sister_vector) {
    sub("_\\d+$", "", sister_vector)
  })
  tree <- read.tree(rooted_tree_path)
  
  # --- Prepare information required by treePL ---
  
  # a. Identify ingroup and outgroup tips.
  all_tips <- tree$tip.label
  get_family_from_tip <- function(tip_label) {
    parts <- str_split(tip_label, "_")[[1]]
    if (length(parts) >= 3) return(parts[3])
    return(NA)
  }
  tip_families <- sapply(all_tips, get_family_from_tip, USE.NAMES = TRUE)
  
  original_tips <- tree$tip.label
  
  families <- sapply(original_tips, extract_family, USE.NAMES = F)
  target_family_name <- names(which.max(table(families)))
  
  NA_tips <- original_tips[families == ""]
  ingroup_tips <- original_tips[families == target_family_name]
  # Get sister family names (without divergence times) to identify the outgroup.
  sister_family_names <- sister_group_family1[[target_family_name]]
  outgroup_tips <- original_tips[families %in% sister_family_names]
  
  outgroup_tips_mmg <- outgroup_tips[grepl("^MMG_", outgroup_tips)]
  # Select only one MMG outgroup tip; other non-MMG or redundant MMG outgroup tips are pruned.
  outgroup_tips_mmg <- sample(outgroup_tips_mmg, 1)
  
  if (is.na(outgroup_tips_mmg) || length(outgroup_tips_mmg) == 0) {
    message(paste0(basename(rooted_tree_path), " doesn't have a suitable MMG outgroup for dating."))
    return(invisible(NULL))
  }
  
  outgroup_tips_mmg_family <- as.character(extract_family(outgroup_tips_mmg))
  # Identify tips to be removed: non-MMG outgroup tips and tips without family info.
  tips_to_delete <- setdiff(outgroup_tips, outgroup_tips_mmg)
  tips_to_delete <- unique(c(tips_to_delete, NA_tips))
  
  if (length(tips_to_delete) > 0) {
    tree <- drop.tip(tree, tips_to_delete)
  }
  
  # # b. Critical validation: Confirm that the ingroup is monophyletic (ensures correct rooting).
  # library(phytools)
  # if (!is.monophyletic(tree, ingroup_tips)) {
  #   warning(paste("  - Ingroup", target_family_name, "is NOT monophyletic. Rooting is incorrect for treePL. Skipping."))
  #   next
  # }
  # cat("  - Ingroup confirmed monophyletic. Rooting is valid for treePL.\n")
  
  # Extract the divergence time (root age) for the calibration.
  # Find the matching sister group entry.
  matched_item <- sister_group_family[[target_family_name]]
  divergence_time <- as.numeric(sub(".*_(\\d+)$", "\\1", matched_item[grep(outgroup_tips_mmg_family, matched_item)]))
  
  if (is.na(divergence_time) || length(divergence_time) == 0) {
    warning("  - Could not extract divergence time. Skipping.")
    return(invisible(NULL))
  }
  cat("  - Using root age calibration:", divergence_time, "Mya\n")
  cat("  - Number of sites for treePL:", num_sites, "\n")
  
  # Set a safe, non-zero minimum branch length acceptable for treePL.
  # Values like 1e-5 or 1e-6 are common choices.
  SAFE_MINIMUM_LENGTH <- 1e-5
  # Find all branches with lengths below this safe threshold.
  indices_to_fix <- which(tree$edge.length < SAFE_MINIMUM_LENGTH)
  if (length(indices_to_fix) > 0) {
    # Set all such branches to the safe minimum value.
    tree$edge.length[indices_to_fix] <- SAFE_MINIMUM_LENGTH
    cat("  - Found and elevated", length(indices_to_fix), "tiny/zero-length branches to a safe minimum value of", SAFE_MINIMUM_LENGTH, ".\n")
  } else {
    cat("  - No tiny/zero-length branches found below the threshold. The tree is clean.\n")
  }
  
  # Write the modified tree to a temporary file for treePL to use.
  fixed_tree_path_win <- tempfile(pattern = "dating", tmpdir = dirname(rooted_tree_path), fileext = ".nwk")
  write.tree(tree, fixed_tree_path_win)
  
  # The treePL config file will be created temporarily.
  config_file_path_win <- tempfile(pattern = "treepl_config_", tmpdir = dirname(rooted_tree_path), fileext = ".txt")
  on.exit(file.remove(config_file_path_win), add = TRUE) # Ensure cleanup on function exit.
  
  # Define the output path for the dated tree from treePL.
  output_dated_filename <- paste0(basename(sub("\\.nwk$", "", rooted_tree_path)), "_dated_treePL.nwk")
  output_dated_path_win <- file.path(dirname(rooted_tree_path), output_dated_filename)
  
  # Convert Windows paths to WSL paths for the config file.
  tree_file_path_wsl <- convert_path_to_wsl(fixed_tree_path_win)
  output_dated_path_wsl <- convert_path_to_wsl(output_dated_path_win)
  
  # Define the MRCA for the root calibration using ingroup and outgroup tips.
  # Use only MMG-derived tips for defining the root calibration.
  mrca_ingroup_tip <- ingroup_tips[grep("^MMG_", ingroup_tips)]
  mrca_outgroup_tip <- outgroup_tips_mmg
  
  # Assemble the content for the treePL configuration file.
  config_content <- c(
    paste("treefile =", tree_file_path_wsl),
    paste("outfile =", output_dated_path_wsl),
    paste("numsites =", num_sites),
    "",
    "# --- Calibration ---",
    "# Define the root node as the MRCA of the ingroup and outgroup",
    paste0("mrca = root_calibration ", paste(mrca_ingroup_tip, collapse=" "), " ", mrca_outgroup_tip),
    "# Constrain the root age (min and max are the same for a fixed age)",
    paste0("min = root_calibration ", divergence_time),
    paste0("max = root_calibration ", divergence_time),
    paste0("nthreads = 10"),
    "",
    "# --- Optimization Options ---",
    "prime          # Perform a primary analysis to find the best smoothing parameter",
    "#cv            # Uncomment to run cross-validation",
    "thorough       # Use a more thorough optimization routine",
    ""
  )
  
  # Write the configuration file.
  writeLines(config_content, config_file_path_win)
  cat("  - treePL config file generated at:", config_file_path_win, "\n")
  
  # Define the WSL path for the configuration file.
  config_file_path_wsl <- convert_path_to_wsl(config_file_path_win)
  
  # Prepare arguments for the system2 call. The command is 'wsl', 
  # and its arguments are the path to the treePL executable and the path to the config file.
  treepl_args <- c(
    treepl_executable_path_wsl,  # First argument passed to wsl
    config_file_path_wsl         # Second argument passed to wsl
  )
  
  cat("  - Running command: wsl", paste(treepl_args, collapse = " "), "\n")
  
  # Execute the command via system2. R handles quoting and spaces correctly.
  exit_code <- system2("wsl", args = treepl_args, stdout = TRUE, stderr = TRUE)
  
  # Check the exit status of the command.
  status_attr <- attr(exit_code, "status")
  if (!is.null(status_attr) && status_attr != 0) {
    warning("  - treePL process failed with exit status:", status_attr)
    print(exit_code) # Print error messages from treePL.
    return(invisible(NULL))
  }
  cat("  - treePL process completed.\n")
  
  # --- Post-processing and Saving ---
  if (file.exists(output_dated_path_win)) {
    cat("  - Dated tree found at:", output_dated_path_win, "\n")
    
    # Read the dated tree generated by treePL.
    dated_tree <- read.tree(output_dated_path_win)
    
    # Resolve polytomies and replace any remaining zero-length branches, 
    # which can represent rapid speciation events.
    if (!is.binary(dated_tree)) {
      dated_tree <- multi2di(dated_tree)
      cat("  - Tree converted to binary.\n")
    }
    zero_branches <- which(dated_tree$edge.length == 0)
    if (length(zero_branches) > 0) {
      dated_tree$edge.length[zero_branches] <- 1e-6 # Assign a very small branch length.
      cat("  - Replaced", length(zero_branches), "zero-length branches.\n")
    }
    
    # Overwrite the file with the finalized, cleaned tree.
    write.tree(dated_tree, file = output_dated_path_win)
    cat("  - Successfully finalized and saved dated tree to:", basename(output_dated_path_win), "\n")
  } else {
    warning("  - treePL execution seemed to succeed, but the output file was not found:", output_dated_path_win)
  }
}

# # Execute the dating function for the backbone tree.
# dating_tree(rooted_tree_path = rooted_tree_path,
#             sister_group_family = sister_group_family,
#             num_sites = 11171
# )

# backbone_dated_path <- "C:/Users/16575/Documents/Rdocuments/test_20250703/Cleridae_test_20250718/Cleridae_1_backbone_rooted_dated_treePL.nwk"
# reconstruction_path <- "C:/Users/16575/Documents/Rdocuments/test_20250703/Cleridae_test_20250718/Cleridae_1_reconstruction_rooted.nwk"

fix_backbone_dating_to_reconstruction <- function(backbone_dated_path, reconstruction_path, plot_result = TRUE, support_threshold = 50) {
  # backbone_dated_path: Path to the dated backbone tree.
  # reconstruction_path: Path to the un-dated reconstruction tree.
  # This function uses ape::chronos() to calibrate the reconstruction tree with node ages from the backbone tree.
  # It employs a Penalized Likelihood method.
  
  library(ape)
  
  # 2. Set file paths.
  backbone_dated_tree_path <- backbone_dated_path
  reconstruction_undating_tree_path <- reconstruction_path
  output_dir <- dirname(reconstruction_undating_tree_path)
  
  # 3. Read the tree files.
  print("Step 1: Loading tree files...")
  backbone_tree <- read.tree(backbone_dated_tree_path)
  recon_tree <- read.tree(reconstruction_undating_tree_path)
  
  if (any(backbone_tree$edge.length <=1e-8)) {
    cat("  - Found and fixed tiny/zero-length branches in the backbone tree.\n")
    backbone_tree$edge.length[backbone_tree$edge.length <= 1e-8] <- 1e-8
  }
  
  if (any(recon_tree$edge.length <= 1e-8)) {
    cat("  - Found and fixed tiny/zero-length branches in the reconstruction tree.\n")
    recon_tree$edge.length[recon_tree$edge.length <= 1e-8] <- 1e-6
  }
  
  # 4. Prepare the calibration data frame for chronos().
  print("Step 2: Preparing calibration data for chronos()...")
  
  # a) Get node age information from the dated backbone tree.
  backbone_ages_map <- setNames(ape::branching.times(backbone_tree),
                                (length(backbone_tree$tip.label) + 1):(length(backbone_tree$tip.label) + backbone_tree$Nnode))
  if (any(backbone_ages_map < 0)) {
    num_neg_nodes <- sum(backbone_ages_map < 0)
    cat(paste("  - Found", num_neg_nodes, "node(s) with negative age in backbone_tree, likely due to numerical precision. Fixing to 1e-6.\n"))
    backbone_ages_map[backbone_ages_map < 0] <- 1e-6
  }
  # b) Create an empty data frame for calibrations.
  # The required format for chronos is: node, age.min, age.max, soft.bounds
  calib_df <- data.frame()
  
  # c) Iterate through nodes of the backbone tree, find their correspondents in the reconstruction tree, and record their ages.
  for (node_id_backbone in names(backbone_ages_map)) {
    
    tips_subset <- ape::extract.clade(backbone_tree, as.numeric(node_id_backbone))$tip.label
    
    if (length(tips_subset) > 1 && all(tips_subset %in% recon_tree$tip.label)) {
      
      # Find the MRCA of these tips in the RECONSTRUCTION tree.
      node_id_recon <- ape::getMRCA(recon_tree, tips_subset)
      
      # If a corresponding node is found...
      if (!is.null(node_id_recon)) {
        
        # Get the age of the node from the backbone tree.
        node_age <- backbone_ages_map[node_id_backbone]
        
        # Add this information to the calibration data frame.
        # We set age.min and age.max to the same value for a fixed calibration point.
        new_row <- data.frame(node = node_id_recon,
                              age.min = node_age,
                              age.max = node_age,
                              soft.bounds = FALSE) # soft.bounds=FALSE indicates a hard constraint.
        calib_df <- rbind(calib_df, new_row)
      }
    }
  }
  
  # Remove duplicate calibration points (can occur if multiple backbone nodes map to the same reconstruction node).
  calib_df <- calib_df[!duplicated(calib_df$node), ]
  
  print(paste("Successfully prepared", nrow(calib_df), "calibration points for chronos()."))

  print("Step: Diagnosing and filtering conflicting calibration points...")
  
  # Conflicts can only exist if there is more than one calibration point.
  if (nrow(calib_df) > 1) {
    # For efficient lookup, create a named vector of node ages.
    nodes_to_check <- calib_df$node
    ages <- calib_df$age.min
    names(ages) <- nodes_to_check
    
    # Initialize a vector to store all nodes involved in conflicts.
    conflicting_nodes <- c()
    
    # Iterate through each calibration point, treating it as a potential ancestor.
    for (i in 1:length(nodes_to_check)) {
      ancestor_node <- nodes_to_check[i]
      ancestor_age <- ages[as.character(ancestor_node)]
      
      # Find all descendant nodes for the current ancestor in the reconstruction tree.
      # Requires the 'phangorn' package. Ensure it's installed.
      descendant_nodes <- phangorn::Descendants(recon_tree, ancestor_node, type = "all")
      
      # We are only interested in descendants that are also in our calibration list.
      calibrated_descendants <- intersect(descendant_nodes, nodes_to_check)
      
      # If any such descendants are found...
      if (length(calibrated_descendants) > 0) {
        # Get the ages of these calibrated descendants.
        descendant_ages <- ages[as.character(calibrated_descendants)]
        
        # Check for the paradox: "descendant is older than the ancestor".
        if (any(descendant_ages > ancestor_age)) {
          # If a paradox is found, record both the ancestor and all offending descendant nodes.
          conflicting_descendants <- calibrated_descendants[descendant_ages > ancestor_age]
          conflicting_nodes <- c(conflicting_nodes, ancestor_node, conflicting_descendants)
        }
      }
    }
    
    # Get the unique list of all nodes that need to be removed.
    conflicting_nodes <- unique(conflicting_nodes)
    
    if (length(conflicting_nodes) > 0) {
      print(paste("  - Found", length(conflicting_nodes), "conflicting calibration points. They will be removed."))
      # Remove all conflicting nodes from the calibration dataframe to get a "clean" version.
      calib_df_cleaned <- calib_df[!calib_df$node %in% conflicting_nodes, ]
    } else {
      # If no conflicts are found, use the original calibration dataframe.
      calib_df_cleaned <- calib_df
    }
    
    print(paste("  - Retaining", nrow(calib_df_cleaned), "non-conflicting calibration points for dating."))
  } else {
    # If there is only one or zero calibration points, no conflicts are possible.
    calib_df_cleaned <- calib_df 
  }
  
  # 5. Run the chronos() function to perform dating.
  print("Step 3: Running chronos() dating analysis... This may take some time.")
  
  # First, prepare the "fallback" calibration dataframe with soft bounds.
  calib_df_soft <- calib_df_cleaned
  calib_df_soft$age.min <- calib_df_cleaned$age.min * 0.95
  calib_df_soft$age.max <- calib_df_cleaned$age.max * 1.05
  calib_df_soft$soft.bounds <- TRUE
  
  dated_tree_chronos <- tryCatch({
    
    # --- Attempt A: First, try dating with strict hard constraints ---
    cat("  - Attempt 1: Running chronos with hard constraints...\n")
    ape::chronos(recon_tree, lambda = 1, model = "correlated", calibration = calib_df_cleaned)
    
  }, error = function(e) {
    
    # --- Attempt B (Fallback): If Attempt A fails, execute this block ---
    cat("  - WARNING: Dating with hard constraints failed. The error was:\n")
    cat(paste("    ", e$message, "\n"))
    cat("  - Attempt 2 (Fallback): Now trying again with flexible soft constraints...\n")
    
    # Use the prepared dataframe with soft constraints.
    ape::chronos(recon_tree, lambda = 1, model = "correlated", calibration = calib_df_soft)
    
  }) # End of tryCatch block

  print("Dating analysis completed!")
  
  dated_recon_tree<-dated_tree_chronos
  
  # ---  Transfer Bootstrap Values from Backbone to Reconstruction Tree ---
  print("Step 4: Transferring bootstrap support values from backbone...")
  
  if (is.null(dated_recon_tree$node.label)) dated_recon_tree$node.label <- rep("", dated_recon_tree$Nnode)
  if (is.null(backbone_tree$node.label)) stop("Backbone tree is missing bootstrap values (node labels).")
  
  support_origin <- rep("barcode", dated_recon_tree$Nnode) 
  
  for (i in 1:backbone_tree$Nnode) {
    backbone_node_id <- length(backbone_tree$tip.label) + i
    tips_subset <- ape::extract.clade(backbone_tree, backbone_node_id)$tip.label
    
    if (length(tips_subset) > 1 && all(tips_subset %in% dated_recon_tree$tip.label)) {
      recon_node_id <- ape::getMRCA(dated_recon_tree, tips_subset)
      
      if (!is.null(recon_node_id)) {
        backbone_bs_value <- backbone_tree$node.label[i]
        
        recon_node_label_idx <- recon_node_id - length(dated_recon_tree$tip.label)
        dated_recon_tree$node.label[recon_node_label_idx] <- backbone_bs_value
        support_origin[recon_node_label_idx] <- "mitogenome"
      }
    }
  }
  print("Bootstrap transfer completed.")
  
  # --- NEW: Optional and Enhanced Plotting ---
  if (plot_result) {
    print("Step 5: Generating detailed, annotated plot for Supplementary Material...")
    
    # --- Part A: Differentiate branch types (backbone vs. barcode-only) ---
    
    # This assumes your mitogenome-derived tips have a unique identifier, e.g., "MMG_" prefix.
    # !!! IMPORTANT: You may need to adapt this line to match your tip naming convention.
    backbone_tip_names <- dated_recon_tree$tip.label[grepl("^MMG_", dated_recon_tree$tip.label)]
    
    # Create a vector to store the color for each edge (branch)
    edge_colors <- rep("red3", nrow(dated_recon_tree$edge)) # Default color for barcode-only branches
    
    # Iterate through each edge to determine its type
    for (i in 1:nrow(dated_recon_tree$edge)) {
      daughter_node <- dated_recon_tree$edge[i, 2]
      
      # Get all descendant tips for the current branch
      descendant_tips <- NULL
      if (daughter_node <= length(dated_recon_tree$tip.label)) {
        # The branch leads directly to a single tip
        descendant_tips <- dated_recon_tree$tip.label[daughter_node]
      } else {
        # The branch leads to an internal node
        descendant_tips <- ape::extract.clade(dated_recon_tree, daughter_node)$tip.label
      }
      
      # If ANY descendant tip is from the original mitogenome backbone, 
      # then this branch is considered a backbone branch.
      if (any(descendant_tips %in% backbone_tip_names)) {
        edge_colors[i] <- "black" # Color for backbone branches
      }
    }
    
    # --- Part B: Prepare node support labels and colors ---
    
    # Use different colors for support values based on their origin
    node_support_colors <- ifelse(support_origin == "mitogenome", "blue", "darkgreen")
    
    # Prepare bootstrap values for plotting (only those above the threshold)
    bs_values <- as.numeric(dated_recon_tree$node.label)
    bs_to_plot <- ifelse(bs_values >= support_threshold, bs_values, "")
    
    # --- Part C: Generate the Plot (Adaptive Sizing Version) ---
    print("Step 5: Generating final annotated plot with adaptive sizing...")
    
    # --- START OF ADAPTIVE SIZING LOGIC ---
    
    # 1. Calculate a dynamic height based on the number of tips.
    #    This allocates a small amount of vertical space (e.g., 0.1 inches) to each tip.
    #    A base height (e.g., 10 inches) is added to ensure smaller trees are also readable.
    plot_height <- 10 + (Ntip(dated_recon_tree) * 0.1)
    
    # 2. Calculate a dynamic width based on the tree's depth (total time).
    #    This allocates horizontal space for the time axis.
    tree_depth <- max(ape::branching.times(dated_recon_tree))
    plot_width <- 10 + (tree_depth / 10) # Allocate 1 inch for every 10 million years.
    
    # 3. Calculate a dynamic font size (cex) to prevent labels from being too large or small.
    #    This formula reduces the font size as the number of tips increases.
    font_cex <- max(0.2, 1 - (Ntip(dated_recon_tree) / 1000))
    
    output_pdf_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(reconstruction_path)), "_annotated_detailed.pdf"))

    # Use the dynamically calculated dimensions.
    # The 'while' loop is good practice to close any stray open devices.
    while (!is.null(dev.list())) {
      dev.off()
    }
    pdf(output_pdf_path, width = plot_width, height = plot_height)
    
    # Plot the tree using the dynamic font size
    plot(dated_recon_tree, 
         type = "phylogram",
         show.tip.label = TRUE,
         edge.color = edge_colors,
         edge.width = 0.7,
         cex = font_cex)  # Use the dynamic font size
    
    axisPhylo()
    
    title(main = "Time-Calibrated Phylogeny with Tiered Node Support and Branch Provenance", 
          sub = "Branches: Black=Backbone, Red=Barcode-only | Support (BS >= 50%): Blue=Mito, Green=Barcode",
          cex.main = 1.2, cex.sub = 1.0, line = 2)
    
    # Use a slightly smaller font for node labels
    nodelabels(text = bs_to_plot, 
               frame = "none", 
               cex = font_cex * 0.8, # Base node label size on the tip label size
               col = node_support_colors)
    
    dev.off()
    print(paste("Detailed adaptive-sized plot saved to:", output_pdf_path))
  }
  
  
  # Save the final dated tree to a file.
  output_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(reconstruction_path)), "_dated_chronos.nwk"))
  ape::write.tree(dated_recon_tree, file = output_path)
  
  print(paste("Final dated tree has been saved to:", output_path))
  
}

# Execute the function to date the reconstruction tree using the backbone's ages.
# fix_backbone_dating_to_reconstruction(backbone_dated_path = backbone_dated_path,
#                                       reconstruction_path = reconstruction_path)



# --- Step 4: Perform Biogeographical Analysis ---
# This step includes analyses with PastML and BioGeoBEARS.
# The required inputs are a dated tree and a table of tip biogeographic states.

# =============================================================================
# cal_transition_probability_tab (formerly cal_transition_rate_tab)
# =============================================================================
#
# IMPORTANT METHODOLOGICAL NOTE:
# ------------------------------
# This function calculates "Transition Probability Contributions" from PastML
# marginal probability output. This is NOT equivalent to the "Transition Rate"
# calculated from BioGeoBEARS BSM (Biogeographical Stochastic Mapping).
#
# Key differences:
#   - PastML: Originally designed for pathogen (bacteria/virus) phylogeography
#             with short time scales and rapid generational turnover.
#   - PastML does NOT perform stochastic mapping or discrete event counting.
#   - This function approximates transition tendencies using:
#       Contribution = P(parent in state A) × P(child in state B)
#   - This is a PROBABILITY-BASED APPROXIMATION, not an event count.
#
# For rigorous Flux/Rate analysis (especially for macroevolutionary studies),
# use BioGeoBEARS BSM results instead.
#
# Output interpretation:
#   - "Total Contribution": Sum of probability products across all branches
#   - "Normalized Contribution": Total divided by tree branch length
#     (This is for normalization purposes only, NOT a true per-lineage rate)
#
# =============================================================================

cal_transition_probability_tab <- function(named_tree_path, probabilities_filepath) {
  library(pheatmap)

  cat("\n")
  cat("###########################################################################\n")
  cat("#                                                                         #\n")
  cat("#  PastML Transition Probability Contribution Analysis                    #\n")
  cat("#                                                                         #\n")
  cat("###########################################################################\n")
  cat("\n")
  cat("NOTE: This analysis uses marginal probabilities from PastML to estimate\n")
  cat("      transition tendencies. This is a PROBABILITY-BASED APPROXIMATION,\n")
  cat("      NOT equivalent to event-based rates from BioGeoBEARS BSM.\n")
  cat("      For rigorous Flux/Rate analysis, use BioGeoBEARS BSM results.\n")
  cat("\n")

  # --- 1. Load Data ---
  # Read the time-calibrated phylogenetic tree.
  phy_tree <- ape::read.tree(named_tree_path)

  # Read the node state probabilities output by PastML.
  # Assumes the first column of 'node_probabilities.tab' is the node/tip label,
  # and subsequent columns are probabilities for each geographic area.
  node_probs_df <- read.delim(probabilities_filepath, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

  # Define geographic areas (must match the order in the node_probs_df columns).
  realms <- colnames(node_probs_df)[2:ncol(node_probs_df)]

  # --- 2. Calculate Total Probability Contribution for Each Transition Type ---
  # Initialize a matrix to store the total probability contributions for all transitions.
  # Rows represent the 'from' state, and columns represent the 'to' state.
  num_realms <- length(realms)
  total_transition_contributions <- matrix(0,
                                           nrow = num_realms,
                                           ncol = num_realms,
                                           dimnames = list(from = realms, to = realms))

  # Iterate over every branch (edge) in the tree.
  # tree$edge is a two-column matrix: column 1 is the parent node, column 2 is the child node.
  for (i in 1:nrow(phy_tree$edge)) {
    parent_node_idx <- phy_tree$edge[i, 1]
    daughter_node_idx <- phy_tree$edge[i, 2]

    # --- Logic to retrieve state probabilities ---
    # This helper function maps an 'ape' node index to a node label
    # and retrieves the corresponding probability vector from the PastML output table.
    get_probs_by_ape_idx <- function(node_idx, tree, probs_df, realm_cols, node_id_col = "node") {
      node_label <- ""
      if (node_idx <= Ntip(tree)) {
        # It's a tip node.
        node_label <- tree$tip.label[node_idx]
      } else {
        # It's an internal node.
        # Ape node labels for nodes Ntip+1 to Ntip+Nnode are stored at indices 1 to Nnode.
        if (!is.null(tree$node.label)) {
          node_label <- tree$node.label[node_idx - Ntip(tree)]
        } else { # Fallback if internal nodes have no labels (use numeric index).
          node_label <- as.character(node_idx)
        }
      }
      probs_vec <- as.numeric(probs_df[probs_df[[node_id_col]] == node_label, realm_cols])
      if (length(probs_vec) == 0) {
        warning(paste("Probability data not found for node", node_label, "(index:", node_idx, "). Returning zero probability."))
        return(rep(0, length(realm_cols)))
      }
      return(probs_vec)
    }

    P_parent <- get_probs_by_ape_idx(parent_node_idx, phy_tree, node_probs_df, realms)
    P_daughter <- get_probs_by_ape_idx(daughter_node_idx, phy_tree, node_probs_df, realms)

    # If probability data for either node is missing, skip this branch.
    if (all(P_parent == 0) || all(P_daughter == 0)) {
      next
    }

    # Calculate and accumulate the contribution to each transition type.
    for (from_idx in 1:num_realms) {
      for (to_idx in 1:num_realms) {
        if (from_idx != to_idx) { # Consider only actual transitions (state changes).
          # Core logic: The contribution of this branch to a given transition (from -> to)
          # is the probability of the parent being in the 'from' state multiplied by
          # the probability of the daughter being in the 'to' state.
          # NOTE: This is a probability approximation, NOT an event count.
          contribution <- P_parent[from_idx] * P_daughter[to_idx]
          total_transition_contributions[from_idx, to_idx] <- total_transition_contributions[from_idx, to_idx] + contribution
        }
      }
    }
  }

  # --- 3. Calculate Normalized Contribution (for comparison purposes only) ---
  # Get the total branch length of the tree (unit: million years).
  total_branch_length_myrs <- sum(phy_tree$edge.length)

  if (total_branch_length_myrs == 0) {
    stop("Error: Total tree branch length is zero. Cannot normalize. Ensure the tree is a time-calibrated phylogram.")
  }

  # Normalized contribution = total probability contribution / total branch length.
  # NOTE: This is NOT a true "rate" - it's just normalized for comparison across trees of different sizes.
  normalized_contribution_matrix <- total_transition_contributions / total_branch_length_myrs

  # --- 4. Output Results ---
  cat("\n--- Results ---\n\n")

  cat("Total Transition Probability Contribution Matrix:\n")
  cat("(Sum of P_parent × P_daughter across all branches)\n\n")
  print(round(total_transition_contributions, 4))

  cat(paste("\nTotal phylogenetic branch length (Myr):", round(total_branch_length_myrs, 2), "\n"))

  cat("\nNormalized Transition Probability Contribution Matrix:\n")
  cat("(Total contribution / branch length - for cross-tree comparison only)\n")
  cat("WARNING: This is NOT equivalent to BioGeoBEARS transition rate!\n\n")
  print(round(normalized_contribution_matrix, 6))

  # --- 5. Save Results to CSV ---
  csv_total_save_path <- sub("\\_probabilities.tab$", "_transition_probability_total.csv", probabilities_filepath)
  csv_normalized_save_path <- sub("\\_probabilities.tab$", "_transition_probability_normalized.csv", probabilities_filepath)

  write.csv(as.data.frame(total_transition_contributions), file = csv_total_save_path, quote = FALSE)
  write.csv(as.data.frame(normalized_contribution_matrix), file = csv_normalized_save_path, quote = FALSE)
  cat(paste("\nTotal contribution matrix saved to:", csv_total_save_path, "\n"))
  cat(paste("Normalized contribution matrix saved to:", csv_normalized_save_path, "\n"))

  # --- 6. Visualization ---
  save_path <- sub("\\_probabilities.tab$", "_transition_probability_heatmap.png", probabilities_filepath)

  # Sort matrix rows and columns alphabetically for consistent comparison with BioGeoBEARS
  sorted_names <- sort(rownames(normalized_contribution_matrix))
  normalized_contribution_matrix <- normalized_contribution_matrix[sorted_names, sorted_names]

  # Convert numbers to scientific notation for display on the heatmap.
  numbers_sci <- formatC(normalized_contribution_matrix, format = "e", digits = 2)

  # Create heatmap with clear labeling
  pheatmap_obj <- pheatmap::pheatmap(
    normalized_contribution_matrix,
    display_numbers = numbers_sci,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "PastML: Normalized Transition Probability Contribution\n(NOT event-based rate - for visualization only)",
    fontsize_main = 11,
    filename = save_path,
    width = 8,
    height = 7,
    silent = TRUE
  )
  cat(paste("Heatmap saved to:", save_path, "\n"))

  cat("\n")
  cat("###########################################################################\n")
  cat("  REMINDER: For rigorous transition rate analysis, use BioGeoBEARS BSM.\n")
  cat("  PastML probability contributions are approximations only.\n")
  cat("###########################################################################\n")
  cat("\n")
}

# Keep old function name as alias for backward compatibility
cal_transition_rate_tab <- cal_transition_probability_tab


pastml_process_tree <- function(tree_path,
                                geography_csv = NULL) {

  library(tidyverse)
  library(ggtree)
  library(ggplot2)
  library(stringr)

  # =============================================================================
  # Step D: PastML Analysis
  # =============================================================================
  #
  # Input:
  #   - tree_path: Path to a clean, ready-to-analyze tree
  #   - geography_csv: Path to PastML-compatible CSV (from prepare_geography_for_pastml)
  #                    Format: ID, realm (one realm per tip, empty allowed)
  #                    If NULL, will attempt to extract from tip labels (legacy mode)
  #
  # Note:
  #   The tree should already be pruned to target family if needed.
  #   Use prune_tree_to_target_family() before calling this function.
  #
  # =============================================================================

  cat("\n")
  cat("###########################################################################\n")
  cat("#  Step D: PastML Ancestral State Reconstruction                         #\n")
  cat("###########################################################################\n\n")

  # Validate inputs
  if (!file.exists(tree_path)) {
    stop("ERROR: Tree file not found: ", tree_path)
  }

  # Create a dedicated directory for the tree's analysis outputs.
  tree <- tree_path
  win_tree_dir <- paste0(tools::file_path_sans_ext(tree), "_pastml")
  dir.create(win_tree_dir, showWarnings = FALSE)

  # Get the path to PastML within WSL.
  # Use pipe() with bash -ic (interactive) to load .bashrc which contains PATH
  pastml_path <- ""
  tryCatch({
    con <- pipe('wsl bash -ic "which pastml"', "r")
    lines <- readLines(con, warn = FALSE)
    close(con)
    # Find line containing pastml path (starts with /)
    for (line in lines) {
      line <- trimws(line)
      if (grepl("^/.*pastml$", line)) {
        pastml_path <- line
        break
      }
    }
  }, error = function(e) {
    pastml_path <<- ""
  })

  if (nchar(pastml_path) == 0) {
    stop("PastML not found in WSL. Please install with: pip install pastml")
  }

  cat("PastML path:", pastml_path, "\n")

  # Read tree
  use.tree <- read.tree(tree)
  use.tree.tip <- use.tree$tip.label
  tree_name <- tools::file_path_sans_ext(basename(tree))
  cat("Tree:", tree_path, "\n")
  cat("Number of tips:", length(use.tree.tip), "\n")

  # --- Load geography data ---
  if (!is.null(geography_csv) && file.exists(geography_csv)) {
    # Use provided CSV file
    cat("Geography CSV:", geography_csv, "\n")

    geog_raw <- utils::read.csv(geography_csv, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(geog_raw) <- tolower(trimws(colnames(geog_raw)))

    # Identify columns
    if ("id" %in% colnames(geog_raw)) {
      id_col <- "id"
    } else {
      id_col <- colnames(geog_raw)[1]
    }
    if ("realm" %in% colnames(geog_raw)) {
      realm_col <- "realm"
    } else if ("state" %in% colnames(geog_raw)) {
      realm_col <- "state"
    } else {
      realm_col <- colnames(geog_raw)[2]
    }

    # Create tips_realm vector (for PastML: take only FIRST realm if multiple)
    tips_realm <- setNames(rep("", length(use.tree.tip)), use.tree.tip)
    for (i in seq_len(nrow(geog_raw))) {
      tip_id <- geog_raw[[id_col]][i]
      if (tip_id %in% use.tree.tip) {
        realm_str <- geog_raw[[realm_col]][i]
        # Split by ; or , and take only the first realm (PastML requires single state)
        realms <- trimws(unlist(strsplit(as.character(realm_str), "[;,]")))
        realms <- realms[realms != "" & !is.na(realms) & tolower(realms) != "unk"]
        tips_realm[tip_id] <- if (length(realms) > 0) realms[1] else ""
      }
    }

  } else {
    # Legacy mode: extract from tip labels
    cat("WARNING: No geography CSV provided. Extracting from tip labels (legacy mode).\n")
    cat("         Consider using prepare_geography_for_pastml() for better control.\n\n")

    extract_realm_local <- function(tip) {
      parts <- str_split_fixed(tip, "_", 4)
      parts[,4]
    }
    tips_realm <- sapply(use.tree.tip, extract_realm_local)
  }

  # Map realm codes to single characters if needed
  if (!all(nchar(tips_realm[tips_realm != "" & !is.na(tips_realm)]) == 1)) {
    realm_map <- c(
      "NA" = "N", "NT" = "T", "PN" = "P", "PA" = "A",
      "SA" = "S", "SJ" = "J", "AT" = "F", "MA" = "M",
      "IM" = "I", "AA" = "U", "OC" = "O"
    )

    unique_realms <- unique(as.character(tips_realm))
    missing_realms <- setdiff(unique_realms, names(realm_map))
    missing_realms <- missing_realms[missing_realms != "" & !is.na(missing_realms)]
    extra_map <- setNames(letters[seq_along(missing_realms)], missing_realms)
    full_map <- c(realm_map, extra_map)

    tips_realm_mapped <- sapply(tips_realm, function(x) {
      if (x == "" || is.na(x)) "" else if (x %in% names(full_map)) full_map[[x]] else x
    })
    tips_realm <- tips_realm_mapped
  }

  # Save location file for PastML
  location_path <- paste0(win_tree_dir, "/matched_location.csv")
  tips_df <- data.frame(
    ID = names(tips_realm),
    realm = as.vector(tips_realm),
    stringsAsFactors = FALSE
  )
  write.csv(tips_df, file = location_path, row.names = FALSE)
  cat("Location file saved to:", location_path, "\n")

  # Statistics
  n_with_realm <- sum(tips_realm != "" & !is.na(tips_realm))
  n_without_realm <- sum(tips_realm == "" | is.na(tips_realm))
  cat("Tips with realm:", n_with_realm, "\n")
  cat("Tips without realm:", n_without_realm, "\n\n")

  safe_feature <- "realm"

  # Save tree for PastML (no pruning needed - tree is already clean)
  tree_save_path <- file.path(win_tree_dir, basename(tree))
  write.tree(use.tree, file = tree_save_path)

  # Convert paths to WSL format for command-line execution.
  wsl_location <- convert_path_to_wsl(location_path)
  wsl_tree <- convert_path_to_wsl(tree_save_path)
  wsl_work_dir <- convert_path_to_wsl(win_tree_dir)
  
  # PastML's 'ALL' method is a batch process but doesn't provide AICc values for model comparison.
  # Therefore, we define and run each combination separately.
  # For MPPA, the EFT and F81 models are equivalent, so we only need to run F81.
  all_combinations <- list(
    # MPPA + model
    list(method = "MPPA", model = "JC"),
    list(method = "MPPA", model = "F81"),
    
    # JOINT + model
    list(method = "JOINT", model = "JC"),
    list(method = "JOINT", model = "F81"),
    
    # Parsimony methods (no model)
    list(method = "DELTRAN", model = ""),
    list(method = "ACCTRAN", model = "")
  )
  
  for (combo in all_combinations) {
    method <- combo$method
    model <- combo$model
    suffix <- if (model != "") {
      paste0(method, "_", model)
    } else {
      method
    }
    
    win_html_out <- file.path(win_tree_dir, paste0(tree_name, "_pastml_", suffix, ".html"))
    win_csv_out <- file.path(win_tree_dir, paste0(tree_name, "_pastml_", suffix, ".csv"))
    wsl_html_out <- convert_path_to_wsl(win_html_out)
    wsl_csv_out <- convert_path_to_wsl(win_csv_out)
    
    # Construct the command string based on whether a model is specified.
    if (model != "") {
      cmd <- sprintf(
        'wsl %s --tree "%s" --data "%s" --columns "%s" --id_index 0 --prediction_method %s --data_sep "," --work_dir "%s" --html_mixed "%s" --out_data "%s" --model %s',
        pastml_path, wsl_tree, wsl_location, safe_feature, method, wsl_work_dir, wsl_html_out, wsl_csv_out, model
      )
    } else {
      cmd <- sprintf(
        'wsl %s --tree "%s" --data "%s" --columns "%s" --id_index 0 --prediction_method %s --data_sep "," --work_dir "%s" --html_mixed "%s" --out_data "%s"',
        pastml_path, wsl_tree, wsl_location, safe_feature, method, wsl_work_dir, wsl_html_out, wsl_csv_out
      )
    }
    cat("Running command:\n", cmd, "\n\n")
    
    # Execute the command.
    command_output <- system(cmd, intern = TRUE)
    
    if (method == "MPPA") {
      # Only MPPA provides marginal probabilities for each node, which are needed to calculate transition rates.
      # Rename the output probability and tree files for clarity and consistency.
      probabilities_filepath <- list.files(path = dirname(win_html_out), pattern = "\\.tab$", full.names = TRUE)
      probabilities_filepath <- probabilities_filepath[grep("marginal_probabilities.character", probabilities_filepath)]
      
      probabilities_new <- sub("\\.html$", "_probabilities.tab", win_html_out)
      file.rename(from = probabilities_filepath, to = probabilities_new)
      
      tree_path_tr <- list.files(path = dirname(win_html_out), pattern = "\\.nwk$", full.names = TRUE)
      named_tree_path <- tree_path_tr[grep("named\\.tree_", tree_path_tr)]
      named_tree_path_new <- sub("\\.html$", "_namedtree.nwk", win_html_out)
      file.rename(from = named_tree_path, to = named_tree_path_new)
      
      # Calculate and plot the transition rate matrix.
      cal_transition_rate_tab(named_tree_path = named_tree_path_new, probabilities_filepath = probabilities_new)
    }
    
    # Read the prediction data and reformat the node column.
    pred_data <- read.csv(win_csv_out) %>%
      separate(col = paste0("node.", safe_feature),
               into = c("node", safe_feature),
               sep = "\t") %>%
      mutate(node = as.character(node))
    
    # Read the tree file (with named internal nodes).
    tree_obj <- read.tree(named_tree_path_new)
    
    # Generate a base tree plot.
    p_base <- ggtree(tree_obj, layout = "circular") +
      theme_minimal()
    
    # Merge prediction data with the tree data.
    tree_data <- p_base$data %>%
      left_join(pred_data, by = c("label" = "node")) # Join using the 'label' column.
    
    # Create the final plot with tips and nodes colored by ancestral state.
    p <- ggtree(tree_obj, layout = "circular", aes(color = .data[[safe_feature]])) %<+%
      tree_data +
      geom_tippoint(aes(color = .data[[safe_feature]]), size = 2, alpha = 0.8) +
      theme(legend.position = "right") +
      labs(
        title = paste0("PastML Reconstruction (", suffix, ")")
      )
    
    tree_plot_save <- sub("\\.html$", "_tree_plot.png", win_html_out)
    
    ggsave(
      filename = tree_plot_save,
      plot = p,
      width = 16, height = 12, dpi = 300, units = "in"
    )
  }
  
  # Also run the 'ALL' method for a summary comparison report from PastML.
  wsl_ALL_work_dir <- file.path(wsl_work_dir, "ALL_methods")
  wsl_ALL_html_out <- file.path(wsl_ALL_work_dir, paste0("ALL_", tree_name, ".html"))
  wsl_ALL_csv_out <- file.path(wsl_ALL_work_dir, paste0("ALL_", tree_name, ".csv"))
  cmd_all <- sprintf(
    'wsl %s --tree "%s" --data "%s" --columns "%s" --id_index 0 --prediction_method ALL --data_sep "," --work_dir "%s" --html_mixed "%s" --out_data "%s" ',
    pastml_path,wsl_tree, wsl_location, safe_feature, wsl_ALL_work_dir, wsl_ALL_html_out, wsl_ALL_csv_out
  )

  system(cmd_all, intern = TRUE)
  
  print(paste0(tree_name, " finished"))
}

# Execute the PastML processing function.
# pastml_process_tree(dated_tree_path = dated_tree_path,
#                     discard_other_family = TRUE
# )



# --- BioGeoBEARS Analysis ---
# tree_filepath <- dated_tree_path
# # These file paths below are examples and will be overwritten by the function's arguments or auto-generated inputs.
# tip_states_filepath <- "C:\\Users\\16575\\Documents\\Rdocuments\\test_20250703\\matched_location_biogeobears.csv"
# timeperiods_filepath <- "C:\\Users\\16575\\Documents\\Rdocuments\\my_timeperiods_20mya.txt"
# dispersal_multipliers_filepath <- "C:\\Users\\16575\\Documents\\Rdocuments\\my_dispersal_multipliers_20mya.txt"

run_biogeobears_pipeline <- function(tree_filepath,
                                     geography_csv = NULL,  # Path to BioGeoBEARS-compatible CSV (from prepare_geography_for_biogeobears)
                                     output_tag = "biogeobears",
                                     num_cores = 12,
                                     models_to_run = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"),
                                     # The two files below are for time-stratified analysis.
                                     timeperiods_filepath = NULL,
                                     dispersal_multipliers_filepath = NULL,
                                     user_max_range_size = 3,
                                     min_brlen = 1e-6,
                                     plot_width_px = 6000,
                                     plot_height_px = 36000,
                                     plot_res_dpi = 400
) {
  # =============================================================================
  # Step E: BioGeoBEARS Analysis
  # =============================================================================
  #
  # Input:
  #   - tree_filepath: Path to a clean, ready-to-analyze tree
  #   - geography_csv: Path to BioGeoBEARS-compatible CSV (from prepare_geography_for_biogeobears)
  #                    Format: ID, realm (one realm per row, multiple rows per tip allowed)
  #                    If NULL, will attempt to extract from tip labels (legacy mode)
  #
  # Note:
  #   The tree should already be pruned appropriately.
  #   Use prune_tree_to_target_family() and prepare_geography_for_biogeobears() before calling.
  #
  # =============================================================================

  library(BioGeoBEARS)
  library(ape)
  library(dplyr)
  library(tidyr)
  library(tools) # For file path manipulation.
  library(phytools)
  library(stringr)
  library(snow)

  cat("\n")
  cat("###########################################################################\n")
  cat("#  Step E: BioGeoBEARS Biogeographical Analysis                          #\n")
  cat("###########################################################################\n\n")

  # --- 1. Validate inputs and set up paths ---
  # Path for caching temporary files.
  backup_dir <- file.path(dirname(tree_filepath), "backup_dir")
  dir.create(backup_dir, showWarnings = FALSE)

  # Path for final BioGeoBEARS results.
  save_biogeobears_path <- paste0(tools::file_path_sans_ext(tree_filepath), "_biogeobears_results")
  dir.create(save_biogeobears_path, showWarnings = FALSE)

  cat("--- 1. Validating inputs and setting up paths ---\n")
  if (!file.exists(tree_filepath)) stop("ERROR: Tree file not found: ", tree_filepath)
  cat("Tree:", tree_filepath, "\n")

  # Read tree
  use.tree <- read.tree(tree_filepath)
  use.tree.tip <- use.tree$tip.label
  cat("Number of tips:", length(use.tree.tip), "\n")

  # --- Load geography data ---
  tip_states_filepath <- paste0(save_biogeobears_path, "/matched_location.csv")

  if (!is.null(geography_csv) && file.exists(geography_csv)) {
    # Use provided CSV file
    cat("Geography CSV:", geography_csv, "\n")

    geog_raw <- utils::read.csv(geography_csv, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(geog_raw) <- tolower(trimws(colnames(geog_raw)))

    # Identify columns
    if ("id" %in% colnames(geog_raw)) {
      id_col <- "id"
    } else {
      id_col <- colnames(geog_raw)[1]
    }
    if ("realm" %in% colnames(geog_raw)) {
      realm_col <- "realm"
    } else if ("state" %in% colnames(geog_raw)) {
      realm_col <- "state"
    } else {
      realm_col <- colnames(geog_raw)[2]
    }

    # Filter to tips in tree and expand multi-realm entries (split "A;B" format)
    geog_filtered <- geog_raw[geog_raw[[id_col]] %in% use.tree.tip, ]

    # Expand multi-realm cells: each realm becomes a separate row
    tips_df <- data.frame(ID = character(), realm = character(), stringsAsFactors = FALSE)
    for (i in seq_len(nrow(geog_filtered))) {
      tip_id <- geog_filtered[[id_col]][i]
      realm_str <- geog_filtered[[realm_col]][i]
      # Split by ; or , to handle multi-realm format
      realms <- trimws(unlist(strsplit(as.character(realm_str), "[;,]")))
      realms <- realms[realms != "" & !is.na(realms) & tolower(realms) != "unk"]

      if (length(realms) > 0) {
        for (r in realms) {
          tips_df <- rbind(tips_df, data.frame(ID = tip_id, realm = r, stringsAsFactors = FALSE))
        }
      }
    }

    if (nrow(tips_df) == 0) {
      stop("ERROR: No valid realm data found after processing geography CSV.")
    }

    # Map realm codes to single characters if needed
    unique_realms <- unique(tips_df$realm)
    unique_realms <- unique_realms[unique_realms != "" & !is.na(unique_realms)]

    if (length(unique_realms) > 0 && !all(nchar(unique_realms) == 1)) {
      realm_map <- c(
        "NA" = "N", "NT" = "T", "PN" = "P", "PA" = "A",
        "SA" = "S", "SJ" = "J", "AT" = "F", "MA" = "M",
        "IM" = "I", "AA" = "U", "OC" = "O"
      )

      missing_realms <- setdiff(unique_realms, names(realm_map))
      missing_realms <- missing_realms[missing_realms != "" & !is.na(missing_realms)]
      extra_map <- setNames(toupper(letters[seq_along(missing_realms)]), missing_realms)
      full_map <- c(realm_map, extra_map)

      tips_df$realm <- sapply(tips_df$realm, function(x) {
        if (x == "" || is.na(x)) "" else if (x %in% names(full_map)) full_map[[x]] else x
      })
    }

    write.csv(tips_df, file = tip_states_filepath, row.names = FALSE)
    cat("Location file saved to:", tip_states_filepath, "\n")

    # Statistics
    n_unique_tips <- length(unique(tips_df$ID))
    n_multi_realm <- sum(table(tips_df$ID) > 1)
    cat("Tips with geography data:", n_unique_tips, "\n")
    cat("Tips with multiple realms:", n_multi_realm, "\n\n")

  } else {
    # Legacy mode: extract from tip labels
    cat("WARNING: No geography CSV provided. Extracting from tip labels (legacy mode).\n")
    cat("         Consider using prepare_geography_for_biogeobears() for better control.\n\n")

    extract_realm <- function(tip) {
      parts <- str_split_fixed(tip, "_", 4)
      parts[,4]
    }
    tips_realm <- sapply(use.tree.tip, extract_realm)

    # Map realm codes to single characters if needed
    if (!all(nchar(tips_realm[tips_realm != ""]) == 1)) {
      realm_map <- c(
        "NA" = "N", "NT" = "T", "PN" = "P", "PA" = "A",
        "SA" = "S", "SJ" = "J", "AT" = "F", "MA" = "M",
        "IM" = "I", "AA" = "U", "OC" = "O"
      )

      unique_realms <- unique(as.character(tips_realm))
      missing_realms <- setdiff(unique_realms, names(realm_map))
      missing_realms <- missing_realms[missing_realms != ""]
      extra_map <- setNames(letters[seq_along(missing_realms)], missing_realms)
      full_map <- c(realm_map, extra_map)

      tips_realm_mapped <- sapply(tips_realm, function(x) {
        if (x == "") "" else full_map[[x]]
      })
      tips_realm <- tips_realm_mapped
    }

    tips_df <- data.frame(
      ID = names(tips_realm),
      realm = as.vector(tips_realm),
      stringsAsFactors = FALSE
    )
    write.csv(tips_df, file = tip_states_filepath, row.names = FALSE)
  }

  if (!file.exists(tip_states_filepath)) stop("Tip states file not found: ", tip_states_filepath)

  output_dir <- save_biogeobears_path
  tree_basename <- tools::file_path_sans_ext(basename(tree_filepath))
  output_prefix <- file.path(output_dir,output_tag)
  
  print(paste("Input Tree:", tree_filepath))
  print(paste("Input States (auto-generated):", tip_states_filepath))
  print(paste("Output Directory:", output_dir))
  print(paste("Output File Prefix:", output_prefix))
  print("-------------------------------------------------")
  
  # --- 2. Reading input files ---
  print("--- 2. Reading input files ---")
  cs_tree <- ape::read.tree(tree_filepath)
  print(paste("Read tree with", length(cs_tree$tip.label), "tips and", cs_tree$Nnode, "internal nodes."))
  
  # Note: Tree pruning should be done beforehand using prune_tree_to_target_family()
  # This function now expects a pre-pruned tree and geography data.

  # Intelligently detect the separator (Tab or comma) in the states file.
  first_line_states <- readLines(tip_states_filepath, n = 1)
  sep_char <- ifelse(grepl("\t", first_line_states), "\t", ",")
  print(paste("Attempting to read states file with separator:", ifelse(sep_char == "\t", "TAB", "COMMA")))
  cs_states_raw <- tryCatch({
    utils::read.table(tip_states_filepath, header = TRUE, sep = sep_char, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
  }, error = function(e) {
    stop("Error reading tip states file '", tip_states_filepath, "': ", e$message)
  })
  colnames(cs_states_raw) <- c("node", "state")
  # Ensure the required columns ('node', 'state') are present.
  required_cols <- c("node", "state")
  if (!all(required_cols %in% colnames(cs_states_raw))) {
    stop("Tip states file must contain columns named 'node' and 'state'. Found: ", paste(colnames(cs_states_raw), collapse = ", "))
  }
  print("Read tip states file successfully.")
  print("-------------------------------------------------")
  
  # Here we assume the geography data is the master list; prune any tips from the tree that are not in this list.
  tip_labels_orig <- cs_tree$tip.label
  tips_to_delete <- setdiff(tip_labels_orig, cs_states_raw$node)
  if(length(tips_to_delete) > 0) {
    cs_tree <- drop.tip(cs_tree, tips_to_delete)
  }
  
  # --- 3. Processing geography data ---
  print("--- 3. Processing geography data ---")
  tip_labels_orig <- cs_tree$tip.label
  # Filter the state data to match the tips present in the (potentially pruned) tree.
  tip_states_data <- cs_states_raw %>%
    dplyr::select(all_of(required_cols)) %>% # Select only the required columns.
    dplyr::filter(node %in% tip_labels_orig)
  
  tip_states_data <- tip_states_data %>%
    filter(state != "")
  
  # Check for any tree tips that are missing from the state file (should not happen after the previous step, but is a good safeguard).
  missing_tips <- base::setdiff(tip_labels_orig, unique(tip_states_data$node))

  if (length(missing_tips) > 0) {
    # warning("The following tips in the tree have no state data: ", paste(missing_tips, collapse=", "))
    # stop("Cannot proceed with tips missing geographic data. Please provide complete data or prune the tree/data accordingly.")
    # Alternative: Pruning ( uncomment and adapt if needed)
    print(paste("Pruning", length(missing_tips), "tips with missing data from tree..."))
    cs_tree <- ape::drop.tip(cs_tree, tip = missing_tips)
    tip_labels_orig <- cs_tree$tip.label # Update tip labels
    tip_states_data <- tip_states_data %>% dplyr::filter(node %in% tip_labels_orig) # Ensure data matches pruned tree
  }
  
  # Get the final list of unique area names.
  unique_areas <- unique(as.character(tip_states_data$state))
  # Remove empty strings from unique areas if they exist
  unique_areas <- unique_areas[unique_areas != ""]
  num_areas <- length(unique_areas)
  if (num_areas == 0) stop("No valid geographic areas found in the state data for the tree tips.")
  print(paste("Identified", num_areas, "geographic areas:", paste(unique_areas, collapse = ", ")))
  
  # Convert the long-format data frame to the wide-format matrix required by BioGeoBEARS.
  geog_matrix_wide <- tip_states_data %>%
    dplyr::mutate(value = 1, state = as.character(state)) %>%
    dplyr::select(node, state, value) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = state, values_from = value, values_fill = 0, id_cols = node) %>%
    as.data.frame()
  
  # Set row names and ensure matrix rows are ordered to match the tree's tip labels.
  rownames(geog_matrix_wide) <- geog_matrix_wide$node
  if ("node" %in% colnames(geog_matrix_wide)) {
    geog_matrix_wide <- geog_matrix_wide[, -which(colnames(geog_matrix_wide) == "node"), drop = FALSE]
  }
  
  if (!all(tip_labels_orig %in% rownames(geog_matrix_wide))) {
    stop("Mismatch after pivoting: Some tree tips not found as rows in the geography matrix.")
  }
  
  geog_matrix_wide <- geog_matrix_wide[tip_labels_orig, , drop = FALSE]
  geog_matrix_wide <- na.omit(geog_matrix_wide)
  
  if (FALSE) {
    # This block is currently disabled. It was a previous method for handling tips with multiple distributions.
    # The current, preferred method is to list each state on a separate row in the input CSV file.
    # For a taxon with multiple geographic states, provide multiple rows in the input file.
    # Loop through sydysfinal_list and update the geography matrix (currently inactive).
    for (tip in names(sydysfinal_list)) {
      if (tip %in% rownames(geog_matrix_wide)) {
        print(tip)
        regions <- sydysfinal_list[[tip]]
        regions_in_matrix <- intersect(regions, colnames(geog_matrix_wide))
        geog_matrix_wide[tip, regions_in_matrix] <- 1
      }
    }
  }
  
  # Calculate the actual maximum number of areas any single tip occupies.
  actual_max_tip_range <- max(rowSums(geog_matrix_wide), na.rm = TRUE)
  print(paste("Actual maximum number of areas occupied by any single tip:", actual_max_tip_range))
  
  # Determine the max_range_size to be used for the analysis, ensuring it's valid.
  run_max_range_size <- NA
  if (!is.null(user_max_range_size)) {
    if (!is.numeric(user_max_range_size) || user_max_range_size < 1) {
      warning("Invalid user_max_range_size provided. Must be a positive integer.")
      user_max_range_size <- NULL # Reset to trigger default calculation.
    } else if (user_max_range_size < actual_max_tip_range) {
      warning(paste("User-provided max_range_size (", user_max_range_size, ") is less than the actual max tip range (", actual_max_tip_range, "). Adjusting to the minimum required.", sep = ""))
      run_max_range_size <- actual_max_tip_range
    } else {
      run_max_range_size <- user_max_range_size
    }
  }
  # If run_max_range_size is still not set, default to the minimum required size.
  if (is.na(run_max_range_size)) {
    run_max_range_size <- actual_max_tip_range
  }
  print(paste("Setting max_range_size for BioGeoBEARS run to:", run_max_range_size))
  
  # Write the processed geography data to a temporary PHYLIP-formatted file.
  geog_fn_temp <- tempfile(tmpdir = backup_dir, fileext = "_geog.txt")
  num_taxa <- nrow(geog_matrix_wide)
  area_names_actual <- colnames(geog_matrix_wide)
  header_line <- paste(num_taxa, num_areas, paste0("(", paste(area_names_actual, collapse = " "), ")"), sep = "\t")
  cat(header_line, "\n", file = geog_fn_temp)
  
  # This vector will track if any tip labels need to be modified (e.g., spaces replaced).
  tree_tip_labels_to_use <- cs_tree$tip.label
  any_label_changed <- FALSE
  
  print("Writing geography data to temporary file...")
  for (i in 1:num_taxa) {
    taxon_name <- rownames(geog_matrix_wide)[i]
    taxon_name_nospace <- gsub(" ", "_", taxon_name) # Replace spaces with underscores.
    if (taxon_name != taxon_name_nospace) {
      warning(paste("Tip label '", taxon_name, "' converted to '", taxon_name_nospace, "' for geography file.", sep = ""))
      # Update the vector of tip labels to be used in the tree object.
      match_indices <- which(tree_tip_labels_to_use == taxon_name)
      if (length(match_indices) > 0) {
        tree_tip_labels_to_use[match_indices] <- taxon_name_nospace
        any_label_changed <- TRUE
      }
    }
    area_string <- paste(as.numeric(geog_matrix_wide[i, ]), collapse = "")
    cat(paste(taxon_name_nospace, "\t", area_string), "\n", file = geog_fn_temp, append = TRUE)
  }
  print(paste("Temporary geography data written to:", geog_fn_temp))
  print("-------------------------------------------------")

  # --- 4. Processing the Phylogenetic Tree ---
  print("--- 4. Processing phylogenetic tree ---")
  tr_processed <- cs_tree # Start with the tree loaded in the previous step.
  
  # If tip labels were changed (e.g., spaces replaced), update the tree object.
  if (any_label_changed) {
    print("Updating tree tip labels to match modified names (spaces to underscores)...")
    tr_processed$tip.label <- tree_tip_labels_to_use
  }
  
  # Ensure the tree is fully bifurcating (binary).
  if (!ape::is.binary(tr_processed)) {
    print("Resolving polytomies using multi2di()...")
    tr_processed <- ape::multi2di(tr_processed)
  } else {
    print("Tree is binary.")
  }
  # Ensure the tree is rooted.
  if (!ape::is.rooted(tr_processed)) {
    tr_processed <- phangorn::midpoint(tr_processed) # Fallback to midpoint rooting if necessary.
  }
  # Fix branch lengths after resolving polytomies.
  print(paste("Fixing non-positive branch lengths to a minimum of", min_brlen))
  tr_processed$edge.length[tr_processed$edge.length <= 0] <- min_brlen
  
  # Final check to ensure no zero-length branches remain.
  if (any(tr_processed$edge.length <= 0)) stop("FATAL ERROR: Tree fixing failed! Non-positive branch lengths still exist.")

  # --- Pre-check: Validate tree is time-calibrated before time-stratified analysis ---
  tree_age <- max(ape::branching.times(tr_processed))

  if (!is.null(timeperiods_filepath) && !is.null(dispersal_multipliers_filepath)) {
    # Check 1: Tree age is suspiciously small (likely not dated)
    if (tree_age < 1) {
      cat("\n")
      cat("###########################################################################\n")
      cat("#                    WARNING: TREE MAY NOT BE DATED                      #\n")
      cat("###########################################################################\n")
      cat("\n")
      cat("  Tree root age:", round(tree_age, 4), "Ma\n")
      cat("\n")
      cat("  This is suspiciously small. The tree may not be time-calibrated.\n")
      cat("  Time-stratified BioGeoBEARS analysis REQUIRES a dated tree.\n")
      cat("\n")
      cat("  Possible causes:\n")
      cat("    - Tree branch lengths are in substitution units, not time units\n")
      cat("    - Molecular clock dating was not performed\n")
      cat("    - Dating failed silently\n")
      cat("\n")
      cat("  Solutions:\n")
      cat("    1. Run molecular clock dating (e.g., treePL, chronos, BEAST)\n")
      cat("    2. Or skip time-stratified analysis by setting:\n")
      cat("       timeperiods_filepath = NULL\n")
      cat("       dispersal_multipliers_filepath = NULL\n")
      cat("\n")
      cat("###########################################################################\n")
      cat("\n")
      stop("Tree root age (", round(tree_age, 4), " Ma) is too small. Please provide a time-calibrated tree for time-stratified analysis.")
    }

    # Check 2: Tree age vs time periods compatibility
    time_boundaries_check <- as.numeric(readLines(timeperiods_filepath))
    min_time_boundary <- min(time_boundaries_check[time_boundaries_check > 0])

    if (tree_age < min_time_boundary) {
      cat("\n")
      cat("###########################################################################\n")
      cat("#              ERROR: TREE AGE vs TIME PERIODS MISMATCH                  #\n")
      cat("###########################################################################\n")
      cat("\n")
      cat("  Tree root age:", round(tree_age, 2), "Ma\n")
      cat("  Smallest time boundary:", min_time_boundary, "Ma\n")
      cat("\n")
      cat("  The tree is younger than your time stratification boundaries.\n")
      cat("  BioGeoBEARS time-stratified analysis cannot proceed.\n")
      cat("\n")
      cat("  Solutions:\n")
      cat("    1. Use a properly dated tree with appropriate age\n")
      cat("    2. Adjust your time periods file to match tree age\n")
      cat("    3. Skip time-stratified analysis by setting:\n")
      cat("       timeperiods_filepath = NULL\n")
      cat("       dispersal_multipliers_filepath = NULL\n")
      cat("\n")
      cat("###########################################################################\n")
      cat("\n")
      stop("Tree root age (", round(tree_age, 2), " Ma) is younger than the smallest time boundary (", min_time_boundary, " Ma).")
    }

    # Check 3: Negative branching times (indicates serious problem)
    all_bt <- ape::branching.times(tr_processed)
    if (any(all_bt < 0)) {
      n_negative <- sum(all_bt < 0)
      cat("\n")
      cat("###########################################################################\n")
      cat("#              WARNING: NEGATIVE BRANCHING TIMES DETECTED                #\n")
      cat("###########################################################################\n")
      cat("\n")
      cat("  Found", n_negative, "nodes with negative branching times.\n")
      cat("  This indicates the tree has serious problems with branch lengths.\n")
      cat("  The tree may not be ultrametric or properly dated.\n")
      cat("\n")
      cat("###########################################################################\n")
      cat("\n")
      stop("Tree has ", n_negative, " nodes with negative branching times. Please check tree dating.")
    }

    cat("  [OK] Tree age check passed:", round(tree_age, 2), "Ma\n")
  }

  # Proactively fix a known issue in BioGeoBEARS for time-stratified analysis.
  if (TRUE) {
    # An error can occur if a node age in the tree is exactly equal to a time-slice boundary.
    # This block identifies such nodes and slightly adjusts their ages by a trivial amount (e.g., 0.01 Myr).
    
    # --- 4a. Set Parameters for Node Age Correction ---
    time_boundaries <- as.numeric(readLines(timeperiods_filepath))
    # Define a small tolerance value to identify nodes "too close" to a boundary.
    # The BioGeoBEARS default is 1e-5; we use a slightly larger value for safety.
    tolerance <- 1e-4
    # The small amount of time to add/subtract from branches to shift a node's age.
    amount_to_add <- 0.01 # A small value is recommended to avoid significantly altering dates.
    
    # --- 4b. Helper function to adjust node ages ---
    fix_node_age <- function(phy, node_to_fix, amount_to_add = 0.01) {
      
      # --- Safety checks ---
      if (node_to_fix <= length(phy$tip.label)) {
        stop("Error: 'node_to_fix' must be an internal node number.")
      }
      if (node_to_fix == (length(phy$tip.label) + 1)) {
        stop("Error: This function cannot be used on the root node. Use the separate root-handling logic.")
      }
      
      # 1. Find the indices of the parent and two daughter edges.
      parent_edge_idx <- which(phy$edge[, 2] == node_to_fix)
      if (length(parent_edge_idx) == 0) {
        stop(paste("Could not find the parent branch for node", node_to_fix))
      }
      child_edges_idx <- which(phy$edge[, 1] == node_to_fix)
      if (length(child_edges_idx) != 2) {
        stop(paste("Node", node_to_fix, "is not a simple bifurcation and cannot be processed."))
      }
      
      # --- Final safety check ---
      # Ensure daughter branches are long enough to be shortened.
      if (any(phy$edge.length[child_edges_idx] < amount_to_add)) {
        stop(paste("Error: A daughter branch of node", node_to_fix, "is too short to subtract", amount_to_add, ". Try a smaller 'amount_to_add' value."))
      }
      
      cat(sprintf("Correcting node #%d...\n", node_to_fix))
      
      # 2. Perform the "push-pull" operation on branch lengths.
      # Lengthen the parent branch (pushes node deeper in time).
      phy$edge.length[parent_edge_idx] <- phy$edge.length[parent_edge_idx] + amount_to_add
      cat(sprintf("  - Parent branch (edge %d) length increased by %f\n", parent_edge_idx, amount_to_add))
      
      # Shorten the daughter branches (pulls descendants forward in time).
      phy$edge.length[child_edges_idx] <- phy$edge.length[child_edges_idx] - amount_to_add
      cat(sprintf("  - Daughter branches (edges %d, %d) length decreased by %f\n", child_edges_idx[1], child_edges_idx[2], amount_to_add))
      
      return(phy)
    }
    
    # --- 4c. Iterate through all internal nodes to check and correct them ---
    num_tips <- length(tr_processed$tip.label)
    num_nodes <- tr_processed$Nnode
    internal_nodes <- (num_tips + 1):(num_tips + num_nodes)
    
    for (node_num in internal_nodes) {
      
      # Node heights must be recalculated in each iteration as the tree may have been modified.
      all_heights <- nodeHeights(tr_processed)
      max_height <- max(all_heights) # Current total tree height (age of the root).
      
      # Calculate the age of the current node (time before present).
      # Node height is calculated from the tips (time 0).
      # The root node does not appear in edge[,2], so it requires special handling.
      if (node_num == (num_tips + 1)) { # This is the root node.
        node_height <- max_height # The root's height is the max height.
      } else {
        node_height <- all_heights[which(tr_processed$edge[, 2] == node_num), 2]
      }
      node_age <- max_height - node_height
      
      # Check if the node's age is too close to any time boundary.
      distance_to_boundaries <- abs(node_age - time_boundaries)
      
      if (any(distance_to_boundaries < tolerance)) {
        
        # A problematic node was found.
        problem_boundary <- time_boundaries[which.min(distance_to_boundaries)]
        cat(sprintf("  - Found node #%s (age: %f Ma) is too close to boundary %s Ma.\n",
                    node_num, node_age, problem_boundary))
        
        # Correct the node's age.
        if (node_num == (num_tips + 1)) { # Root node handling
          cat("    -> Correcting root node...\n")
          root_edge_indices <- which(tr_processed$edge[, 1] == node_num)
          # Lengthen the two daughter branches stemming from the root.
          tr_processed$edge.length[root_edge_indices] <- tr_processed$edge.length[root_edge_indices] + amount_to_add
        } else { # Internal node handling
          cat("    -> Correcting internal node...\n")
          tr_processed <- fix_node_age(phy = tr_processed, node_to_fix = node_num, amount_to_add = amount_to_add)
        }
        
        # Verify the age after correction.
        new_max_height <- max(nodeHeights(tr_processed))
        cat(sprintf("    -> Tree has been modified. New total height is: %f Ma\n", new_max_height))
      }
    }
  } # End of proactive node age correction.
  
  print("Tree processing complete.")
  tree_age <- max(ape::branching.times(tr_processed))
  print(paste("Current tree age (root depth):", round(tree_age, 2), "Ma"))
  
  # Write the fully processed tree to a temporary file.
  tree_fn_temp <- tempfile(fileext = "_tree.nwk", tmpdir = backup_dir)
  ape::write.tree(tr_processed, file = tree_fn_temp)
  print(paste("Processed tree written to temporary file:", tree_fn_temp))
  print("-------------------------------------------------")
  
  # --- 5. BioGeoBEARS Analysis Setup ---
  print("--- 5. Setting up base BioGeoBEARS run object ---")
  BioGeoBEARS_run_object_base <- BioGeoBEARS::define_BioGeoBEARS_run()
  BioGeoBEARS_run_object_base$trfn <- tree_fn_temp
  BioGeoBEARS_run_object_base$geogfn <- geog_fn_temp
  BioGeoBEARS_run_object_base$max_range_size <- run_max_range_size
  BioGeoBEARS_run_object_base$min_branchlength <- min_brlen
  BioGeoBEARS_run_object_base$include_null_range <- TRUE
  BioGeoBEARS_run_object_base$num_cores_to_use <- num_cores
  BioGeoBEARS_run_object_base$force_sparse <- FALSE
  BioGeoBEARS_run_object_base$on_NaN_error <- -1e50
  BioGeoBEARS_run_object_base$speedup <- TRUE
  BioGeoBEARS_run_object_base$use_optimx <- TRUE
  BioGeoBEARS_run_object_base$min_dist_between_node_and_stratum_line <- 1e-9
  
  # These functions read the input files and validate the base run object.
  # Note: `readfiles` will be called again inside the loop for time-stratified setup.
  BioGeoBEARS_run_object_base <- BioGeoBEARS::readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object_base)
  BioGeoBEARS_run_object_base <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object_base)
  BioGeoBEARS::check_BioGeoBEARS_run(BioGeoBEARS_run_object_base)
  
  # --- 6. Model Comparison Loop ---
  print(paste("--- 6. Starting model comparison for:", paste(models_to_run, collapse = ", "), "---"))
  
  # Create a list to store the results of each model run.
  results_list <- list()
  
  for (model_name in models_to_run) {
    print(paste("--- Setting up and running model:", model_name, "---"))
    
    # Copy the base configuration for the current model.
    run_object <- BioGeoBEARS_run_object_base
    
    # Configure for time-stratified analysis if files are provided.
    if (!is.null(timeperiods_filepath) && !is.null(dispersal_multipliers_filepath)) {
      
      print("Time-stratified analysis is enabled.")
      
      # "Purify" the master dispersal matrix to include only the areas present in the current dataset.
      if (TRUE) {
        # --- Automated Matrix Purification Process ---
        # 1. Read and parse the master dispersal multipliers file, which may contain areas not in this analysis.
        print("Step 1: Reading and parsing the master dispersal multipliers file...")
        master_filepath <- dispersal_multipliers_filepath
        all_lines <- readLines(master_filepath)
        blank_lines_indices <- c(0, which(all_lines == ""), length(all_lines) + 1)
        master_mats_list <- list()
        for (i in 1:(length(blank_lines_indices) - 1)) {
          start_index <- blank_lines_indices[i] + 1
          end_index <- blank_lines_indices[i+1] - 1
          if (start_index > end_index) next
          matrix_block_text <- all_lines[start_index:end_index]
          header <- trimws(strsplit(matrix_block_text[1], "\t")[[1]])
          data_mat <- as.matrix(utils::read.table(text = matrix_block_text[-1], sep = "\t"))
          colnames(data_mat) <- header
          rownames(data_mat) <- header
          master_mats_list[[i]] <- data_mat
        }
        print(paste("Successfully parsed", length(master_mats_list), "master matrices from the file."))
        
        # 2. Subset the master matrices to keep only the relevant areas for this specific dataset.
        print("Step 2: Subsetting matrices to keep only the relevant areas...")
        areas_to_keep <- unique_areas
        subset_mats_list <- list()
        for (master_mat in master_mats_list) {
          if (!all(areas_to_keep %in% colnames(master_mat))) {
            stop("Error: Some areas in 'areas_to_keep' are not found in the master matrix header.")
          }
          subset_mat <- master_mat[areas_to_keep, areas_to_keep]
          subset_mats_list[[length(subset_mats_list) + 1]] <- subset_mat
        }
        
        # 3. Filter time periods and matrices to only those relevant to the tree's age.
        master_time_periods_list <- readLines(timeperiods_filepath)
        filtered_time_periods <- list()
        filtered_dispersal_mats <- list()
        for (i in 1:length(master_time_periods_list)) {
          period <- as.numeric(master_time_periods_list[[i]])
          if (tree_age > period) {
            filtered_time_periods[[length(filtered_time_periods) + 1]] <- period
            filtered_dispersal_mats[[length(filtered_dispersal_mats) + 1]] <- subset_mats_list[[i]]
          } else {
            # Add the first time period that is older than the tree and then stop.
            filtered_time_periods[[length(filtered_time_periods) + 1]] <- period
            filtered_dispersal_mats[[length(filtered_dispersal_mats) + 1]] <- subset_mats_list[[i]]
            break
          }
        }
        
        # 4. Write the purified/subsetted files to temporary locations.
        # Time periods file
        time_periods_filepath_temp <- tempfile(tmpdir = backup_dir, pattern = "timeperiods_", fileext = ".txt")
        time_boundaries <- sort(unique(unlist(filtered_time_periods)))
        time_boundaries_no_zero <- time_boundaries[time_boundaries > 0]
        writeLines(text = as.character(time_boundaries_no_zero), con = time_periods_filepath_temp)
        
        # Dispersal multipliers file
        subset_filepath <- tempfile(pattern = "dispersal_mats_", tmpdir = backup_dir, fileext = ".txt")
        print(paste("Step 3: Writing the subsetted matrices to a new temporary file:", subset_filepath))
        cat("", file = subset_filepath) # Clear file.
        for (current_mat in filtered_dispersal_mats) {
          cat(paste(colnames(current_mat), collapse = "\t"), "\n", file = subset_filepath, append = TRUE)
          write.table(current_mat, file = subset_filepath, append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
          cat("\n", file = subset_filepath, append = TRUE)
        }
        print("Purification of time-stratified inputs complete!")
      }
      
      # Set the paths in the run object to the temporary, purified files.
      run_object$timesfn <- time_periods_filepath_temp
      run_object$dispersal_multipliers_fn <- subset_filepath
      
      # Critical step for time-stratified analysis: section_the_tree() must be run AFTER readfiles_BioGeoBEARS_run().
      run_object <- BioGeoBEARS::readfiles_BioGeoBEARS_run(run_object)
      run_object <- section_the_tree(inputs = run_object, make_master_table = TRUE, plot_pieces = FALSE)
      
    } else {
      # For a non-time-stratified analysis, the files have already been read.
      # This else block is technically redundant given the base setup, but is good practice.
      run_object <- BioGeoBEARS::readfiles_BioGeoBEARS_run(run_object)
    }
    
    # Configure the parameters for the specific model being run.
    if (model_name == "DEC") {
      # DEC is the default model; no changes needed.
    } else if (model_name == "DEC+J") {
      run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
      run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- 0.01
    } else if (model_name == "DIVALIKE") {
      run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
      run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
      run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "free"
      run_object$BioGeoBEARS_model_object@params_table["v", "init"] <- 1.0
    } else if (model_name == "DIVALIKE+J") {
      run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
      run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
      run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "free"
      run_object$BioGeoBEARS_model_object@params_table["v", "init"] <- 1.0
      run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
      run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- 0.01
      run_object$BioGeoBEARS_model_object@params_table["j", "max"] <- 1.99999 # Max is typically 2 for DIVALIKE+J.
    } else if (model_name == "BAYAREALIKE") {
      run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
      run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
      run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
      run_object$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.0
    } else if (model_name == "BAYAREALIKE+J") {
      run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
      run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
      run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
      run_object$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.0
      run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
      run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- 0.01
      run_object$BioGeoBEARS_model_object@params_table["j", "max"] <- 0.99999 # Max is typically 1 for BAYAREALIKE+J.
    } else {
      warning("Model '", model_name, "' is not recognized. Skipping.")
      next
    }
    
    # Run the optimization.
    optim_error <- NULL
    res <- NULL
    # A final check on parameter bounds before running.
    run_object <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(run_object)
    tryCatch({
      res <- BioGeoBEARS::bears_optim_run(run_object)
    }, error = function(e) {
      optim_error <<- e
      print(paste("ERROR during bears_optim_run for model", model_name, ":", e$message))
    })
    
    # Store the results if the run was successful.
    if (!is.null(res)) {
      results_list[[model_name]] <- res
      print(paste("--- Finished model:", model_name, "---"))
    }
  }
  
  # --- 7. Summarize Results and Select the Best Model ---
  print("--- 7. Collecting results and selecting the best model ---")
  if (length(results_list) == 0) {
    stop("No models were successfully run. Halting execution.")
  }
  
  # Get the first valid result object to access shared inputs like the tree file.
  first_valid_res <- results_list[[1]]
  if (is.null(first_valid_res)) {
    stop("All model runs failed, cannot proceed to build summary table.")
  }
  
  tr <- ape::read.tree(first_valid_res$inputs$trfn)
  n_samples <- length(tr$tip.label)
  
  # Create an empty list to store a one-row data frame for each model's results.
  list_for_restable <- list()
  
  for (model_name in names(results_list)) {
    res <- results_list[[model_name]]
    
    lnL <- res$total_loglikelihood
    
    params_table <- res$inputs$BioGeoBEARS_model_object@params_table
    # Count the number of free parameters directly.
    num_params <- sum(params_table$type == "free")
    
    # Extract the final estimated values for d, e, and j from the 'est' column.
    d_est <- params_table["d", "est"]
    e_est <- params_table["e", "est"]
    # Check if 'j' is a free parameter in the current model; if so, extract its estimated value.
    j_est <- ifelse("j" %in% rownames(params_table) && params_table["j", "type"] == "free",
                    params_table["j", "est"],
                    NA)
    
    # Calculate AIC and AICc values.
    aic_val <- calc_AIC_column(LnL_vals = lnL, nparam_vals = num_params)
    aicc_val <- calc_AICc_column(LnL_vals = lnL, nparam_vals = num_params, samplesize = n_samples)
    
    # Construct a single-row data frame with all key information for this model.
    model_summary_row <- data.frame(
      model = model_name,
      LnL = lnL,
      num_params = num_params,
      d = d_est,
      e = e_est,
      j = j_est,
      AIC = aic_val$AIC,
      AICc = aicc_val$AICc,
      stringsAsFactors = FALSE
    )
    
    list_for_restable[[model_name]] <- model_summary_row
  }
  
  # Combine the list of data frames into a single results table.
  restable <- plyr::rbind.fill(list_for_restable)
  
  # Calculate AICc weights for model comparison.
  restable <- restable[order(restable$AICc), ] # Sort the table by AICc in ascending order.
  restable$delta_AICc <- restable$AICc - min(restable$AICc)
  restable$rel_likelihood <- exp(-0.5 * restable$delta_AICc)
  restable$AICc_weight <- restable$rel_likelihood / sum(restable$rel_likelihood)
  
  # Format the output for better readability.
  restable$LnL <- round(restable$LnL, 2)
  restable$AICc <- round(restable$AICc, 2)
  restable$delta_AICc <- round(restable$delta_AICc, 2)
  restable$AICc_weight <- round(restable$AICc_weight, 4)
  
  print("Model Comparison Results:")
  print(restable)
  
  # Save the model comparison table.
  comparison_table_path <- paste0(output_prefix, "_model_comparison.csv")
  write.csv(restable, file = comparison_table_path, row.names = FALSE)
  print(paste("Model comparison table saved to:", comparison_table_path))
  
  # Identify the best model based on the lowest AICc score.
  best_model_name <- restable$model[1]
  best_model_results <- results_list[[best_model_name]]
  print(paste("Best model selected based on AICc:", best_model_name))
  print("-------------------------------------------------")
  
  # Save all model results for future reference.
  all_results_rdata_path <- paste0(output_prefix, "_ALL_MODELS_results.rds")
  saveRDS(results_list, file = all_results_rdata_path)
  
  # --- 8. Save and Process Results for the Best Model ---
  print(paste("--- 8. Saving results for the best model (", best_model_name, ") ---"))
  
  # A. Save the raw result object for the best model.
  results_rdata_path <- paste0(output_prefix, "_BESTMODEL_", best_model_name, "_results.rds")
  print(paste("Saving raw results object to:", results_rdata_path))
  saveRDS(best_model_results, file = results_rdata_path)
  
  # B. Generate the main plot (PNG) using the best model's results.
  plot_png_path <- paste0(output_prefix, "_BESTMODEL_", best_model_name, "_plot.png")
  print(paste("Generating main plot PNG:", plot_png_path))
  
  tryCatch({
    
    # Dynamically set plotting dimensions based on the number of tips.
    num_tips <- n_samples
    width_px <- max(plot_width_px, 10 * num_tips)  # Allocate space for each tip label.
    height_px <- max(plot_height_px, 50 * num_tips) # Ensure visibility for tall trees.
    
    # Set a maximum pixel dimension to prevent memory issues.
    max_px <- 48000
    width_px <- min(ceiling(max_px / 5), width_px)
    height_px <- min(max_px, height_px)
    
    analysis_titletxt_plot <- paste(basename(tree_basename), "(Best Model:", best_model_name, ")")
    
    cat("Cleaning up graphics devices...\n")
    while (dev.cur() > 1) {
      dev.off()
    }
    
    png(filename = plot_png_path, width = width_px, height = height_px, res = plot_res_dpi)
    
    # Use layout() to create a two-panel plot: one for the tree, one for the legend.
    # Adjust the widths ratio as needed (e.g., c(4, 1)).
    layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
    
    # Panel 1: Plot the main tree.
    # Set margins for the first plot.
    par(mar = c(5, 4, 4, 0))
    BioGeoBEARS::plot_BioGeoBEARS_results(
      results_object = best_model_results,
      analysis_titletxt = analysis_titletxt_plot,
      plotwhat = "pie",
      plotlegend = FALSE, # CRITICAL: Do not plot the legend here.
      label.offset = 0.02,
      tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8,
      plotsplits = TRUE,
      show.tip.label = TRUE
    )
    
    # Panel 2: Plot only the legend.
    # Set smaller margins for a compact legend.
    par(mar = c(5, 0, 4, 1))
    BioGeoBEARS::plot_BioGeoBEARS_results(
      results_object = best_model_results,
      plotlegend = TRUE,   # CRITICAL: Plot the legend here.
      skiptree = TRUE,     # CRITICAL: Skip plotting the tree.
      legend_cex = 1.0     # Adjust legend text size if needed.
    )
    dev.off()
    print("Main plot PNG for best model saved.")
    
  }, error = function(e) {
    warning("Failed to generate main plot PNG: ", e$message)
    if (names(dev.cur()) != "null device") {
      try(dev.off(), silent = TRUE)
    }
  })
  
  # Call the next function to perform Biogeographical Stochastic Mapping (BSM).
  biogeobears_transition_matrices(
    best_model_results = best_model_results,
    save_biogeobears_path = save_biogeobears_path,
    comparison_table = restable,
    time_boundaries = time_boundaries_no_zero
  )
  
  # --- 9. Return a list containing all important information ---
  return(invisible(list(
    comparison_table = restable,
    all_results = results_list,
    best_result = best_model_results
  )))
  
}

# # --- Execute the full BioGeoBEARS pipeline ---
# result_info <- run_biogeobears_pipeline(
#   tree_filepath = tree_filepath,
#   # tip_states_filepath is now auto-generated inside the function
#   timeperiods_filepath = timeperiods_filepath,
#   dispersal_multipliers_filepath = dispersal_multipliers_filepath
# )

# =============================================================================
# create_chord_diagram - 创建圆环图展示区域间扩散流向
# =============================================================================
#
# 【目的】
# 创建弦图 (chord diagram) 展示区域间的扩散流向。
# 弦图直观显示哪些区域之间有扩散事件，以及扩散的相对强度。
#
# 【输入参数】
# - transition_matrix: 转换矩阵 (行 = 来源区, 列 = 目标区)
# - title:             图表标题
# - area_colors:       区域颜色 (可选，自动生成默认颜色)
#
# 【依赖】
# 需要 circlize 包。如果未安装会跳过并发出警告。
#
# =============================================================================

create_chord_diagram <- function(transition_matrix,
                                  title = "Biogeographic Dispersal Flow",
                                  area_colors = NULL,
                                  min_value = 0) {

  # Check circlize package
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("Package 'circlize' not installed. Skipping chord diagram.")
    return(invisible(NULL))
  }

  if (!is.matrix(transition_matrix)) {
    transition_matrix <- as.matrix(transition_matrix)
  }

  areanames <- rownames(transition_matrix)
  if (is.null(areanames)) {
    areanames <- colnames(transition_matrix)
  }
  if (is.null(areanames)) {
    areanames <- paste0("Area", 1:nrow(transition_matrix))
  }

  rownames(transition_matrix) <- areanames
  colnames(transition_matrix) <- areanames

  # Filter small values and handle NA
  transition_matrix[is.na(transition_matrix)] <- 0
  if (min_value > 0) {
    transition_matrix[transition_matrix < min_value] <- 0
  }

  # Skip if all zeros

  if (all(transition_matrix == 0)) {
    return(invisible(NULL))
  }

  # Set default colors for biogeographic realms
  if (is.null(area_colors)) {
    default_colors <- c(
      "N" = "#E41A1C", "T" = "#377EB8", "P" = "#4DAF4A", "A" = "#984EA3",
      "S" = "#FF7F00", "J" = "#FFFF33", "F" = "#A65628", "M" = "#F781BF",
      "I" = "#999999", "U" = "#66C2A5", "O" = "#8DA0CB"
    )
    area_colors <- sapply(areanames, function(a) {
      if (a %in% names(default_colors)) default_colors[a] else rainbow(length(areanames))[which(areanames == a)]
    })
    names(area_colors) <- areanames
  }

  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 4)

  tryCatch({
    circlize::chordDiagram(
      transition_matrix,
      grid.col = area_colors,
      transparency = 0.3,
      directional = 1,
      direction.type = c("diffHeight", "arrows"),
      link.arr.type = "big.arrow",
      annotationTrack = c("name", "grid")
    )
    title(main = title, cex.main = 1.0, line = -1)
  }, error = function(e) {
    # Silently handle errors
  })

  circlize::circos.clear()

  return(invisible(transition_matrix))
}


# =============================================================================
# generate_chord_diagrams_set - 生成一组圆环图（全局 + 时间切片）
# =============================================================================
#
# 生成两张图：
# 1. 全局圆环图 (单图)
# 2. 时间切片圆环图 (多个小图组合)
#
# =============================================================================

generate_chord_diagrams_set <- function(global_matrix,
                                         timeslice_list,
                                         metric_type,
                                         savedir,
                                         model_name) {

  if (!requireNamespace("circlize", quietly = TRUE)) {
    cat("  Package 'circlize' not installed. Skipping chord diagrams.\n")
    return(invisible(NULL))
  }

  cat(paste0("\n--- Generating ", metric_type, " chord diagrams ---\n"))

  # --- 1. Global Chord Diagram ---
  global_path <- file.path(savedir, paste0(model_name, "_chord_", metric_type, "_GLOBAL.png"))

  png(filename = global_path, width = 10, height = 10, units = "in", res = 300)
  create_chord_diagram(
    transition_matrix = global_matrix,
    title = paste0("Total Dispersal (d+j) - ", metric_type, " (Global)")
  )
  dev.off()
  cat(paste0("  Global ", metric_type, " chord diagram saved: ", basename(global_path), "\n"))

  # --- 2. Time-Stratified Chord Diagrams (combined into one figure) ---
  n_slices <- length(timeslice_list)
  if (n_slices == 0) {
    cat("  No time slices to plot.\n")
    return(invisible(NULL))
  }

  # Calculate grid layout
  n_cols <- min(3, n_slices)
  n_rows <- ceiling(n_slices / n_cols)

  timeslice_path <- file.path(savedir, paste0(model_name, "_chord_", metric_type, "_TIMESLICES.png"))

  # Set up PNG with appropriate dimensions
  png(filename = timeslice_path,
      width = 6 * n_cols,
      height = 6 * n_rows,
      units = "in",
      res = 200)

  par(mfrow = c(n_rows, n_cols), mar = c(1, 1, 3, 1))

  for (i in seq_along(timeslice_list)) {
    slice_name <- names(timeslice_list)[i]
    slice_matrix <- timeslice_list[[i]]

    # Handle NA matrices (time slices outside tree age)
    if (all(is.na(slice_matrix)) || all(slice_matrix == 0, na.rm = TRUE)) {
      plot.new()
      title(main = paste0(slice_name, "\n(No data)"), cex.main = 0.9)
    } else {
      create_chord_diagram(
        transition_matrix = slice_matrix,
        title = slice_name
      )
    }
  }

  dev.off()
  cat(paste0("  Time-stratified ", metric_type, " chord diagrams saved: ", basename(timeslice_path), "\n"))

  return(invisible(NULL))
}


biogeobears_transition_matrices <- function(best_model_results,
                                            save_biogeobears_path,
                                            comparison_table,
                                            bsm_sims = 100,
                                            time_boundaries) {
  
  # This first line re-defines the save path based on a global variable 'tree_filepath', 
  # which could be problematic. It should ideally use the passed 'save_biogeobears_path'.
  # save_biogeobears_path <- paste0(tools::file_path_sans_ext(tree_filepath), "_biogeobears_results")
  # dir.create(save_biogeobears_path, showWarnings = FALSE)
  
  # --- 7. BioGeoBEARS Biogeographical Stochastic Mapping (BSM) ---
  print("--- 7. Starting Biogeographical Stochastic Mapping (BSM) ---")
  
  # Assign the results object to 'res', the conventional name for BSM function inputs.
  res <- best_model_results
  tree_fn_temp <- res$inputs$trfn
  geog_fn_temp <- res$inputs$geogfn
  
  # Verify that the temporary tree and geography files exist.
  if (!exists("tree_fn_temp") || !file.exists(tree_fn_temp)) {
    stop("Error: Temporary tree file 'tree_fn_temp' not found. BSM requires this file.")
  }
  if (!exists("geog_fn_temp") || !file.exists(geog_fn_temp)) {
    stop("Error: Temporary geography data file 'geog_fn_temp' not found. BSM requires this file.")
  }
  
  tr <- ape::read.tree(tree_fn_temp) # Read the processed tree.
  tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = geog_fn_temp)
  areanames <- colnames(tipranges@df)
  
  print("Area names used in this analysis (areanames):")
  print(areanames)
  
  # Check if area names are single characters, as some BSM functions expect this format.
  if (any(nchar(areanames) > 1)) {
    warning("Warning: Some BSM plotting and summary functions may expect single-character area names. This could cause issues or aesthetic problems in downstream steps.")
  }
  
  # Get the name of the best model for file naming conventions.
  model_name <- comparison_table[1, "model"]
  
  # --- 7.1 Get Inputs for BSM ---
  # This pre-calculates branch likelihoods, which can be slow for large trees.
  BSM_inputs_fn <- file.path(save_biogeobears_path, paste0(model_name, "_BSM_inputs.Rdata"))
  cat("Calculating BSM inputs...\n")
  
  # Ensure the 'res' object has the num_cores_to_use parameter.
  if (is.null(res$inputs$num_cores_to_use)) {
    res$inputs$num_cores_to_use <- 1 
  }
  stochastic_mapping_inputs_list <- BioGeoBEARS::get_inputs_for_stochastic_mapping(res = res)
  save(stochastic_mapping_inputs_list, file = BSM_inputs_fn)
  cat("BSM inputs calculated and saved to:", BSM_inputs_fn, "\n")
  
  # --- 7.2 Run BSM (Generate Stochastic Maps) ---
  # Define BSM parameters.
  nummaps_goal <- bsm_sims # The desired number of stochastic maps.
  maxnum_maps_to_try <- nummaps_goal * 2 # Max attempts to reach the goal.
  maxtries_per_branch <- 40000
  save_after_every_try <- TRUE # Recommended for long runs.
  savedir <- save_biogeobears_path # Directory to save intermediate results.
  seedval <- 1 # For reproducibility, this should be a fixed integer.
  
  # Define output filenames.
  clado_events_fn <- paste0(savedir, "/", model_name, "_clado_events_tables.rds")
  ana_events_fn <- paste0(savedir, "/", model_name, "_ana_events_tables.rds")
  
  runBSMslow <- TRUE # Set to FALSE to load previously saved BSM results.
  
  if (runBSMslow || !file.exists(clado_events_fn) || !file.exists(ana_events_fn)) {
    cat("Running BSM to generate", nummaps_goal, "stochastic maps...\n")
    set.seed(seedval) # Set random seed for reproducibility.
    BSM_output <- BioGeoBEARS::runBSM(res = res,
                                      stochastic_mapping_inputs_list = stochastic_mapping_inputs_list,
                                      maxnum_maps_to_try = maxnum_maps_to_try,
                                      nummaps_goal = nummaps_goal,
                                      maxtries_per_branch = maxtries_per_branch,
                                      save_after_every_try = save_after_every_try,
                                      savedir = savedir,
                                      seedval = seedval,
                                      wait_before_save = 0.01,
                                      master_nodenum_toPrint = 0) # Set to 0 to reduce console output.
    
    # The runBSM function returns an object containing the event tables.
    RES_clado_events_tables <- BSM_output$RES_clado_events_tables
    RES_ana_events_tables <- BSM_output$RES_ana_events_tables
    
    # Save these tables to the specified filenames.
    saveRDS(RES_clado_events_tables, file = clado_events_fn)
    saveRDS(RES_ana_events_tables, file = ana_events_fn)
    
    cat("BSM complete. Cladogenetic event tables saved to:", clado_events_fn, "\n")
    cat("Anagenetic event tables saved to:", ana_events_fn, "\n")
  } else {
    cat("Loading previously saved BSM event tables...\n")
    RES_clado_events_tables <- readRDS(file = clado_events_fn)
    RES_ana_events_tables <- readRDS(file = ana_events_fn)
    # Create the BSM_output list for compatibility with subsequent code.
    BSM_output <- list(RES_clado_events_tables = RES_clado_events_tables,
                       RES_ana_events_tables = RES_ana_events_tables)
  }
  
  # --- 7.3 Summarize BSM Events and Calculate Transition Counts ---
  # Extract the event tables from the BSM output.
  clado_events_tables <- BSM_output$RES_clado_events_tables
  ana_events_tables <- BSM_output$RES_ana_events_tables
  
  # (Optional but recommended) Simulate source areas for dispersal events to make 'from-to' transitions explicit.
  cat("Simulating source areas for dispersal events (this may take some time)...\n")
  BSMs_w_sourceAreas <- BioGeoBEARS::simulate_source_areas_ana_clado(
    res = res,
    clado_events_tables = clado_events_tables,
    ana_events_tables = ana_events_tables,
    areanames = areanames
  )
  
  clado_events_tables_src <- BSMs_w_sourceAreas$clado_events_tables
  ana_events_tables_src <- BSMs_w_sourceAreas$ana_events_tables
  
  # Count all anagenetic and cladogenetic events across the stochastic maps.
  cat("Counting biogeographical events...\n")
  counts_list <- BioGeoBEARS::count_ana_clado_events(clado_events_tables = clado_events_tables_src,
                                                     ana_events_tables = ana_events_tables_src,
                                                     areanames = areanames,
                                                     actual_names = areanames)
  
  # Print a summary of event counts (averaged over all maps).
  cat("\n--- BSM Event Count Summary (mean per BSM) ---\n")
  print(BioGeoBEARS::conditional_format_table(counts_list$summary_counts_BSMs))
  
  cat("\n--- Mean Range-Expansion Events (d-type, from-to) ---\n")
  if (!is.null(counts_list$d_counts_fromto_means)) {
    print(BioGeoBEARS::conditional_format_table(counts_list$d_counts_fromto_means))
  } else {
    print("No d-type events were detected, or the from-to matrix could not be calculated.")
  }
  
  cat("\n--- Mean Founder/Jump Dispersal Events (j-type, from-to) ---\n")
  if (!is.null(counts_list$founder_counts_fromto_means)) {
    print(BioGeoBEARS::conditional_format_table(counts_list$founder_counts_fromto_means))
  } else {
    print("No j-type events were detected, or the from-to matrix could not be calculated.")
  }
  
  # --- 7.4 Calculate Transition Rates ---
  # Total branch length of the phylogenetic tree (units should be in Myr).
  total_branch_length_Myrs <- sum(tr$edge.length)
  cat("\nTotal tree branch length (units match the input tree, typically Myr):", total_branch_length_Myrs, "\n")
  
  # Initialize rate matrices as NULL in case no events of a certain type occurred.
  rates_d_events <- NULL
  rates_j_events <- NULL
  
  # Calculate transition rate matrix for d-type (range expansion) events.
  if (!is.null(counts_list$d_counts_fromto_means) && total_branch_length_Myrs > 0) {
    rates_d_events <- counts_list$d_counts_fromto_means / total_branch_length_Myrs
    cat("\n--- Range-Expansion (d-type) Event Transition Rates (events/Myr, from-to) ---\n")
    print(BioGeoBEARS::conditional_format_table(rates_d_events))
    # Save to file.
    write.table(BioGeoBEARS::conditional_format_table(rates_d_events),
                file = paste0(savedir, "/", model_name, "_transition_rates_d_events.txt"),
                quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
    cat("d-type event transition rates saved to:", paste0(model_name, "_transition_rates_d_events.txt"), "\n")
    
    if (all(rates_d_events == 0)) {
      message("NOTE: All d-type event rates are zero. Skipping heatmap plot.")
    } else {
      rates_mat <- as.matrix(rates_d_events)
      # Sort matrix rows and columns alphabetically for consistent comparison
      sorted_names <- sort(rownames(rates_mat))
      rates_mat <- rates_mat[sorted_names, sorted_names]
      numbers_sci <- matrix(formatC(rates_mat, format = "e", digits = 2),
                            nrow = nrow(rates_mat), ncol = ncol(rates_mat), dimnames = dimnames(rates_mat))
      pheatmap::pheatmap(
        rates_mat,
        display_numbers = numbers_sci,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = "Range-Expansion (d-type) Transition Rate Heatmap (Global Average)\n(Y-axis: From -> X-axis: To)",
        fontsize_main = 12,
        filename = paste0(savedir, "/", model_name, "_transition_rates_d_events.png"),
        width = 8, height = 7, silent = TRUE
      )
    }
  } else {
    cat("\nCould not calculate d-type event transition rates (no events or zero total branch length).\n")
  }
  
  # Calculate transition rate matrix for j-type (founder/jump) events.
  if (!is.null(counts_list$founder_counts_fromto_means) && total_branch_length_Myrs > 0) {
    rates_j_events <- counts_list$founder_counts_fromto_means / total_branch_length_Myrs
    cat("\n--- Founder/Jump (j-type) Event Transition Rates (events/Myr, from-to) ---\n")
    print(BioGeoBEARS::conditional_format_table(rates_j_events))
    # Save to file.
    write.table(BioGeoBEARS::conditional_format_table(rates_j_events),
                file = paste0(savedir, "/", model_name, "_transition_rates_j_events.txt"),
                quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
    cat("j-type event transition rates saved to:", paste0(model_name, "_transition_rates_j_events.txt"), "\n")
    
    if (all(rates_j_events == 0)) {
      message("NOTE: All j-type event rates are zero. Skipping heatmap plot.")
    } else {
      rates_mat <- as.matrix(rates_j_events)
      # Sort matrix rows and columns alphabetically for consistent comparison
      sorted_names <- sort(rownames(rates_mat))
      rates_mat <- rates_mat[sorted_names, sorted_names]
      numbers_sci <- matrix(formatC(rates_mat, format = "e", digits = 2),
                            nrow = nrow(rates_mat), ncol = ncol(rates_mat), dimnames = dimnames(rates_mat))
      pheatmap::pheatmap(
        rates_mat,
        display_numbers = numbers_sci,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = "Founder-Event (j-type) Transition Rate Heatmap (Global Average)\n(Y-axis: From -> X-axis: To)",
        fontsize_main = 12,
        filename = paste0(savedir, "/", model_name, "_transition_rates_j_events.png"),
        width = 8, height = 7, silent = TRUE
      )
    }
  } else {
    cat("\nCould not calculate j-type event transition rates (no events or zero total branch length).\n")
  }
  
  # --- Calculate and visualize total dispersal rates (d+j) ---
  cat("\n--- Calculating total dispersal rates (d+j) ---\n")
  # Ensure at least one of the rate matrices was successfully calculated.
  # Treat a NULL matrix (due to no events) as a zero matrix for addition.
  if (!is.null(rates_d_events) || !is.null(rates_j_events)) {
    
    # Create a zero matrix as a template.
    template_matrix <- if (!is.null(rates_d_events)) rates_d_events else rates_j_events
    total_dispersal_rates <- matrix(0,
                                    nrow = nrow(template_matrix),
                                    ncol = ncol(template_matrix),
                                    dimnames = dimnames(template_matrix))
    # Add d and j rates to the total matrix.
    if (!is.null(rates_d_events)) {
      total_dispersal_rates <- total_dispersal_rates + rates_d_events
    }
    if (!is.null(rates_j_events)) {
      total_dispersal_rates <- total_dispersal_rates + rates_j_events
    }
    
    # Check if the sum matrix is all zeros.
    if (all(total_dispersal_rates == 0)) {
      message("NOTE: Total dispersal rates (d+j) are all zero. Skipping heatmap plot.")
    } else {
      rates_mat <- as.matrix(total_dispersal_rates)
      # Sort matrix rows and columns alphabetically for consistent comparison
      sorted_names <- sort(rownames(rates_mat))
      rates_mat <- rates_mat[sorted_names, sorted_names]
      numbers_sci <- matrix(formatC(rates_mat, format = "e", digits = 2),
                            nrow = nrow(rates_mat), ncol = ncol(rates_mat), dimnames = dimnames(rates_mat))
      pheatmap::pheatmap(
        rates_mat,
        display_numbers = numbers_sci,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = "Total Dispersal (d+j) Transition Rate Heatmap (Global Average)\n(Y-axis: From -> X-axis: To)",
        fontsize = 10,
        filename = paste0(savedir, "/", model_name, "_transition_rates_total_dispersal(d+j).png"),
        width = 8, height = 7, silent = TRUE
      )
      cat("Total dispersal (d+j) rate heatmap saved to:", paste0(savedir, "/", model_name, "_transition_rates_total_dispersal.png"), "\n")
    }
  } else {
    cat("\nNeither d-type nor j-type rates could be calculated; skipping total dispersal calculation.\n")
  }
  
  # --- Calculate per-area local extinction (e-type) rates ---
  
  # Check if the per-area event count table exists.
  if (!is.null(counts_list$e_counts_rectangle) && total_branch_length_Myrs > 0) {
    
    # The `$e_counts_rectangle` is a matrix of (num_sims x num_areas).
    # We calculate the mean of each column to get the average event count per area.
    mean_e_events_per_area <- colMeans(counts_list$e_counts_rectangle)
    
    # Divide the mean event count by total branch length to get the rate for each area.
    rates_e_per_area <- mean_e_events_per_area / total_branch_length_Myrs
    names(rates_e_per_area) <- areanames
    
    # Print the result, which is now a named vector of per-area rates.
    cat("\n--- Mean per-area local extinction (e-type) rates (events/Myr) ---\n")
    print(rates_e_per_area)
    
    # Convert the named vector to a two-column data frame for saving.
    e_rates_df <- data.frame(
      Area = names(rates_e_per_area),
      Extinction_Rate = rates_e_per_area,
      row.names = NULL 
    )
    
    # Save the data frame using write.table.
    write.table(
      e_rates_df,
      file = paste0(savedir, "/", model_name, "_extinction_rates_per_area.txt"),
      quote = FALSE,      # Do not add quotes around text.
      sep = "\t",         # Use tab as the separator.
      row.names = FALSE   # Do not write data frame row numbers.
    )
    
    cat("Per-area e-type event rates saved to:", paste0(savedir, "/", model_name, "_extinction_rates_per_area.txt"), "\n")
    
  } else {
    cat("\nCould not calculate per-area e-type event rates.\n")
  }
  
  # --- (Optional) Recalculate and print the global total rate for validation ---
  # This total rate should equal the sum of the per-area rates calculated above.
  if (!is.null(counts_list$e_totals_list) && total_branch_length_Myrs > 0) {
    rate_e_events_total <- mean(counts_list$e_totals_list) / total_branch_length_Myrs
    cat("\n--- (Validation) Global mean total local extinction (e-type) rate ---\n")
    print(rate_e_events_total)
    cat("(Note: This value should equal the sum of the per-area rates above)\n")
  }
  
  cat("\n--- BSM analysis and transition rate calculation complete. ---\n")
  
  # --- 7.5 (Optional) Further Visualization and Diagnostics ---
  
  # Plot histograms of event counts to assess uncertainty across simulations.
  hist_event_counts_fn <- paste0(savedir, "/", model_name, "_histograms_of_event_counts.pdf")
  BioGeoBEARS::hist_event_counts(
    counts_list = counts_list,
    pdffn = hist_event_counts_fn # Use the 'pdffn' argument to specify the filename.
  )
  # Interpretation Guide:
  # A narrow histogram indicates a robust conclusion.
  # A wide histogram indicates high uncertainty.
  
  # The frequency distributions for the total number of key biogeographical events, 
  # as inferred across 100 Bayesian Stochastic Maps (BSM), are shown in Figure X. 
  # Each panel displays a histogram for a specific event type (founder-event, vicariance, subset sympatry, narrow sympatry, and anagenetic dispersal). 
  # The x-axis indicates the total count of an event per simulated history, and the y-axis represents the number of simulations (frequency) that inferred that count. 
  # The shape of each distribution reveals the uncertainty in the frequency of that event; 
  # narrow peaks indicate high consistency and low uncertainty in the inferred number of events, whereas wide distributions suggest high uncertainty.
  
  cat("Event count histograms saved to:", hist_event_counts_fn, "\n")
  
  # Check consistency between ML ancestral state probabilities and BSM averages.
  # Requires the 'MultinomialCI' package.
  library(MultinomialCI)
  check_ML_vs_BSM_fn <- paste0(savedir, "/", sanitize_filename(model_name), "_ML_vs_BSM_check.pdf")

  # Ensure that the cladogenetic event tables list contains valid data frames.
  if (length(clado_events_tables_src) > 0 && inherits(clado_events_tables_src[[1]], "data.frame")) {
    check_ML_vs_BSM(res = res,
                    clado_events_tables = clado_events_tables_src,
                    model_name = sanitize_filename(model_name),
                    tr = tr,
                    plot_each_node = FALSE,
                    linreg_plot = TRUE,
                    MultinomialCI = TRUE)

  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", main = "Could not generate ML vs BSM plot:\n'clado_events_tables_src' is invalid")
    text(1, 1, "'clado_events_tables_src' may be empty or incorrectly formatted.")
    warning("Could not generate ML vs BSM plot; 'clado_events_tables_src' may be empty or incorrectly formatted.")
  }
  #cat("Comparison plot of ML ancestral states vs. BSM averages saved to:", check_ML_vs_BSM_fn, "\n")
  
  # Construct the full path to the original file, which was generated in the current working directory.
  source_file <- file.path(getwd(), paste0(sanitize_filename(model_name), "_ML_vs_BSM.pdf"))
  
  # Define the final destination path and filename where you want to store the plot.
  # (This is the same as the 'check_ML_vs_BSM_fn' variable from your previous code).
  destination_file <- check_ML_vs_BSM_fn
  
  # Before trying to move the file, it's best practice to check if it was successfully created.
  if (file.exists(source_file)) {
    # Execute the move and rename operation.
    file.rename(from = source_file, to = destination_file)
    # Print a success message, which is now accurate, confirming the final location.
    cat("Comparison plot of ML ancestral states vs. BSM averages moved and saved to:", destination_file, "\n")
  } else {
    # If the source file wasn't found, issue a warning.
    warning("BioGeoBEARS::check_ML_vs_BSM did not seem to create the expected output file at: ", source_file)
  }
  
  # --- Time-Stratified Event Rate Calculation ---

  # Create a list of time periods with upper and lower bounds.
  my_time_periods <- lapply(seq_along(time_boundaries), function(i) {
    upper_bound <- time_boundaries[i]
    lower_bound <- if (i == 1) 0 else time_boundaries[i - 1]
    return(c(upper_bound, lower_bound))
  })

  # --- NEW: Calculate total branch length within each time slice ---
  # This is needed to compute per-lineage rate (as opposed to flux)
  cat("\nCalculating branch lengths per time slice...\n")

  # Get node heights (distance from root) for all nodes
  node_heights <- ape::nodeHeights(tr)
  # Convert to "time before present" (age)
  tree_height <- max(node_heights)
  # node_ages: column 1 = parent age, column 2 = child age
  node_ages <- tree_height - node_heights

  # Calculate total branch length within each time slice
  branch_length_per_slice <- numeric(length(my_time_periods))
  names(branch_length_per_slice) <- sapply(my_time_periods, function(p) paste(p[1], "-", p[2], "Ma"))

  for (ts_index in 1:length(my_time_periods)) {
    slice_upper <- my_time_periods[[ts_index]][1]  # older bound
    slice_lower <- my_time_periods[[ts_index]][2]  # younger bound

    total_length_in_slice <- 0

    # Iterate through each edge
    for (edge_idx in 1:nrow(tr$edge)) {
      parent_age <- node_ages[edge_idx, 1]
      child_age <- node_ages[edge_idx, 2]

      # Calculate overlap between this edge and the time slice
      # Edge spans from parent_age (older) to child_age (younger)
      overlap_start <- min(parent_age, slice_upper)
      overlap_end <- max(child_age, slice_lower)

      overlap_length <- overlap_start - overlap_end
      if (overlap_length > 0) {
        total_length_in_slice <- total_length_in_slice + overlap_length
      }
    }

    branch_length_per_slice[ts_index] <- total_length_in_slice
  }

  cat("Branch lengths per time slice (Ma):\n")
  print(round(branch_length_per_slice, 2))
  
  # Initialize lists to store event counts for each time slice, separately for d and j events.
  d_timeslice_counts_list <- list()
  j_timeslice_counts_list <- list()
  
  for (i in 1:length(my_time_periods)) {
    zero_matrix <- matrix(0, nrow = length(areanames), ncol = length(areanames), dimnames = list(areanames, areanames))
    d_timeslice_counts_list[[i]] <- zero_matrix
    j_timeslice_counts_list[[i]] <- zero_matrix
  }
  
  names(d_timeslice_counts_list) <- sapply(my_time_periods, function(p) paste(p[1], "-", p[2], "Ma"))
  names(j_timeslice_counts_list) <- sapply(my_time_periods, function(p) paste(p[1], "-", p[2], "Ma"))
  
  # Iterate through each BSM simulation to tally events within each time slice.
  for (i in 1:length(ana_events_tables_src)) {
    
    # Tally d-type (anagenetic) events.
    events_table <- ana_events_tables_src[[i]]
    if (is.data.frame(events_table) && nrow(events_table) > 0) {
      d_events <- events_table[which(events_table$event_type == "d"), ]
      if (nrow(d_events) > 0) {
        for (k in 1:nrow(d_events)) {
          event <- d_events[k, ]
          from_area <- event$ana_dispersal_from
          to_area <- event$dispersal_to
          event_time <- event$abs_event_time
          
          if (!is.na(from_area) && !is.na(to_area) && from_area %in% areanames && to_area %in% areanames) {
            for (ts_index in 1:length(my_time_periods)) {
              time_slice <- my_time_periods[[ts_index]]
              if (event_time >= time_slice[2] && event_time < time_slice[1]) {
                d_timeslice_counts_list[[ts_index]][from_area, to_area] <- d_timeslice_counts_list[[ts_index]][from_area, to_area] + 1
                break
              }
            }
          }
        }
      }
    }
    
    # Tally j-type (cladogenetic) events.
    clado_table <- clado_events_tables_src[[i]]
    if (is.data.frame(clado_table) && nrow(clado_table) > 0) {
      j_events <- clado_table[which(clado_table$clado_event_type == "founder (j)"), ]
      if (nrow(j_events) > 0) {
        for (k in 1:nrow(j_events)) {
          event <- j_events[k, ]
          from_area <- event$clado_dispersal_from
          to_area <- event$clado_dispersal_to
          event_time <- event$time_bp
          
          if (!is.na(from_area) && !is.na(to_area) && from_area %in% areanames && to_area %in% areanames) {
            for (ts_index in 1:length(my_time_periods)) {
              time_slice <- my_time_periods[[ts_index]]
              if (event_time >= time_slice[2] && event_time < time_slice[1]) {
                j_timeslice_counts_list[[ts_index]][from_area, to_area] <- j_timeslice_counts_list[[ts_index]][from_area, to_area] + 1
                break
              }
            }
          }
        }
      }
    }
  }
  
  # Calculate the average number of events per time slice.
  avg_d_timeslice_list <- lapply(d_timeslice_counts_list, function(mat) mat / bsm_sims)
  avg_j_timeslice_list <- lapply(j_timeslice_counts_list, function(mat) mat / bsm_sims)
  # Calculate the total average by summing the d and j matrices.
  avg_total_timeslice_list <- mapply("+", avg_d_timeslice_list, avg_j_timeslice_list, SIMPLIFY = FALSE)

  # --- Calculate FLUX (events per unit time) ---
  # Flux = Count / Duration
  # This measures "how busy" a time period was (absolute activity level)
  durations_Ma <- sapply(my_time_periods, function(p) p[1] - p[2])

  avg_d_timeslice_flux_list <- mapply("/", avg_d_timeslice_list, durations_Ma, SIMPLIFY = FALSE)
  avg_j_timeslice_flux_list <- mapply("/", avg_j_timeslice_list, durations_Ma, SIMPLIFY = FALSE)
  avg_total_timeslice_flux_list <- mapply("/", avg_total_timeslice_list, durations_Ma, SIMPLIFY = FALSE)

  # --- Calculate RATE (events per unit branch length) ---
  # Rate = Count / Total Branch Length in Slice
  # This measures "per-lineage dispersal propensity" (intrinsic rate)

  # Handle zero branch lengths to avoid Inf/NaN
  # Replace zeros with NA for safe division, then convert Inf/NaN to NA in results
  safe_branch_length <- branch_length_per_slice
  zero_slice_indices <- which(safe_branch_length == 0 | is.na(safe_branch_length))

  if (length(zero_slice_indices) > 0) {
    zero_slice_names <- names(branch_length_per_slice)[zero_slice_indices]
    cat("\nWARNING: The following time slices have zero branch length (outside tree age or no lineages):\n")
    cat("  ", paste(zero_slice_names, collapse = ", "), "\n")
    cat("  Rate values for these slices will be set to NA.\n\n")
  }

  # Custom safe division function that returns NA instead of Inf/NaN
  safe_divide_matrix <- function(mat, divisor) {
    if (is.na(divisor) || divisor == 0) {
      # Return a matrix of NAs with the same dimensions
      result <- mat
      result[,] <- NA
      return(result)
    } else {
      return(mat / divisor)
    }
  }

  avg_d_timeslice_rate_list <- mapply(safe_divide_matrix, avg_d_timeslice_list, branch_length_per_slice, SIMPLIFY = FALSE)
  avg_j_timeslice_rate_list <- mapply(safe_divide_matrix, avg_j_timeslice_list, branch_length_per_slice, SIMPLIFY = FALSE)
  avg_total_timeslice_rate_list <- mapply(safe_divide_matrix, avg_total_timeslice_list, branch_length_per_slice, SIMPLIFY = FALSE)

  # --- Save results ---
  # Save FLUX results
  timeslice_flux_results <- list(
    d_flux = avg_d_timeslice_flux_list,
    j_flux = avg_j_timeslice_flux_list,
    total_flux = avg_total_timeslice_flux_list,
    durations_Ma = durations_Ma,
    time_periods = my_time_periods
  )
  flux_fn <- file.path(savedir, paste0(model_name, "_avg_timeslice_FLUX_results.rds"))
  saveRDS(timeslice_flux_results, flux_fn)
  cat("\nFlux results saved to:", flux_fn, "\n")

  # Save RATE results
  timeslice_rate_results <- list(
    d_rate = avg_d_timeslice_rate_list,
    j_rate = avg_j_timeslice_rate_list,
    total_rate = avg_total_timeslice_rate_list,
    branch_length_per_slice = branch_length_per_slice,
    time_periods = my_time_periods
  )
  rate_fn <- file.path(savedir, paste0(model_name, "_avg_timeslice_RATE_results.rds"))
  saveRDS(timeslice_rate_results, rate_fn)
  cat("Rate results saved to:", rate_fn, "\n")

  # --- Print FLUX results ---
  cat("\n\n")
  cat("###########################################################################\n")
  cat("#                    TIME-STRATIFIED FLUX RESULTS                        #\n")
  cat("#        (events per Myr - measures absolute activity level)             #\n")
  cat("###########################################################################\n\n")

  cat("--- (1) D-TYPE Flux per Time Slice (events/Myr) ---\n")
  print(lapply(avg_d_timeslice_flux_list, function(mat) round(mat, 4)))

  cat("\n--- (2) J-TYPE Flux per Time Slice (events/Myr) ---\n")
  print(lapply(avg_j_timeslice_flux_list, function(mat) round(mat, 4)))

  cat("\n--- (3) TOTAL (d+j) Flux per Time Slice (events/Myr) ---\n")
  print(lapply(avg_total_timeslice_flux_list, function(mat) round(mat, 4)))

  # --- Print RATE results ---
  cat("\n\n")
  cat("###########################################################################\n")
  cat("#                    TIME-STRATIFIED RATE RESULTS                        #\n")
  cat("#     (events per branch-length-Myr - measures per-lineage propensity)   #\n")
  cat("###########################################################################\n\n")

  cat("--- (1) D-TYPE Rate per Time Slice (events/branch-Myr) ---\n")
  print(lapply(avg_d_timeslice_rate_list, function(mat) round(mat, 6)))

  cat("\n--- (2) J-TYPE Rate per Time Slice (events/branch-Myr) ---\n")
  print(lapply(avg_j_timeslice_rate_list, function(mat) round(mat, 6)))

  cat("\n--- (3) TOTAL (d+j) Rate per Time Slice (events/branch-Myr) ---\n")
  print(lapply(avg_total_timeslice_rate_list, function(mat) round(mat, 6)))
  
  # --- Visualize All Time-Stratified FLUX and RATE Matrices ---

  # A. Global Setup for Plotting
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(grid)

  # B. A more flexible and reusable plotting function for heatmaps
  # Added 'legend_label' parameter to distinguish Flux vs Rate
  plot_heatmap <- function(mat, subplot_title, scale_mode = "global", color_limits = NULL,
                           legend_label = "Value") {

    # Sort matrix rows and columns alphabetically for consistent comparison
    sorted_names <- sort(rownames(mat))
    mat <- mat[sorted_names, sorted_names]

    current_limits <- NULL
    if (scale_mode == "local") {
      local_max <- if (length(mat) > 0 && max(mat, na.rm = TRUE) > 0) max(mat, na.rm = TRUE) else 1
      current_limits <- c(0, local_max)
    } else {
      current_limits <- color_limits
    }

    df <- reshape2::melt(as.matrix(mat))
    colnames(df) <- c("From", "To", "Value")

    ggplot(df, aes(x = To, y = From, fill = Value)) +
      geom_tile(color = "grey80") +
      geom_text(aes(label = ifelse(Value > 0, sprintf("%.1e", Value), "0.0e+00")), size = 3) +
      scale_fill_gradient(
        low = "lightblue", high = "darkred",
        limits = current_limits,
        name = legend_label
      ) +
      ggtitle(subplot_title) +
      labs(x = "To (Destination)", y = "From (Source)") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10, face = "bold")
      )
  }

  # Helper function to generate a complete set of heatmaps for one metric type
  generate_heatmap_set <- function(d_list, j_list, total_list, metric_type, legend_label, savedir, model_name) {

    cat(paste0("\n--- Generating ", metric_type, " heatmaps ---\n"))

    # Calculate global color scale
    all_values <- unlist(c(d_list, j_list))
    val_max <- if (length(all_values) > 0 && max(all_values, na.rm = TRUE) > 0) max(all_values, na.rm = TRUE) else 1
    global_color_limits <- c(0, val_max)

    # --- D-TYPE ---
    titles_d <- names(d_list)

    # Global Scale
    plots_d_global <- mapply(plot_heatmap, d_list, titles_d,
                             MoreArgs = list(scale_mode = "global", color_limits = global_color_limits, legend_label = legend_label),
                             SIMPLIFY = FALSE)
    plot_grid_d_global <- gridExtra::grid.arrange(
      grobs = plots_d_global, ncol = 3,
      top = grid::textGrob(paste0("Time-Stratified D-type ", metric_type, " (Global Scale)"),
                           gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    save_path <- file.path(savedir, paste0(model_name, "_Time-Stratified_", metric_type, "_D_events_GLOBAL.png"))
    ggsave(filename = save_path, plot = plot_grid_d_global, width = 20, height = 5 * ceiling(length(plots_d_global)/3), dpi = 300, limitsize = FALSE)
    cat(paste0("  D-type ", metric_type, " (Global) saved: ", basename(save_path), "\n"))

    # Independent Scale
    plots_d_local <- mapply(plot_heatmap, d_list, titles_d,
                            MoreArgs = list(scale_mode = "local", legend_label = legend_label),
                            SIMPLIFY = FALSE)
    plot_grid_d_local <- gridExtra::grid.arrange(
      grobs = plots_d_local, ncol = 3,
      top = grid::textGrob(paste0("Time-Stratified D-type ", metric_type, " (Independent Scales)"),
                           gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    save_path <- file.path(savedir, paste0(model_name, "_Time-Stratified_", metric_type, "_D_events_INDEPENDENT.png"))
    ggsave(filename = save_path, plot = plot_grid_d_local, width = 20, height = 5 * ceiling(length(plots_d_local)/3), dpi = 300, limitsize = FALSE)
    cat(paste0("  D-type ", metric_type, " (Independent) saved: ", basename(save_path), "\n"))

    # --- J-TYPE ---
    if (max(unlist(j_list), na.rm = TRUE) > 0) {
      titles_j <- names(j_list)

      # Global Scale
      plots_j_global <- mapply(plot_heatmap, j_list, titles_j,
                               MoreArgs = list(scale_mode = "global", color_limits = global_color_limits, legend_label = legend_label),
                               SIMPLIFY = FALSE)
      plot_grid_j_global <- gridExtra::grid.arrange(
        grobs = plots_j_global, ncol = 3,
        top = grid::textGrob(paste0("Time-Stratified J-type ", metric_type, " (Global Scale)"),
                             gp = grid::gpar(fontsize = 16, fontface = "bold"))
      )
      save_path <- file.path(savedir, paste0(model_name, "_Time-Stratified_", metric_type, "_J_events_GLOBAL.png"))
      ggsave(filename = save_path, plot = plot_grid_j_global, width = 20, height = 5 * ceiling(length(plots_j_global)/3), dpi = 300, limitsize = FALSE)
      cat(paste0("  J-type ", metric_type, " (Global) saved: ", basename(save_path), "\n"))

      # Independent Scale
      plots_j_local <- mapply(plot_heatmap, j_list, titles_j,
                              MoreArgs = list(scale_mode = "local", legend_label = legend_label),
                              SIMPLIFY = FALSE)
      plot_grid_j_local <- gridExtra::grid.arrange(
        grobs = plots_j_local, ncol = 3,
        top = grid::textGrob(paste0("Time-Stratified J-type ", metric_type, " (Independent Scales)"),
                             gp = grid::gpar(fontsize = 16, fontface = "bold"))
      )
      save_path <- file.path(savedir, paste0(model_name, "_Time-Stratified_", metric_type, "_J_events_INDEPENDENT.png"))
      ggsave(filename = save_path, plot = plot_grid_j_local, width = 20, height = 5 * ceiling(length(plots_j_local)/3), dpi = 300, limitsize = FALSE)
      cat(paste0("  J-type ", metric_type, " (Independent) saved: ", basename(save_path), "\n"))

    } else {
      cat(paste0("  All J-type ", metric_type, " values are zero. Skipping J-type plots.\n"))
    }

    # --- TOTAL (d+j) ---
    titles_total <- names(total_list)

    # Global Scale
    plots_total_global <- mapply(plot_heatmap, total_list, titles_total,
                                 MoreArgs = list(scale_mode = "global", color_limits = global_color_limits, legend_label = legend_label),
                                 SIMPLIFY = FALSE)
    plot_grid_total_global <- gridExtra::grid.arrange(
      grobs = plots_total_global, ncol = 3,
      top = grid::textGrob(paste0("Time-Stratified Total (d+j) ", metric_type, " (Global Scale)"),
                           gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    save_path <- file.path(savedir, paste0(model_name, "_Time-Stratified_", metric_type, "_TOTAL_GLOBAL.png"))
    ggsave(filename = save_path, plot = plot_grid_total_global, width = 20, height = 5 * ceiling(length(plots_total_global)/3), dpi = 300, limitsize = FALSE)
    cat(paste0("  Total ", metric_type, " (Global) saved: ", basename(save_path), "\n"))

    # Independent Scale
    plots_total_local <- mapply(plot_heatmap, total_list, titles_total,
                                MoreArgs = list(scale_mode = "local", legend_label = legend_label),
                                SIMPLIFY = FALSE)
    plot_grid_total_local <- gridExtra::grid.arrange(
      grobs = plots_total_local, ncol = 3,
      top = grid::textGrob(paste0("Time-Stratified Total (d+j) ", metric_type, " (Independent Scales)"),
                           gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    save_path <- file.path(savedir, paste0(model_name, "_Time-Stratified_", metric_type, "_TOTAL_INDEPENDENT.png"))
    ggsave(filename = save_path, plot = plot_grid_total_local, width = 20, height = 5 * ceiling(length(plots_total_local)/3), dpi = 300, limitsize = FALSE)
    cat(paste0("  Total ", metric_type, " (Independent) saved: ", basename(save_path), "\n"))
  }

  # --- Generate FLUX heatmaps ---
  generate_heatmap_set(
    d_list = avg_d_timeslice_flux_list,
    j_list = avg_j_timeslice_flux_list,
    total_list = avg_total_timeslice_flux_list,
    metric_type = "FLUX",
    legend_label = "Flux\n(events/Myr)",
    savedir = savedir,
    model_name = model_name
  )

  # --- Generate RATE heatmaps ---
  generate_heatmap_set(
    d_list = avg_d_timeslice_rate_list,
    j_list = avg_j_timeslice_rate_list,
    total_list = avg_total_timeslice_rate_list,
    metric_type = "RATE",
    legend_label = "Rate\n(events/\nbranch-Myr)",
    savedir = savedir,
    model_name = model_name
  )

  cat("\n--- All FLUX and RATE heatmaps generated successfully ---\n")

  # --- 7.5.1 Generate Chord Diagrams (Flux and Rate) ---
  cat("\n--- Generating chord diagrams ---\n")

  # Calculate global matrices for d, j, and total
  # For FLUX: use counts / total duration
  global_d_flux_matrix <- Reduce("+", avg_d_timeslice_list) / sum(durations_Ma)
  global_j_flux_matrix <- Reduce("+", avg_j_timeslice_list) / sum(durations_Ma)
  global_total_flux_matrix <- Reduce("+", avg_total_timeslice_list) / sum(durations_Ma)

  # For RATE: use counts / total branch length
  global_d_rate_matrix <- Reduce("+", avg_d_timeslice_list) / total_branch_length_Myrs
  global_j_rate_matrix <- Reduce("+", avg_j_timeslice_list) / total_branch_length_Myrs
  global_total_rate_matrix <- Reduce("+", avg_total_timeslice_list) / total_branch_length_Myrs

  # --- FLUX Chord Diagrams (6 diagrams: d, j, total × global/timeslice) ---

  # d-event FLUX
  generate_chord_diagrams_set(
    global_matrix = global_d_flux_matrix,
    timeslice_list = avg_d_timeslice_flux_list,
    metric_type = "FLUX_d-events",
    savedir = savedir,
    model_name = model_name
  )

  # j-event FLUX
  generate_chord_diagrams_set(
    global_matrix = global_j_flux_matrix,
    timeslice_list = avg_j_timeslice_flux_list,
    metric_type = "FLUX_j-events",
    savedir = savedir,
    model_name = model_name
  )

  # Total (d+j) FLUX
  generate_chord_diagrams_set(
    global_matrix = global_total_flux_matrix,
    timeslice_list = avg_total_timeslice_flux_list,
    metric_type = "FLUX_total",
    savedir = savedir,
    model_name = model_name
  )

  # --- RATE Chord Diagrams (6 diagrams: d, j, total × global/timeslice) ---

  # d-event RATE
  generate_chord_diagrams_set(
    global_matrix = global_d_rate_matrix,
    timeslice_list = avg_d_timeslice_rate_list,
    metric_type = "RATE_d-events",
    savedir = savedir,
    model_name = model_name
  )

  # j-event RATE
  generate_chord_diagrams_set(
    global_matrix = global_j_rate_matrix,
    timeslice_list = avg_j_timeslice_rate_list,
    metric_type = "RATE_j-events",
    savedir = savedir,
    model_name = model_name
  )

  # Total (d+j) RATE
  generate_chord_diagrams_set(
    global_matrix = global_total_rate_matrix,
    timeslice_list = avg_total_timeslice_rate_list,
    metric_type = "RATE_total",
    savedir = savedir,
    model_name = model_name
  )

  cat("\n--- All chord diagrams generated successfully (12 diagrams total) ---\n")

  # --- 7.6 Generate all temporal trend plots ---
  cat("\n--- Generating temporal trend plots ---\n")
  
  tree_age <- max(branching.times(tr))
  binning_size <- 10 # Set a uniform bin size for time-based event counting.
  
  # --- Main Figure ---
  save_data_path_main<-file.path(savedir, paste0(model_name, "_main_Immigration_vs_Emigration.rds"))
  main_plot <- plot_immigration_emigration_trends(
    clado_events_tables = clado_events_tables_src,
    ana_events_tables = ana_events_tables_src,
    areanames = areanames, tree_age = tree_age, bin_size = binning_size,
    save_data_path = save_data_path_main
  )
  ggsave(file.path(savedir, paste0(model_name, "_main_Immigration_vs_Emigration.png")), 
         plot = main_plot, width = 12, height = 8, dpi = 300)
  cat("Main trend plot (Immigration vs. Emigration) saved.\n")
  
  # --- Supplementary Figure 1 ---
  save_data_path_supp1<-file.path(savedir, paste0(model_name, "_supp_Immigration_Components.rds"))
  supp_plot1 <- plot_immigration_components(
    clado_events_tables = clado_events_tables_src,
    ana_events_tables = ana_events_tables_src,
    areanames = areanames, tree_age = tree_age, bin_size = binning_size,
    save_data_path = save_data_path_supp1
  )
  ggsave(file.path(savedir, paste0(model_name, "_supp_Immigration_Components.png")), 
         plot = supp_plot1, width = 12, height = 8, dpi = 300)
  cat("Supplementary plot (Immigration Components d vs j) saved.\n")
  
  # --- Supplementary Figure 2 ---
  save_data_path_supp2<-file.path(savedir, paste0(model_name, "_supp_Emigration_Components.rds"))
  supp_plot2 <- plot_emigration_components(
    clado_events_tables = clado_events_tables_src,
    ana_events_tables = ana_events_tables_src,
    areanames = areanames, tree_age = tree_age, bin_size = binning_size,
    save_data_path = save_data_path_supp2
  )
  ggsave(file.path(savedir, paste0(model_name, "_supp_Emigration_Components.png")), 
         plot = supp_plot2, width = 12, height = 8, dpi = 300)
  cat("Supplementary plot (Emigration Components d vs j) saved.\n")

  save_data_path_details<-file.path(savedir, paste0(model_name, "_transition_details.rds"))
  family_name <- sub("_.*", "", basename(savedir))
  extract_and_save_inter_realm_events(clado_events_tables = clado_events_tables_src, 
                                      ana_events_tables  = ana_events_tables_src,
                                      family_name=family_name,
                                      tree_age =tree_age, 
                                      bin_size = binning_size,
                                      save_path = save_data_path_details)

  cat("\nBioGeoBEARS visualization finished. All 6 plot files generated.\n")
  
  
}

check_ML_vs_BSM <- function(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE){
  # We assessed the robustness of the Maximum Likelihood (ML) ancestral state reconstructions by comparing them with the results from a Bayesian Stochastic Mapping (BSM) analysis. 
  # Using the check_ML_vs_BSM function in the R package BioGeoBEARS, we generated a scatter plot of the ML marginal probabilities versus the mean state probabilities averaged across 100 stochastic maps for every state at each node in the phylogeny. 
  # The concordance between these two estimation methods was quantified using the R-squared (R²) value of a linear regression. A high R² value indicates that the ML estimates are a reliable representation of the mean posterior probabilities from the BSMs.
  
  # Get tree if needed
  if (is.null(tr))
  {
    #tr = read.tree(res$inputs$trfn)
    tr = check_trfn(trfn=res$inputs$trfn)
  } # END if (is.null(tr))
  
  # Determine if a stratified analysis if needed
  if (is.null(res$inputs$stratified))
  {
    if (is.numeric(res$inputs$timeperiods) == TRUE)
    {
      res$inputs$stratified = TRUE
    } else {
      res$inputs$stratified = FALSE
    }# END if (is.numeric(BioGeoBEARS_run_object$timeperiods) == TRUE)
  } # END if (is.null(res$inputs$stratified))
  stratified = res$inputs$stratified
  
  
  pdffn = paste0(model_name, "_ML_vs_BSM.pdf")
  pdf(file=pdffn, height=6, width=6)
  
  x_all = NULL
  y_all = NULL
  numtips = length(tr$tip.label)
  intnodenums = (numtips+1):(numtips+tr$Nnode)
  
  MLstateprobs = res$ML_marginal_prob_each_state_at_branch_top_AT_node
  numstates = ncol(MLstateprobs)
  BSMstates_summary = calc_BSM_mean_node_states(clado_events_tables, tr, numstates)
  state_counts_1based = BSMstates_summary$state_counts_1based
  meanBSMprobs = BSMstates_summary$meanBSMprobs
  sumBSMcounts = BSMstates_summary$sumBSMcounts
  
  node_history_samples_allNodes = NULL
  cat("\nCalculating BSM means for node #:", sep="")
  
  for (nodenum in intnodenums)
  {
    cat(nodenum, " ", sep="")
    
    node_history_samples = NULL
    for (i in 1:length(clado_events_tables))
    {
      clado_events_table = clado_events_tables[[i]]
      if (stratified == TRUE)
      {
        TF1 = clado_events_table$SUBnode.type != "tip"
        TF2 = clado_events_table$node.type != "tip"
        TF = (TF1+TF2)==2
      } else {
        TF = clado_events_table$node.type != "tip"
      }
      clado_events_table = clado_events_table[TF,]
      TF = clado_events_table$node == nodenum
      node_history_sample = clado_events_table[TF,]
      node_history_samples = rbind(node_history_samples, node_history_sample)
    } # END for (i in 1:length(clado_events_tables))
    head(node_history_samples)
    class(node_history_samples)
    dim(node_history_samples)
    
    # Save
    node_history_samples_allNodes = rbind(node_history_samples_allNodes, node_history_samples)
    
    node_history_samples$clado_event_txt
    node_history_samples$sampled_states_AT_nodes
    
    
    cbind(node_history_samples$sampled_states_AT_nodes, node_history_samples$sampled_states_AT_brbots)
    table(node_history_samples$sampled_states_AT_nodes)
    table(node_history_samples$sampled_states_AT_brbots)
    
    x=round(res$ML_marginal_prob_each_state_at_branch_top_AT_node[nodenum,],3)
    
    num_root_state = table(node_history_samples$sampled_states_AT_nodes)
    fract_root_state = num_root_state / sum(num_root_state)
    obs_BSM_probs_at_root = rep(0, length(x))
    obs_BSM_probs_at_root[as.numeric(names(fract_root_state))] = fract_root_state
    
    x=round(res$ML_marginal_prob_each_state_at_branch_top_AT_node[nodenum,],3)
    y=obs_BSM_probs_at_root
    x_all = c(x_all, x)
    y_all = c(y_all, y)
    
    if (plot_each_node == TRUE)
    {
      plot(x,y, xlim=c(0,1), ylim=c(0,1), xlabel="ML marginal probs", ylabel="BSM probs")
      segments(0,0,1,1)
      title(nodenum)
    } # END if (plot_each_node == TRUE)
  } # END for (nodenum in 20:37)
  cat("...done.\n\n")
  
  # Plot allnodes
  if (linreg_plot == FALSE)
  {
    plot(x_all, y_all, xlim=c(0,1), ylim=c(0,1), xlabel="ML marginal probs", ylabel="BSM probs")
    segments(0,0,1,1)
    title("all nodes")
  } # END if (linreg_plot == FALSE)
  
  if (linreg_plot == TRUE)
  {
    linear_regression_plot(x=x_all, y=y_all, tmppch=1, xlabel="State probabilities under ML model", ylabel="State probabilities as mean of BSMs", xlim=c(0,1), ylim=c(0,1))
    # Multiple R-squared:  0.9655,	Adjusted R-squared:  0.9655 
    title(paste0(model_name, ":\nML state probs vs. mean of BSMs"))
  } # END if (linreg_plot == TRUE)
  
  if (MultinomialCI == TRUE)
  {
    # Try doing confidence intervals
    #require(MultinomialCI)
    sumBSMcounts2 = sumBSMcounts[-(1:numtips),]
    
    # Confident intervals for a particular node
    # CI95 = multinomialCI(x=sumBSMcounts2[1,], alpha=0.05, verbose=TRUE)
    # CI95
    
    # Add multinomial CI95s to plot
    cat("\nAdding CI95 segments for: ", sep="")
    for (i in 1:tr$Nnode)
    {
      cat(i+numtips, " ", sep="")
      tmpcounts = sumBSMcounts2[i,]
      CI95 = MultinomialCI::multinomialCI(x=tmpcounts, alpha=0.05, verbose=FALSE)
      #xvals = tmpcounts / sum(tmpcounts)
      xvals = MLstateprobs[-(1:numtips),][i,]
      segments(x0=xvals, x1=xvals, y0=CI95[,1], y1=CI95[,2])
    } # for (1 in 1:numnodes)
  } # END if (MultinomialCI == TRUE)
  cat("...done.\n\n")
  
  # Close PDF and open for viewer
  dev.off()
  #cmdstr = paste("open ", pdffn)
  #system(cmdstr)
}


# Function 1 (Main Figure): Plot temporal dynamics of Immigration vs. Emigration
plot_immigration_emigration_trends <- function(clado_events_tables, 
                                               ana_events_tables, 
                                               areanames,
                                               tree_age,
                                               bin_size = 10,
                                               save_data_path) {
  
  time_breaks <- seq(0, ceiling(tree_age / bin_size) * bin_size, by = bin_size)
  time_labels <- paste(time_breaks[-length(time_breaks)], time_breaks[-1], sep="-")
  num_sims <- length(clado_events_tables)
  
  results_per_sim <- list()
  for (i in 1:num_sims) {
    clado_table <- clado_events_tables[[i]]
    ana_table <- ana_events_tables[[i]]
    
    events_in_list <- list()
    if(is.data.frame(ana_table) && nrow(ana_table) > 0) {
      d_in_subset <- ana_table %>% filter(!is.na(event_type) & event_type == "d") %>% select(time=abs_event_time, area=dispersal_to)
      if(nrow(d_in_subset) > 0) events_in_list[['d']] <- d_in_subset
    }
    if(is.data.frame(clado_table) && nrow(clado_table) > 0) {
      j_in_subset <- clado_table %>% filter(!is.na(clado_event_type) & clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_to)
      if(nrow(j_in_subset) > 0) events_in_list[['j']] <- j_in_subset
    }
    
    if (length(events_in_list) > 0) {
      imm_counts <- bind_rows(events_in_list) %>%
        mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
        filter(!is.na(time_bin)) %>%
        group_by(time_bin, area) %>% 
        summarise(count = dplyr::n(), .groups = 'drop') %>% 
        mutate(event_type = "Immigration")
    } else {
      imm_counts <- data.frame(time_bin=factor(), area=character(), count=integer(), event_type=character())
    }
    
    events_out_list <- list()
    if(is.data.frame(ana_table) && nrow(ana_table) > 0) {
      d_out_subset <- ana_table %>% filter(!is.na(event_type) & event_type == "d") %>% select(time=abs_event_time, area=ana_dispersal_from)
      if(nrow(d_out_subset) > 0) events_out_list[['d']] <- d_out_subset
    }
    if(is.data.frame(clado_table) && nrow(clado_table) > 0) {
      j_out_subset <- clado_table %>% filter(!is.na(clado_event_type) & clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_from)
      if(nrow(j_out_subset) > 0) events_out_list[['j']] <- j_out_subset
    }
    
    if(length(events_out_list) > 0) {
      emi_counts <- bind_rows(events_out_list) %>%
        mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
        filter(!is.na(time_bin) & !is.na(area)) %>%
        group_by(time_bin, area) %>% 
        summarise(count = dplyr::n(), .groups = 'drop') %>% 
        mutate(event_type = "Emigration")
    } else {
      emi_counts <- data.frame(time_bin=factor(), area=character(), count=integer(), event_type=character())
    }
    
    sim_results <- bind_rows(imm_counts, emi_counts)
    if(nrow(sim_results) > 0) {
      sim_results$sim_id <- i
      results_per_sim[[i]] <- sim_results
    }
  }
  
  all_sim_results <- bind_rows(results_per_sim)
  
  full_grid <- expand.grid(time_bin = time_labels, area = areanames, 
                           event_type = c("Immigration", "Emigration"),
                           sim_id = 1:num_sims, stringsAsFactors = FALSE)
  
  summary_data <- full_grid %>%
    left_join(all_sim_results, by = c("time_bin", "area", "event_type", "sim_id")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  final_summary <- summary_data %>%
    group_by(time_bin, area, event_type) %>%
    summarise(mean_events = mean(count),
              lower_q = quantile(count, 0.05),
              upper_q = quantile(count, 0.95), .groups = 'drop') %>%
    mutate(time_mid = as.numeric(sub("-.*", "", time_bin)) + (bin_size / 2))
  
  saveRDS(final_summary,file = save_data_path)

  p <- ggplot(final_summary, aes(x = time_mid, y = mean_events, color = event_type, fill = event_type)) +
    geom_ribbon(aes(ymin = lower_q, ymax = upper_q), alpha = 0.2, linetype = 0) +
    geom_line(size = 1) +
    facet_wrap(~ area, scales = "fixed") +
    scale_color_manual(name = "Process", values = c("Immigration" = "darkblue", "Emigration" = "darkred")) +
    scale_fill_manual(name = "Process", values = c("Immigration" = "lightblue", "Emigration" = "lightcoral")) +
    labs(title = "Temporal Dynamics of Immigration vs. Emigration by Region",
         subtitle = paste("Averaged over", num_sims, "BSM scenarios in", bin_size, "Ma bins"),
         x = "Time (Million years ago)", y = "Average Number of Events per Time Bin") +
    scale_x_reverse() + theme_bw() + theme(legend.position = "bottom")
  
  return(p)
}

# Function 2 (Supplementary Figure 1): Plot components of Immigration events (d vs. j)
plot_immigration_components <- function(clado_events_tables, 
                                        ana_events_tables, 
                                        areanames,
                                        tree_age,
                                        bin_size = 10,
                                        save_data_path) {
  
  time_breaks <- seq(0, ceiling(tree_age / bin_size) * bin_size, by = bin_size)
  time_labels <- paste(time_breaks[-length(time_breaks)], time_breaks[-1], sep="-")
  num_sims <- length(clado_events_tables)
  
  results_per_sim <- list()
  for (i in 1:num_sims) {
    clado_table <- clado_events_tables[[i]]
    ana_table <- ana_events_tables[[i]]
    
    events_in_list <- list()
    if(is.data.frame(ana_table) && nrow(ana_table) > 0) {
      d_in_subset <- ana_table %>% filter(!is.na(event_type) & event_type == "d") %>% select(time=abs_event_time, area=dispersal_to) %>% mutate(event_type="d (Range Expansion)")
      if(nrow(d_in_subset) > 0) events_in_list[['d']] <- d_in_subset
    }
    if(is.data.frame(clado_table) && nrow(clado_table) > 0) {
      j_in_subset <- clado_table %>% filter(!is.na(clado_event_type) & clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_to) %>% mutate(event_type="j (Founder Event)")
      if(nrow(j_in_subset) > 0) events_in_list[['j']] <- j_in_subset
    }
    
    if (length(events_in_list) > 0) {
      sim_results <- bind_rows(events_in_list) %>%
        mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
        filter(!is.na(time_bin)) %>%
        group_by(time_bin, area, event_type) %>% 
        summarise(count = dplyr::n(), .groups = 'drop') 
      
      if(nrow(sim_results) > 0) {
        sim_results$sim_id <- i
        results_per_sim[[i]] <- sim_results
      }
    }
  }
  
  all_sim_results <- bind_rows(results_per_sim)
  
  full_grid <- expand.grid(time_bin = time_labels, area = areanames, 
                           event_type = c("d (Range Expansion)", "j (Founder Event)"),
                           sim_id = 1:num_sims, stringsAsFactors = FALSE)
  
  summary_data <- full_grid %>%
    left_join(all_sim_results, by = c("time_bin", "area", "event_type", "sim_id")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  final_summary <- summary_data %>%
    group_by(time_bin, area, event_type) %>%
    summarise(mean_events = mean(count),
              lower_q = quantile(count, 0.05),
              upper_q = quantile(count, 0.95), .groups = 'drop') %>%
    mutate(time_mid = as.numeric(sub("-.*", "", time_bin)) + (bin_size / 2))
  
  saveRDS(final_summary,file = save_data_path)
  
  p <- ggplot(final_summary, aes(x = time_mid, y = mean_events, color = event_type, fill = event_type)) +
    geom_ribbon(aes(ymin = lower_q, ymax = upper_q), alpha = 0.2, linetype = 0) +
    geom_line(size = 1) +
    facet_wrap(~ area, scales = "free_y") +
    labs(title = "Components of Immigration Events by Region",
         subtitle = paste("Averaged over", num_sims, "BSM scenarios in", bin_size, "Ma bins"),
         x = "Time (Million years ago)", y = "Average Number of Events per Time Bin",
         color = "Immigration Type", fill = "Immigration Type") +
    scale_x_reverse() + theme_bw() + theme(legend.position = "bottom")
  
  return(p)
}

# Function 3 (Supplementary Figure 2): Plot components of Emigration events (d vs. j)
plot_emigration_components <- function(clado_events_tables, 
                                       ana_events_tables, 
                                       areanames,
                                       tree_age,
                                       bin_size = 10,
                                       save_data_path) {
  
  time_breaks <- seq(0, ceiling(tree_age / bin_size) * bin_size, by = bin_size)
  time_labels <- paste(time_breaks[-length(time_breaks)], time_breaks[-1], sep="-")
  num_sims <- length(clado_events_tables)
  
  results_per_sim <- list()
  for (i in 1:num_sims) {
    clado_table <- clado_events_tables[[i]]
    ana_table <- ana_events_tables[[i]]
    
    events_out_list <- list()
    if(is.data.frame(ana_table) && nrow(ana_table) > 0) {
      d_out_subset <- ana_table %>% filter(!is.na(event_type) & event_type == "d") %>% select(time=abs_event_time, area=ana_dispersal_from) %>% mutate(event_type="d (Range Expansion)")
      if(nrow(d_out_subset) > 0) events_out_list[['d']] <- d_out_subset
    }

    if(is.data.frame(clado_table) && nrow(clado_table) > 0) {
      j_out_subset <- clado_table %>% filter(!is.na(clado_event_type) & clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_from) %>% mutate(event_type="j (Founder Event)")
      if(nrow(j_out_subset) > 0) events_out_list[['j']] <- j_out_subset
    }
    
    if (length(events_out_list) > 0) {
      sim_results <- bind_rows(events_out_list) %>%
        mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
        filter(!is.na(time_bin) & !is.na(area)) %>%
        group_by(time_bin, area, event_type) %>% 
        summarise(count = dplyr::n(), .groups = 'drop') 
      
      if(nrow(sim_results) > 0) {
        sim_results$sim_id <- i
        results_per_sim[[i]] <- sim_results
      }
    }
  }
  
  all_sim_results <- bind_rows(results_per_sim)
  
  full_grid <- expand.grid(time_bin = time_labels, area = areanames, 
                           event_type = c("d (Range Expansion)", "j (Founder Event)"),
                           sim_id = 1:num_sims, stringsAsFactors = FALSE)
  
  summary_data <- full_grid %>%
    left_join(all_sim_results, by = c("time_bin", "area", "event_type", "sim_id")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  final_summary <- summary_data %>%
    group_by(time_bin, area, event_type) %>%
    summarise(mean_events = mean(count),
              lower_q = quantile(count, 0.05),
              upper_q = quantile(count, 0.95), .groups = 'drop') %>%
    mutate(time_mid = as.numeric(sub("-.*", "", time_bin)) + (bin_size / 2))
  
  saveRDS(final_summary,file = save_data_path)
  
  p <- ggplot(final_summary, aes(x = time_mid, y = mean_events, color = event_type, fill = event_type)) +
    geom_ribbon(aes(ymin = lower_q, ymax = upper_q), alpha = 0.2, linetype = 0) +
    geom_line(size = 1) +
    facet_wrap(~ area, scales = "free_y") +
    labs(title = "Components of Emigration Events by Region",
         subtitle = paste("Averaged over", num_sims, "BSM scenarios in", bin_size, "Ma bins"),
         x = "Time (Million years ago)", y = "Average Number of Events per Time Bin",
         color = "Emigration Type", fill = "Emigration Type") +
    scale_x_reverse() + theme_bw() + theme(legend.position = "bottom")
  
  return(p)
}



# Load necessary R packages
library(dplyr)
library(tidyr)

#' @title Extract and save detailed inter-realm dispersal events for a single family
#'
#' @description This function takes the BSM output for a single family, extracts 
#' all anagenetic ('d') and cladogenetic ('j') dispersal events, formats them 
#' into a tidy dataframe with all necessary details, and saves it to a file.
#'
#' @param clado_events_tables The RES_clado_events_tables from a BSM output.
#' @param ana_events_tables The RES_ana_events_tables from a BSM output.
#' @param family_name A character string identifying the family (e.g., "Coccinellidae"). 
#'                    This is ESSENTIAL for downstream comparative analysis.
#' @param tree_age The root age of the tree, used to define time bins.
#' @param bin_size The duration of each time slice in millions of years.
#' @param save_path The file path where the output dataframe (.rds) will be saved.

extract_and_save_inter_realm_events <- function(clado_events_tables, 
                                                ana_events_tables,
                                                family_name,
                                                tree_age, 
                                                bin_size,
                                                save_path) {
  
  # --- 1. Define time bins ---
  time_breaks <- seq(0, ceiling(tree_age / bin_size) * bin_size, by = bin_size)
  time_labels <- paste(time_breaks[-length(time_breaks)], time_breaks[-1], sep="-")
  num_sims <- length(clado_events_tables)
  
  # --- 2. Iterate through each BSM simulation to extract the event list ---
  events_from_all_sims <- lapply(1:num_sims, function(i) {
    clado_table <- clado_events_tables[[i]]
    ana_table <- ana_events_tables[[i]]
    
    # A list to store 'd' and 'j' events for the current simulation
    events_this_sim <- list()
    
    # Extract 'd' events (Anagenetic)
    if (is.data.frame(ana_table) && nrow(ana_table) > 0) {
      events_this_sim[['d']] <- ana_table %>%
        filter(!is.na(event_type) & event_type == "d") %>%
        select(event_time = abs_event_time,
               source_realm = ana_dispersal_from,
               target_realm = dispersal_to) %>%
        mutate(event_type = "d")
    }

    # Extract 'j' events (Cladogenetic)
    if (is.data.frame(clado_table) && nrow(clado_table) > 0) {
      events_this_sim[['j']] <- clado_table %>%
        filter(!is.na(clado_event_type) & clado_event_type == "founder (j)") %>%
        select(event_time = time_bp,
               source_realm = clado_dispersal_from,
               target_realm = clado_dispersal_to) %>%
        mutate(event_type = "j")
    }
    
    # Combine all events for the current simulation
    if (length(events_this_sim) > 0) {
      bind_rows(events_this_sim) %>% mutate(map_num = i)
    } else {
      NULL
    }
  })
  
  # --- 3. Aggregate, clean, and format the data ---
  final_df <- bind_rows(events_from_all_sims) %>%
    # Clean up invalid events (e.g., missing source or target)
    filter(!is.na(source_realm), !is.na(target_realm), 
           source_realm != "", target_realm != "",
           source_realm != "NA", target_realm != "NA") %>%
    # Add family identifier and assign time bins
    mutate(
      family = family_name,
      time_bin = cut(event_time, breaks = time_breaks, labels = time_labels, right = FALSE)
    ) %>%
    filter(!is.na(time_bin)) %>%
    # Reorder columns to the final desired format
    select(family, map_num, event_type, event_time, time_bin, source_realm, target_realm)
  
  # --- 4. Save the data ---
  if (nrow(final_df) > 0) {
    cat("Saving", nrow(final_df), "inter-realm events for", family_name, "to:", save_path, "\n")
    saveRDS(final_df, file = save_path)
  } else {
    warning("No valid inter-realm dispersal events found for family: ", family_name)
  }
  
  # --- 5. (Optional) Return the dataframe invisibly ---
  invisible(final_df)
}



# #test
# # Path to the FASTA file containing all barcode sequences. This includes barcodes 
# # for the target family, distant groups, and sister families. For taxa with full 
# # mitogenomes, their corresponding barcode regions are also included in this alignment.
# allsequences_path <- "C:\\Users\\16575\\Documents\\Rdocuments\\aligned_fin_k13k14_barcod_renamed_20250815.fasta"
# # Path to the FASTA file containing the complete mitogenome sequences (supermatrix). 
# # This file includes all mitogenomes used in the study, encompassing the distant and sister groups.
# mitogenomes_path <- "C:\\Users\\16575\\Documents\\Rdocuments\\5_Supermatrix_13PCG_NT_cleaned_rename_20250822.fasta"
# 
# # Reference to the list of sister groups for each family (defined previously).
# sister_group_family  
# 
# # Reference to the list of distant groups. These will serve as outgroups for rooting, 
# # while sister groups will be part of the ingroup for dating calibrations.
# distantly_group
# 
# #Length of mitogenomes
# mito_length<-11000
# 
# #Time-Stratified File
# timeperiods_filepath = "C:\\Users\\16575\\Documents\\Rdocuments\\my_timeperiods_20mya.txt"
# dispersal_multipliers_filepath = "C:\\Users\\16575\\Documents\\Rdocuments\\my_dispersal_multipliers_20mya.txt"
# 
# threads<-14


tree_biogeography_pipeline<-function(allsequences_path,
                                     mitogenomes_path,
                                     mito_length,
                                     sister_group_family,
                                     distantly_group,
                                     threads=10,
                                     construction_model="GTR+G",
                                     timeperiods_filepath = timeperiods_filepath,
                                     dispersal_multipliers_filepath = dispersal_multipliers_filepath){
  output_dir <- dirname(allsequences_path)
  base_name <- sub("\\.fasta$", "", basename(allsequences_path))

  # Define output file names.
  backbone_tree_path <- file.path(output_dir, paste0(base_name, "_backbone.nwk"))
  final_tree_path <- file.path(output_dir, paste0(base_name, "_reconstruction.nwk"))
  backbone_rooted_file <- sub("\\.nwk$", "_rooted.nwk", backbone_tree_path)
  final_rooted_file <- sub("\\.nwk$", "_rooted.nwk", final_tree_path)
  backbone_dated_file <- sub("\\.nwk$", "_dated_treePL.nwk", backbone_rooted_file)
  dated_tree_path <- file.path(dirname(final_rooted_file), paste0(tools::file_path_sans_ext(basename(final_rooted_file)), "_dated_chronos.nwk"))
  
  construction_tree(nwk_file = NULL,
                    sister_group_family = sister_group_family,
                    distantly_group = distantly_group,
                    allsequences_path =allsequences_path,
                    mitogenomes_path = mitogenomes_path,
                    threads = threads,
                    construction_model=construction_model
  )
  

  
  rooting_tree(undated_tree = backbone_tree_path,
               distantly_group = distantly_group,
               sister_group_family = sister_group_family)
  rooting_tree(undated_tree = final_tree_path,
               distantly_group = distantly_group,
               sister_group_family = sister_group_family)
  

  
  
  dating_tree(rooted_tree_path = backbone_rooted_file,
              sister_group_family = sister_group_family,
              num_sites = mito_length
  )
  fix_backbone_dating_to_reconstruction(backbone_dated_path = backbone_dated_file,
                                        reconstruction_path = final_rooted_file)
  
  
  
  
  
  pastml_process_tree(tree_path = dated_tree_path,
                      geography_csv = NULL  # Will use legacy mode (extract from tip labels)
  )
  result_info <- run_biogeobears_pipeline(
    tree_filepath = dated_tree_path,
    geography_csv = NULL,  # Will use legacy mode (extract from tip labels)
    timeperiods_filepath = timeperiods_filepath,
    dispersal_multipliers_filepath = dispersal_multipliers_filepath,
    num_cores = threads
  )
  
  
  
  pipeline_results <- list(
    paths = list(
      backbone_unrooted = backbone_tree_path,
      reconstruction_unrooted = final_tree_path,
      backbone_rooted = backbone_rooted_file,
      reconstruction_rooted = final_rooted_file,
      backbone_dated = backbone_dated_file,
      final_dated_tree = dated_tree_path
    ),
    biogeobears_analysis = result_info
  )
  
  cat("\n--- Pipeline Finished: All results are collected. ---\n")
  
  return(pipeline_results)
}


# sister_group_family <- list()
# sister_group_family$Staphylinidae<-c("Silphidae_135","Leiodidae_185","Agyrtidae_185","Hydraenidae_185","Ptiliidae_185")
# distantly_group<-list()
# distantly_group$Staphylinidae<- c("Scarabaeidae", "Tenebrionidae")   
# allsequences_path<-"C:\\Users\\AD\\Documents\\visit_nhm\\test_20250913\\test_staphylinidae_barcodes.fasta"
# mitogenomes_path<-"C:\\Users\\AD\\Documents\\visit_nhm\\test_20250913\\test_staphylinidae_mitogenomes.fasta"
# timeperiods_filepath = "C:\\Users\\AD\\Documents\\visit_nhm\\test_20250913\\my_timeperiods_20mya.txt"
# dispersal_multipliers_filepath = "C:\\Users\\AD\\Documents\\visit_nhm\\test_20250913\\my_dispersal_multipliers_20mya.txt"

# results<-tree_biogeography_pipeline(allsequences_path=allsequences_path,
#                            mitogenomes_path=mitogenomes_path,
#                            mito_length=11000,
#                            sister_group_family=sister_group_family,
#                            distantly_group=distantly_group,
#                            threads = 12,
#                            timeperiods_filepath = timeperiods_filepath,
#                            dispersal_multipliers_filepath = dispersal_multipliers_filepath,
#                            construction_model = "GTR+G")




