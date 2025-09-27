#20250927
# The pipeline consists of the following steps:
# construction_tree: Constructs the phylogenetic tree.
# rooting_tree: Roots the tree generated in the previous step.
# dating_tree: Dates the backbone tree. While it's possible to date the reconstruction tree directly,
# for the reconstruction tree, we opt to use the dating results from the backbone via the fix_backbone_dating_to_reconstruction function.
# biogeography_analysis_pastml: Performs biogeographical analysis using PastML.
# biogeography_analysis_biogeobears: Performs biogeographical analysis using BioGeoBEARS.

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
  parts[length(parts) - 1]  # The second-to-last element is the taxon name.
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
  
  # Root the tree using the cleaned outgroup. resolve.root = TRUE ensures a bifurcating root, which is standard practice.
  rooted_tree <- root(pruned_tree, outgroup = my_outgroups, resolve.root = TRUE)
  
  # Write the final rooted tree to a new file.
  file_dir  <- dirname(use_nwk_tree)
  file_base <- basename(use_nwk_tree)
  file_base_clean <- sanitize_filename(file_base)
  use_nwk_tree <- file.path(file_dir, file_base_clean)
  rooted_file <-ifelse(grepl("\\.nwk$", use_nwk_tree),
                       sub("\\.nwk$", "_rooted.nwk", use_nwk_tree),
                       paste0(use_nwk_tree, "_rooted.nwk"))
  
  write.tree(rooted_tree, file = rooted_file)
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

fix_backbone_dating_to_reconstruction <- function(backbone_dated_path, reconstruction_path) {
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
  
  # 4. Prepare the calibration data frame for chronos().
  print("Step 2: Preparing calibration data for chronos()...")
  
  # a) Get node age information from the dated backbone tree.
  backbone_ages_map <- setNames(ape::branching.times(backbone_tree),
                                (length(backbone_tree$tip.label) + 1):(length(backbone_tree$tip.label) + backbone_tree$Nnode))
  
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
  
  # 5. Run the chronos() function to perform dating.
  print("Step 3: Running chronos() dating analysis... This may take some time.")
  
  # lambda = 1 is a common value for the smoothing parameter.
  # model = "correlated" specifies the commonly used correlated rates model.
  dated_tree_chronos <- ape::chronos(recon_tree, lambda = 1, model = "correlated", calibration = calib_df)
  
  print("Dating analysis completed!")
  
  # 6. Visualize and save the results.
  print("Step 4: Plotting and saving the final dated tree...")
  
  # Plot the time-calibrated tree with a time axis.
  plot(dated_tree_chronos, cex = 0.5)
  axisPhylo()
  title("Final Tree Dated with ape::chronos()")
  
  # Save the final dated tree to a file.
  output_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(reconstruction_path)), "_dated_chronos.nwk"))
  ape::write.tree(dated_tree_chronos, file = output_path)
  
  print(paste("Final dated tree has been saved to:", output_path))
  
}

# Execute the function to date the reconstruction tree using the backbone's ages.
# fix_backbone_dating_to_reconstruction(backbone_dated_path = backbone_dated_path,
#                                       reconstruction_path = reconstruction_path)





# --- Step 4: Perform Biogeographical Analysis ---
# This step includes analyses with PastML and BioGeoBEARS.
# The required inputs are a dated tree and a table of tip biogeographic states.

cal_transition_rate_tab <- function(named_tree_path, probabilities_filepath) {
  library(pheatmap)
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
          contribution <- P_parent[from_idx] * P_daughter[to_idx]
          total_transition_contributions[from_idx, to_idx] <- total_transition_contributions[from_idx, to_idx] + contribution
        }
      }
    }
  }
  
  # --- 3. Calculate Transition Rates ---
  # Get the total branch length of the tree (unit: million years).
  total_branch_length_myrs <- sum(phy_tree$edge.length)
  
  if (total_branch_length_myrs == 0) {
    stop("Error: Total tree branch length is zero. Cannot calculate rates. Ensure the tree is a time-calibrated phylogram.")
  }
  
  # Transition rate matrix = total probability contribution / total branch length.
  transition_rate_matrix_per_myr <- total_transition_contributions / total_branch_length_myrs
  
  # --- 4. Output Results ---
  print("Total Transition Probability Contribution Matrix:")
  print(round(total_transition_contributions, 4))
  
  print(paste("Total phylogenetic branch length (Myr):", total_branch_length_myrs))
  
  print("Transition Rate Matrix (Avg. Probability Contribution / Myr):")
  print(round(transition_rate_matrix_per_myr, 6))
  
  # Visualization (e.g., heatmap).
  save_path <- sub("\\_probabilities.tab$", "_transition_plots.png", probabilities_filepath)
  # Convert numbers to scientific notation for display on the heatmap.
  numbers_sci <- formatC(transition_rate_matrix_per_myr, format = "e", digits = 2)
  
  pheatmap_obj <- pheatmap::pheatmap(
    transition_rate_matrix_per_myr,
    display_numbers = numbers_sci,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Transition Rates (events/Myr) from PastML (Y-axis: From -> X-axis: To)",
    fontsize_main = 12,
    # Specify filename and dimensions directly here.
    filename = save_path,
    width = 8,  # inches
    height = 7, # inches
    # silent = TRUE prevents the plot from being displayed in the R graphics device.
    silent = TRUE
  )
  print(paste0(save_path, " finished"))
}


pastml_process_tree <- function(dated_tree_path,
                                discard_other_family = TRUE) {
  
  library(tidyverse)
  library(ggtree)
  library(ggplot2)
  library(stringr)
  # dated_tree_path: Path to the dated tree from the previous step.
  # location_path: Path to the tip states file (ID, realm). This is now auto-generated.
  # discard_other_family: If TRUE, removes tips not belonging to the main target family.
  # pastml_path: Absolute path to the PastML executable in WSL.
  
  # Create a dedicated directory for the tree's analysis outputs.
  tree <- dated_tree_path
  win_tree_dir <- paste0(tools::file_path_sans_ext(tree), "_pastml")
  dir.create(win_tree_dir, showWarnings = FALSE)
  
  # Get the path to PastML within WSL.
  pastml_path_raw <- system("wsl bash -ic 'which pastml'", intern = TRUE)
  pastml_path <- tail(pastml_path_raw, 1)
  # Print the located path for verification.
  print(pastml_path)
  
  # Auto-generate the location data file from tip labels.
  use.tree <- read.tree(tree)
  use.tree.tip <- use.tree$tip.label
  extract_realm <- function(tip) {
    parts <- str_split_fixed(tip, "_", 4)
    parts[length(parts)] # The last element is the biogeographic realm.
  }
  tips_realm <- sapply(use.tree.tip, extract_realm)
  location_path <- paste0(win_tree_dir, "/matched_location.csv")
  
  # If realm codes are not single characters, map them.
  if (!all(nchar(tips_realm[tips_realm != ""]) == 1)) {
    realm_map <- c(
      "NA" = "N",  # Nearctic
      "NT" = "T",  # Neotropical
      "PN" = "P",  # Panamanian
      "PA" = "A",  # Palaearctic
      "SA" = "S",  # Saharo-Arabian
      "SJ" = "J",  # Sino-Japanese
      "AT" = "F",  # Afrotropical
      "MA" = "M",  # Madagascan
      "IM" = "I",  # Indomalayan
      "AA" = "U",  # Australasian
      "OC" = "O"   # Oceania
    )
    
    # Dynamically map any realms not in the predefined list.
    unique_realms <- unique(as.character(tips_realm))
    missing_realms <- setdiff(unique_realms, names(realm_map))
    missing_realms <- missing_realms[missing_realms != ""]
    extra_map <- setNames(letters[seq_along(missing_realms)], missing_realms)
    full_map <- c(realm_map, extra_map)
    
    # Apply the mapping.
    tips_realm_mapped <- sapply(tips_realm, function(x) {
      if (x == "") {
        "" # Preserve empty strings.
      } else {
        full_map[[x]]
      }
    })
    tips_realm <- tips_realm_mapped
  }
  
  # Convert the named vector to a data frame and save as CSV.
  tips_df <- data.frame(
    ID = names(tips_realm),
    realm = as.vector(tips_realm),
    stringsAsFactors = FALSE
  )
  write.csv(tips_df, file = location_path, row.names = FALSE)
  
  safe_feature <- "realm"
  
  # Process the tree to ensure only the target family is included.
  use.tree <- read.tree(tree)
  tree_name <- tools::file_path_sans_ext(basename(tree))
  if (discard_other_family) {
    all_tips <- use.tree$tip.label
    families <- sapply(all_tips, extract_family)
    target_family <- names(which.max(table(families)))
    # Identify tips that do not belong to the target family (or are NA).
    tips_to_drop <- names(families)[families != target_family | is.na(families)]
    # Ensure all tips to be dropped actually exist in the tree to prevent errors.
    tips_to_drop <- intersect(tips_to_drop, all_tips)
    # Prune the tips.
    use.tree <- drop.tip(use.tree, tips_to_drop)
  }
  # Save the pruned tree to the working directory.
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
                                     output_tag = "biogeobears",
                                     num_cores = 12,
                                     models_to_run = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"),
                                     discard_other_family = TRUE, # If TRUE, prune tips from other families to focus the analysis.
                                     sydysfinal_list = NULL, # Legacy parameter for handling taxa with multiple geographic records; now handled by data format.
                                     # The two files below are for time-stratified analysis.
                                     timeperiods_filepath = "C:\\Users\\16575\\Documents\\Rdocuments\\my_timeperiods_20mya.txt",
                                     dispersal_multipliers_filepath = "C:\\Users\\16575\\Documents\\Rdocuments\\my_dispersal_multipliers_20mya.txt",
                                     user_max_range_size = 3,
                                     min_brlen = 1e-6,
                                     plot_width_px = 6000,
                                     plot_height_px = 36000,
                                     plot_res_dpi = 400
) {
  library(BioGeoBEARS)
  library(ape)
  library(dplyr)
  library(tidyr)
  library(tools) # For file path manipulation.
  library(phytools)
  library(stringr)
  library(snow)
  
  # dated_tree_path: The dated tree from the previous step.
  # location_path: Tip states file (ID, realm).
  # discard_other_family: If TRUE, non-target families (used for dating/rooting) are pruned.
  
  # --- 1. Validate inputs and set up paths ---
  # Path for caching temporary files.
  backup_dir <- file.path(dirname(tree_filepath), "backup_dir")
  dir.create(backup_dir, showWarnings = FALSE)
  
  # Path for final BioGeoBEARS results.
  save_biogeobears_path <- paste0(tools::file_path_sans_ext(tree_filepath), "_biogeobears_results")
  dir.create(save_biogeobears_path, showWarnings = FALSE)
  
  print("--- 1. Validating inputs and setting up paths ---")
  if (!file.exists(tree_filepath)) stop("Tree file not found: ", tree_filepath)
  
  # Auto-generate the location data file from tip labels.
  use.tree <- read.tree(tree_filepath)
  use.tree.tip <- use.tree$tip.label
  extract_realm <- function(tip) {
    parts <- str_split_fixed(tip, "_", 4)
    parts[length(parts)] # The last element is the biogeographic realm.
  }
  tips_realm <- sapply(use.tree.tip, extract_realm)
  tip_states_filepath <- paste0(save_biogeobears_path, "/matched_location.csv")
  
  # If realm codes are not single characters, map them to single characters.
  if (!all(nchar(tips_realm[tips_realm != ""]) == 1)) {
    realm_map <- c(
      "NA" = "N",  # Nearctic
      "NT" = "T",  # Neotropical
      "PN" = "P",  # Panamanian
      "PA" = "A",  # Palaearctic
      "SA" = "S",  # Saharo-Arabian
      "SJ" = "J",  # Sino-Japanese
      "AT" = "F",  # Afrotropical
      "MA" = "M",  # Madagascan
      "IM" = "I",  # Indomalayan
      "AA" = "U",  # Australasian
      "OC" = "O"   # Oceania
    )
    
    unique_realms <- unique(as.character(tips_realm))
    # Find any realms not in the predefined map.
    missing_realms <- setdiff(unique_realms, names(realm_map))
    missing_realms <- missing_realms[missing_realms != ""]
    # Assign letters (a, b, c...) to any missing realms.
    extra_map <- setNames(letters[seq_along(missing_realms)], missing_realms)
    # Create the final, complete map.
    full_map <- c(realm_map, extra_map)
    # Apply the mapping.
    tips_realm_mapped <- sapply(tips_realm, function(x) {
      if (x == "") {
        "" # Preserve empty strings.
      } else {
        full_map[[x]]
      }
    })
    tips_realm <- tips_realm_mapped
  }
  
  # Convert the named vector to a data frame for saving.
  tips_df <- data.frame(
    ID = names(tips_realm),
    realm = as.vector(tips_realm),
    stringsAsFactors = FALSE
  )
  write.csv(tips_df, file = tip_states_filepath, row.names = FALSE)
  
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
  
  # Optionally prune non-target families from the tree.
  if (discard_other_family) {
    all_tips <- cs_tree$tip.label
    families <- sapply(all_tips, extract_family)
    target_family <- names(which.max(table(families)))
    # Identify tips that do not belong to the target family (or are NA).
    tips_to_drop <- names(families)[families != target_family | is.na(families)]
    # Ensure all tips to be dropped actually exist in the tree to prevent errors.
    tips_to_drop <- intersect(tips_to_drop, all_tips)
    # Prune the tips.
    cs_tree <- drop.tip(cs_tree, tips_to_drop)
  }
  
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
  if(length(tips_to_delete > 0)) {
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
    } else if (model_name == "DIVALIKE+J") {
      run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
      run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
      run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "free"
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
      d_events <- events_table[events_table$event_type == "d", ]
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
      j_events <- clado_table[clado_table$clado_event_type == "founder (j)", ]
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
  
  # Convert "average event counts" to "average event rates" by dividing by the duration of each time slice.
  durations_Ma <- sapply(my_time_periods, function(p) p[1] - p[2])
  
  avg_d_timeslice_rate_list <- mapply("/", avg_d_timeslice_list, durations_Ma, SIMPLIFY = FALSE)
  avg_j_timeslice_rate_list <- mapply("/", avg_j_timeslice_list, durations_Ma, SIMPLIFY = FALSE)
  avg_total_timeslice_rate_list <- mapply("/", avg_total_timeslice_list, durations_Ma, SIMPLIFY = FALSE)
  
  # Print the final time-stratified rate matrices.
  cat("\n\n--- FINAL TIME-STRATIFIED RESULTS (AS RATES) ---\n\n")
  cat("--- (1) Average D-TYPE Event Rate Matrices per Time Slice (events/Myr) ---\n")
  print(lapply(avg_d_timeslice_rate_list, function(mat) round(mat, 4)))
  
  cat("\n--- (2) Average J-TYPE Event Rate Matrices per Time Slice (events/Myr) ---\n")
  print(lapply(avg_j_timeslice_rate_list, function(mat) round(mat, 4)))
  
  cat("\n--- (3) Average TOTAL (d+j) Event Rate Matrices per Time Slice (events/Myr) ---\n")
  print(lapply(avg_total_timeslice_rate_list, function(mat) round(mat, 4)))
  
  # --- Visualize All Time-Stratified Rate Matrices ---
  
  # A. Global Setup for Plotting
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(grid)
  
  # B. A more flexible and reusable plotting function for heatmaps
  plot_heatmap <- function(mat, subplot_title, scale_mode = "global", color_limits = NULL) {
    
    # --- NEW LOGIC for flexible scaling ---
    current_limits <- NULL
    if (scale_mode == "local") {
      # For local scale, calculate limits from the current matrix 'mat'
      local_max <- if (length(mat) > 0 && max(mat, na.rm = TRUE) > 0) max(mat, na.rm = TRUE) else 1
      current_limits <- c(0, local_max)
    } else {
      # For global scale, use the externally provided 'color_limits'
      current_limits <- color_limits
    }
    
    df <- reshape2::melt(as.matrix(mat))
    colnames(df) <- c("From", "To", "Rate")
    
    ggplot(df, aes(x = To, y = From, fill = Rate)) +
      geom_tile(color = "grey80") +
      geom_text(aes(label = ifelse(Rate > 0, sprintf("%.1e", Rate), "0.0e+00")), size = 3) +
      scale_fill_gradient(
        low = "lightblue", high = "darkred",
        limits = current_limits, # <-- Use the determined limits
        name = "Rate\n(events/Myr)"
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
  
  # =========================================================================
  # Generate and Save All Plots (Both Global and Independent Scales)
  # =========================================================================
  
  # --- C. Calculate the SINGLE, GLOBAL color scale ---
  # This part remains the same.
  all_rates <- unlist(c(avg_d_timeslice_rate_list, avg_j_timeslice_rate_list))
  rate_min <- 0
  rate_max <- if (length(all_rates) > 0 && max(all_rates, na.rm = TRUE) > 0) max(all_rates, na.rm = TRUE) else 1
  global_color_limits <- c(rate_min, rate_max)
  
  
  # --- D. Plot and Save: D-TYPE Event Rates ---
  titles_d <- names(avg_d_timeslice_rate_list)
  
  # D.1: Global Scale for D-events
  plots_d_global <- mapply(plot_heatmap, avg_d_timeslice_rate_list, titles_d, 
                           MoreArgs = list(scale_mode = "global", color_limits = global_color_limits), SIMPLIFY = FALSE)
  plot_grid_d_global <- gridExtra::grid.arrange(
    grobs = plots_d_global, ncol = 3,
    top = grid::textGrob("Time-Stratified Range-Expansion (d-type) Rates (Global Scale)", gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )
  save_path_d_global <- file.path(savedir, paste0(model_name, "_Time-Stratified_rates_D_events_GLOBAL_scale.png"))
  ggsave(filename = save_path_d_global, plot = plot_grid_d_global, width = 20, height = 5 * ceiling(length(plots_d_global)/3), dpi = 300, limitsize = FALSE)
  cat("D-type (Global Scale) plot saved to:", save_path_d_global, "\n")
  
  # D.2: Independent Scale for D-events
  plots_d_local <- mapply(plot_heatmap, avg_d_timeslice_rate_list, titles_d, 
                          MoreArgs = list(scale_mode = "local"), SIMPLIFY = FALSE)
  plot_grid_d_local <- gridExtra::grid.arrange(
    grobs = plots_d_local, ncol = 3,
    top = grid::textGrob("Time-Stratified Range-Expansion (d-type) Rates (Independent Scales)", gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )
  save_path_d_local <- file.path(savedir, paste0(model_name, "_Time-Stratified_rates_D_events_INDEPENDENT_scale.png"))
  ggsave(filename = save_path_d_local, plot = plot_grid_d_local, width = 20, height = 5 * ceiling(length(plots_d_local)/3), dpi = 300, limitsize = FALSE)
  cat("D-type (Independent Scale) plot saved to:", save_path_d_local, "\n")
  
  
  # --- E. Plot and Save: J-TYPE Event Rates ---
  if (max(unlist(avg_j_timeslice_rate_list), na.rm = TRUE) > 0) {
    titles_j <- names(avg_j_timeslice_rate_list)
    
    # E.1: Global Scale for J-events
    plots_j_global <- mapply(plot_heatmap, avg_j_timeslice_rate_list, titles_j, 
                             MoreArgs = list(scale_mode = "global", color_limits = global_color_limits), SIMPLIFY = FALSE)
    plot_grid_j_global <- gridExtra::grid.arrange(
      grobs = plots_j_global, ncol = 3,
      top = grid::textGrob("Time-Stratified Founder-Event (j-type) Rates (Global Scale)", gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    save_path_j_global <- file.path(savedir, paste0(model_name, "_Time-Stratified_rates_J_events_GLOBAL_scale.png"))
    ggsave(filename = save_path_j_global, plot = plot_grid_j_global, width = 20, height = 5 * ceiling(length(plots_j_global)/3), dpi = 300, limitsize = FALSE)
    cat("J-type (Global Scale) plot saved to:", save_path_j_global, "\n")
    
    # E.2: Independent Scale for J-events
    plots_j_local <- mapply(plot_heatmap, avg_j_timeslice_rate_list, titles_j, 
                            MoreArgs = list(scale_mode = "local"), SIMPLIFY = FALSE)
    plot_grid_j_local <- gridExtra::grid.arrange(
      grobs = plots_j_local, ncol = 3,
      top = grid::textGrob("Time-Stratified Founder-Event (j-type) Rates (Independent Scales)", gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    save_path_j_local <- file.path(savedir, paste0(model_name, "_Time-Stratified_rates_J_events_INDEPENDENT_scale.png"))
    ggsave(filename = save_path_j_local, plot = plot_grid_j_local, width = 20, height = 5 * ceiling(length(plots_j_local)/3), dpi = 300, limitsize = FALSE)
    cat("J-type (Independent Scale) plot saved to:", save_path_j_local, "\n")
    
  } else {
    cat("All J-type event rates are zero. Plotting was skipped for J-events.\n")
  }
  
  
  # --- F. Plot and Save: TOTAL (d+j) Event Rates ---
  titles_total <- names(avg_total_timeslice_rate_list)
  
  # F.1: Global Scale for Total-events
  plots_total_global <- mapply(plot_heatmap, avg_total_timeslice_rate_list, titles_total, 
                               MoreArgs = list(scale_mode = "global", color_limits = global_color_limits), SIMPLIFY = FALSE)
  plot_grid_total_global <- gridExtra::grid.arrange(
    grobs = plots_total_global, ncol = 3,
    top = grid::textGrob("Time-Stratified Total Dispersal (d+j) Rates (Global Scale)", gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )
  save_path_total_global <- file.path(savedir, paste0(model_name, "_Time-Stratified_rates_TOTAL_dispersal_GLOBAL_scale.png"))
  ggsave(filename = save_path_total_global, plot = plot_grid_total_global, width = 20, height = 5 * ceiling(length(plots_total_global)/3), dpi = 300, limitsize = FALSE)
  cat("Total Dispersal (Global Scale) plot saved to:", save_path_total_global, "\n")
  
  # F.2: Independent Scale for Total-events
  plots_total_local <- mapply(plot_heatmap, avg_total_timeslice_rate_list, titles_total, 
                              MoreArgs = list(scale_mode = "local"), SIMPLIFY = FALSE)
  plot_grid_total_local <- gridExtra::grid.arrange(
    grobs = plots_total_local, ncol = 3,
    top = grid::textGrob("Time-Stratified Total Dispersal (d+j) Rates (Independent Scales)", gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )
  save_path_total_local <- file.path(savedir, paste0(model_name, "_Time-Stratified_rates_TOTAL_dispersal_INDEPENDENT_scale.png"))
  ggsave(filename = save_path_total_local, plot = plot_grid_total_local, width = 20, height = 5 * ceiling(length(plots_total_local)/3), dpi = 300, limitsize = FALSE)
  cat("Total Dispersal (Independent Scale) plot saved to:", save_path_total_local, "\n")
  
  # --- 7.6 Generate all temporal trend plots ---
  cat("\n--- Generating temporal trend plots ---\n")
  
  tree_age <- max(branching.times(tr))
  binning_size <- 10 # Set a uniform bin size for time-based event counting.
  
  # --- Main Figure ---
  main_plot <- plot_immigration_emigration_trends(
    clado_events_tables = clado_events_tables_src,
    ana_events_tables = ana_events_tables_src,
    areanames = areanames, tree_age = tree_age, bin_size = binning_size
  )
  ggsave(file.path(savedir, paste0(model_name, "_main_Immigration_vs_Emigration.png")), 
         plot = main_plot, width = 12, height = 8, dpi = 300)
  cat("Main trend plot (Immigration vs. Emigration) saved.\n")
  
  # --- Supplementary Figure 1 ---
  supp_plot1 <- plot_immigration_components(
    clado_events_tables = clado_events_tables_src,
    ana_events_tables = ana_events_tables_src,
    areanames = areanames, tree_age = tree_age, bin_size = binning_size
  )
  ggsave(file.path(savedir, paste0(model_name, "_supp_Immigration_Components.png")), 
         plot = supp_plot1, width = 12, height = 8, dpi = 300)
  cat("Supplementary plot (Immigration Components d vs j) saved.\n")
  
  # --- Supplementary Figure 2 ---
  supp_plot2 <- plot_emigration_components(
    clado_events_tables = clado_events_tables_src,
    ana_events_tables = ana_events_tables_src,
    areanames = areanames, tree_age = tree_age, bin_size = binning_size
  )
  ggsave(file.path(savedir, paste0(model_name, "_supp_Emigration_Components.png")), 
         plot = supp_plot2, width = 12, height = 8, dpi = 300)
  cat("Supplementary plot (Emigration Components d vs j) saved.\n")
  

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
                                               bin_size = 10) {
  
  # Define time bins based on tree age and desired bin size.
  time_breaks <- seq(0, ceiling(tree_age / bin_size) * bin_size, by = bin_size)
  time_labels <- paste(time_breaks[-length(time_breaks)], time_breaks[-1], sep="-")
  num_sims <- length(clado_events_tables)
  
  # Iterate through each BSM simulation to tally events.
  results_per_sim <- list()
  for (i in 1:num_sims) {
    clado_table <- clado_events_tables[[i]]
    ana_table <- ana_events_tables[[i]]
    
    # a. Tally Immigration events (d_in + j_in).
    # Immigration is defined as any dispersal event where a region is the DESTINATION.
    d_in <- if(nrow(ana_table) > 0) ana_table %>% filter(event_type == "d") %>% select(time=abs_event_time, area=dispersal_to) else data.frame()
    j_in <- if(nrow(clado_table) > 0) clado_table %>% filter(clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_to) else data.frame()
    imm_counts <- bind_rows(d_in, j_in) %>%
      mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
      filter(!is.na(time_bin)) %>%
      group_by(time_bin, area) %>% summarise(count = n(), .groups = 'drop') %>%
      mutate(event_type = "Immigration")
    
    # b. Tally Emigration events (d_out + j_out).
    # Emigration is defined as any dispersal event where a region is the SOURCE.
    d_out <- if(nrow(ana_table) > 0) ana_table %>% filter(event_type == "d") %>% select(time=abs_event_time, area=ana_dispersal_from) else data.frame()
    j_out <- if(nrow(clado_table) > 0) clado_table %>% filter(clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_from) else data.frame()
    emi_counts <- bind_rows(d_out, j_out) %>%
      mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
      filter(!is.na(time_bin) & !is.na(area)) %>% # Source area cannot be NA.
      group_by(time_bin, area) %>% summarise(count = n(), .groups = 'drop') %>%
      mutate(event_type = "Emigration")
    
    # Combine results for the current simulation.
    sim_results <- bind_rows(imm_counts, emi_counts)
    if(nrow(sim_results) > 0) {
      sim_results$sim_id <- i
      results_per_sim[[i]] <- sim_results
    }
  }
  
  # Aggregate results from all simulations.
  all_sim_results <- bind_rows(results_per_sim)
  
  # Create a full grid of all possible combinations to handle zero-count cases correctly.
  full_grid <- expand.grid(time_bin = time_labels, area = areanames, 
                           event_type = c("Immigration", "Emigration"),
                           sim_id = 1:num_sims, stringsAsFactors = FALSE)
  
  # Join simulation results to the full grid and replace NAs with 0.
  summary_data <- full_grid %>%
    left_join(all_sim_results, by = c("time_bin", "area", "event_type", "sim_id")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  # Calculate the mean and 5%-95% quantiles for each group.
  final_summary <- summary_data %>%
    group_by(time_bin, area, event_type) %>%
    summarise(mean_events = mean(count),
              lower_q = quantile(count, 0.05),
              upper_q = quantile(count, 0.95), .groups = 'drop') %>%
    # Convert time bin factor to a numeric value (midpoint) for plotting.
    mutate(time_mid = as.numeric(sub("-.*", "", time_bin)) + (bin_size / 2))
  
  # Generate the plot using ggplot2.
  p <- ggplot(final_summary, aes(x = time_mid, y = mean_events, color = event_type, fill = event_type)) +
    geom_ribbon(aes(ymin = lower_q, ymax = upper_q), alpha = 0.2, linetype = 0) + # Shaded 90% quantile interval
    geom_line(size = 1) + # Mean trend line
    facet_wrap(~ area, scales = "fixed") + # Facet by area with a fixed y-axis for direct comparison
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
                                        bin_size = 10) {
  
  time_breaks <- seq(0, ceiling(tree_age / bin_size) * bin_size, by = bin_size)
  time_labels <- paste(time_breaks[-length(time_breaks)], time_breaks[-1], sep="-")
  num_sims <- length(clado_events_tables)
  
  results_per_sim <- list()
  for (i in 1:num_sims) {
    clado_table <- clado_events_tables[[i]]
    ana_table <- ana_events_tables[[i]]
    
    # Tally d-type and j-type immigration events separately.
    d_in <- if(nrow(ana_table) > 0) ana_table %>% filter(event_type == "d") %>% select(time=abs_event_time, area=dispersal_to) %>% mutate(event_type="d (Range Expansion)") else data.frame()
    j_in <- if(nrow(clado_table) > 0) clado_table %>% filter(clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_to) %>% mutate(event_type="j (Founder Event)") else data.frame()
    
    sim_results <- bind_rows(d_in, j_in) %>%
      mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
      filter(!is.na(time_bin)) %>%
      group_by(time_bin, area, event_type) %>% summarise(count = n(), .groups = 'drop')
    
    if(nrow(sim_results) > 0) {
      sim_results$sim_id <- i
      results_per_sim[[i]] <- sim_results
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
  
  p <- ggplot(final_summary, aes(x = time_mid, y = mean_events, color = event_type, fill = event_type)) +
    geom_ribbon(aes(ymin = lower_q, ymax = upper_q), alpha = 0.2, linetype = 0) +
    geom_line(size = 1) +
    facet_wrap(~ area, scales = "free_y") + # Using a free y-axis is often better for component plots
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
                                       bin_size = 10) {
  
  time_breaks <- seq(0, ceiling(tree_age / bin_size) * bin_size, by = bin_size)
  time_labels <- paste(time_breaks[-length(time_breaks)], time_breaks[-1], sep="-")
  num_sims <- length(clado_events_tables)
  
  results_per_sim <- list()
  for (i in 1:num_sims) {
    clado_table <- clado_events_tables[[i]]
    ana_table <- ana_events_tables[[i]]
    
    # Tally d-type and j-type emigration events separately.
    d_out <- if(nrow(ana_table) > 0) ana_table %>% filter(event_type == "d") %>% select(time=abs_event_time, area=ana_dispersal_from) %>% mutate(event_type="d (Range Expansion)") else data.frame()
    j_out <- if(nrow(clado_table) > 0) clado_table %>% filter(clado_event_type == "founder (j)") %>% select(time=time_bp, area=clado_dispersal_from) %>% mutate(event_type="j (Founder Event)") else data.frame()
    
    sim_results <- bind_rows(d_out, j_out) %>%
      mutate(time_bin = cut(time, breaks = time_breaks, labels = time_labels, right = FALSE)) %>%
      filter(!is.na(time_bin) & !is.na(area)) %>%
      group_by(time_bin, area, event_type) %>% summarise(count = n(), .groups = 'drop')
    
    if(nrow(sim_results) > 0) {
      sim_results$sim_id <- i
      results_per_sim[[i]] <- sim_results
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
  
  p <- ggplot(final_summary, aes(x = time_mid, y = mean_events, color = event_type, fill = event_type)) +
    geom_ribbon(aes(ymin = lower_q, ymax = upper_q), alpha = 0.2, linetype = 0) +
    geom_line(size = 1) +
    facet_wrap(~ area, scales = "free_y") + # Using a free y-axis is often better for component plots
    labs(title = "Components of Emigration Events by Region",
         subtitle = paste("Averaged over", num_sims, "BSM scenarios in", bin_size, "Ma bins"),
         x = "Time (Million years ago)", y = "Average Number of Events per Time Bin",
         color = "Emigration Type", fill = "Emigration Type") +
    scale_x_reverse() + theme_bw() + theme(legend.position = "bottom")
  
  return(p)
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
  
  
  
  
  
  pastml_process_tree(dated_tree_path = dated_tree_path,
                      discard_other_family = TRUE
  )
  result_info <- run_biogeobears_pipeline(
    tree_filepath = dated_tree_path,
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

