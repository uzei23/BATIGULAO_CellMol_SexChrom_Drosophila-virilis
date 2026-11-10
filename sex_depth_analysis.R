
# ============================================================
# DROSOPHILA VIRILIS SEX-CHROMOSOME DOSAGE ANALYSIS
# ============================================================

# 1. SELECT BOTH FILES IN ANY ORDER
cat("Select your Galaxy coverage file OR chromosome_groups.csv.\n")
file1 <- file.choose()

cat("Select the remaining file.\n")
file2 <- file.choose()

# 2. IDENTIFY THE FILES AUTOMATICALLY
read_first_line <- function(path) {
  readLines(path, n = 1, warn = FALSE)
}

first1 <- read_first_line(file1)
first2 <- read_first_line(file2)

if (grepl("meandepth", first1, ignore.case = TRUE)) {
  coverage_path <- file1
  groups_path <- file2
} else if (grepl("meandepth", first2, ignore.case = TRUE)) {
  coverage_path <- file2
  groups_path <- file1
} else {
  stop("Cannot identify the Galaxy coverage file. Check that you selected the .tabular file.")
}

# 3. CREATE RESULTS FOLDER
results_dir <- file.path(dirname(groups_path), "results")
dir.create(results_dir, showWarnings = FALSE)

# 4. READ COVERAGE FILE
coverage <- read.delim(
  coverage_path,
  header = TRUE,
  sep = "\t",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(coverage) <- trimws(names(coverage))
names(coverage) <- sub("^#", "", names(coverage))

if (!all(c("rname", "meandepth") %in% names(coverage))) {
  stop(paste("Incorrect coverage columns:", paste(names(coverage), collapse = ", ")))
}

# 5. READ CHROMOSOME GROUPS CSV
groups <- read.csv(
  groups_path,
  stringsAsFactors = FALSE
)

names(groups) <- trimws(names(groups))

if (!all(c("rname", "group") %in% names(groups))) {
  stop("Chromosome CSV must contain rname and group columns.")
}

# 6. MATCH CHROMOSOMES
coverage$meandepth <- as.numeric(coverage$meandepth)

analysis_data <- merge(
  groups,
  coverage[, c("rname", "meandepth")],
  by = "rname"
)

if (nrow(analysis_data) != nrow(groups)) {
  stop("Chromosome names do not match.")
}

analysis_data <- analysis_data[
  match(groups$rname, analysis_data$rname),
]

analysis_data$chromosome <- c("X", "2", "3", "4", "5", "6")

autosomes <- analysis_data[analysis_data$group == "A", ]
x_chromosome <- analysis_data[analysis_data$group == "X", ]

if (nrow(autosomes) != 5 || nrow(x_chromosome) != 1) {
  stop("Expected five autosomes and one X chromosome.")
}

# 7. CALCULATE DEPTH
autosomal_mean <- mean(autosomes$meandepth)
x_depth <- x_chromosome$meandepth
x_to_a_ratio <- x_depth / autosomal_mean

analysis_data$relative_depth <-
  analysis_data$meandepth / autosomal_mean

# 8. SENSITIVITY CHECK
autosomes_without_chr6 <- autosomes[
  autosomes$rname != "NC_091548.1", ]

autosomal_mean_without_chr6 <-
  mean(autosomes_without_chr6$meandepth)

x_to_a_without_chr6 <-
  x_depth / autosomal_mean_without_chr6

# 9. INTERPRETATION
if (x_to_a_ratio <= 0.65) {
  inference <- "Consistent with one X chromosome relative to two autosomal copies."
} else if (x_to_a_ratio >= 0.85) {
  inference <- "Consistent with two X chromosomes relative to two autosomal copies."
} else {
  inference <- "Ambiguous X-chromosome dosage."
}

# 10. SAVE TABLE
write.csv(
  analysis_data,
  file.path(results_dir, "chromosome_relative_depth.csv"),
  row.names = FALSE
)

# 11. SAVE GRAPH
png(
  file.path(results_dir, "chromosome_depth.png"),
  width = 1200,
  height = 800,
  res = 120
)

barplot(
  analysis_data$relative_depth,
  names.arg = analysis_data$chromosome,
  col = ifelse(analysis_data$group == "X", "steelblue", "gray70"),
  main = "Relative Chromosome Depth in Drosophila virilis",
  xlab = "Chromosome",
  ylab = "Relative Depth (Autosomal Mean = 1)",
  ylim = c(0, max(analysis_data$relative_depth) * 1.2)
)

abline(h = 1, lty = 2, col = "red")
abline(h = 0.5, lty = 2, col = "blue")

legend(
  "topright",
  legend = c("X chromosome", "Autosomes", "Expected 0.5", "Expected 1.0"),
  fill = c("steelblue", "gray70", NA, NA),
  lty = c(NA, NA, 2, 2),
  col = c(NA, NA, "blue", "red"),
  bty = "n"
)

dev.off()

# 12. SAVE REPORT
sink(file.path(results_dir, "sex_inference.txt"))

cat("SEX-CHROMOSOME DOSAGE ANALYSIS\n")
cat("Species: Drosophila virilis\n\n")
cat("WGS Run: SRR6399436\n")
cat("WGS BioSample: SAMN08216748\n")
cat("Recorded WGS sample sex: Male\n")
cat("Strain: 15010-1051.87\n\n")

cat("Assembly: GCF_030788295.1\n")
cat("Assembly BioSample: SAMN35709873\n")
cat("Recorded assembly donor sex: Male\n")
cat("The assembly and WGS BioSamples are different.\n\n")

cat("CHROMOSOME DEPTH RESULTS\n")
print(analysis_data[, c(
  "chromosome", "group", "meandepth", "relative_depth"
)])

cat("\nMean autosomal depth:", autosomal_mean, "\n")
cat("X chromosome depth:", x_depth, "\n")
cat("X:A ratio:", x_to_a_ratio, "\n")

cat("\nSENSITIVITY CHECK EXCLUDING CHROMOSOME 6\n")
cat("X:A ratio:", x_to_a_without_chr6, "\n")

cat("\nINTERPRETATION\n")
cat(inference, "\n")

cat("\nBIOLOGICAL NOTE\n")
cat("Drosophila sex determination depends on the X-to-autosome ratio.\n")
cat("The Y chromosome is primarily required for male fertility.\n")

cat("\nLIMITATIONS\n")
cat("1. Assembly and WGS BioSamples are different.\n")
cat("2. Mapping bias and repetitive sequences may affect read depth.\n")
cat("3. Chromosome 6 has lower coverage than the other autosomes.\n")
cat("4. Relative depth does not independently confirm biological sex.\n")
cat("5. Interpretation thresholds are simplified classroom guidelines.\n")

sink()

# 13. DISPLAY RESULTS
cat("\n====================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("====================================\n")

cat("X:A ratio:", round(x_to_a_ratio, 3), "\n")
cat("X:A ratio excluding chromosome 6:",
    round(x_to_a_without_chr6, 3), "\n")
cat("Interpretation:", inference, "\n")
cat("Results folder:", results_dir, "\n")

results_dir <- "C:/Users/Uzei/Documents/Sex_Dosage_Drosophila_virilis/data/results"

cat("CHECKING RESULTS FOLDER...\n\n")

cat("Folder exists:", dir.exists(results_dir), "\n\n")

cat("Files saved:\n")
print(list.files(results_dir, full.names = TRUE))

cat("\nFile sizes:\n")
print(file.info(list.files(results_dir, full.names = TRUE))$size)



# 1. SET FILE PATHS
project_dir <- "C:/Users/Uzei/Documents/Sex_Dosage_Drosophila_virilis"

coverage_path <- file.path(
  project_dir, "data",
  "Galaxy16-[Samtools coverage on dataset 15].tabular"
)

groups_path <- file.path(
  project_dir, "data",
  "chromosome_groups.csv"
)

results_dir <- file.path(project_dir, "data", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# 2. READ THE TWO INPUT FILES
coverage <- read.delim(
  coverage_path,
  header = TRUE,
  sep = "\t",
  comment.char = "",
  check.names = FALSE
)

names(coverage) <- sub("^#", "", trimws(names(coverage)))

groups <- read.csv(groups_path, stringsAsFactors = FALSE)

# 3. MATCH CHROMOSOMES
coverage$meandepth <- as.numeric(coverage$meandepth)

analysis_data <- merge(
  groups,
  coverage[, c("rname", "meandepth")],
  by = "rname"
)

analysis_data <- analysis_data[
  match(groups$rname, analysis_data$rname),
]

analysis_data$chromosome <- c("X", "2", "3", "4", "5", "6")

# 4. CALCULATE RELATIVE DEPTH
autosomes <- analysis_data[analysis_data$group == "A", ]
x_chromosome <- analysis_data[analysis_data$group == "X", ]

autosomal_mean <- mean(autosomes$meandepth)
analysis_data$relative_depth <- analysis_data$meandepth / autosomal_mean

x_depth <- x_chromosome$meandepth
x_to_a_ratio <- x_depth / autosomal_mean

# 5. SAVE CSV
csv_path <- file.path(results_dir, "chromosome_relative_depth.csv")

write.csv(
  analysis_data,
  csv_path,
  row.names = FALSE
)

# 6. CREATE GRAPH
png(
  filename = file.path(results_dir, "chromosome_depth.png"),
  width = 1200,
  height = 800,
  res = 120,
  bg = "white"
)

par(mar = c(5, 5, 4, 2))

barplot(
  height = analysis_data$relative_depth,
  names.arg = analysis_data$chromosome,
  col = ifelse(analysis_data$group == "X", "steelblue", "gray70"),
  main = "Relative Chromosome Depth in Drosophila virilis",
  xlab = "Chromosome",
  ylab = "Relative Depth (Autosomal Mean = 1)",
  ylim = c(0, max(analysis_data$relative_depth) * 1.2)
)

abline(h = 1, lty = 2, col = "red", lwd = 2)
abline(h = 0.5, lty = 2, col = "blue", lwd = 2)

legend(
  "topright",
  legend = c("X chromosome", "Autosomes", "Expected 0.5", "Expected 1.0"),
  fill = c("steelblue", "gray70", NA, NA),
  lty = c(NA, NA, 2, 2),
  col = c(NA, NA, "blue", "red"),
  bty = "n"
)

dev.off()

# 7. VERIFY SAVED FILES
cat("\nCHROMOSOME DEPTH RESULTS\n")
print(analysis_data[, c("chromosome", "group", "meandepth", "relative_depth")])

cat("\nX:A ratio:", round(x_to_a_ratio, 3), "\n")

cat("\nSAVED FILES:\n")
print(file.info(c(
  csv_path,
  file.path(results_dir, "chromosome_depth.png")
))[, "size"])

# 8. OPEN RESULTS FOLDER
shell.exec(normalizePath(results_dir))
