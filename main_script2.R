# Windows powershell compatibility
options(
  error = function() {
    traceback(2)
    flush.console()
    q(status = 1, save = "no")
  }
)


cat("[", Sys.time(), "] Beginning main_script2.R\n")



set.seed(17)
options(future.globals.maxSize = 2 * 1024^3)  # 2GB

global_dir <- "."
gtf_file    <- file.path(global_dir, "FB_dmel-all-r6.58.gtf")
source("aux2.R")



# (older version included more parameters and values)
datasets <- c("5D_GAA","10D_GAA","10D_-AA","18D_GAA")
params <- list(
  lower_bound = 100,
  retain      = 500,
  significance_threshold = 0.01,
  nmads       = 3,
  propRandom  = 0.5
)

# Define unwanted features (e.g. mitochondrial)
# This should be computed per dataset, as genes captured may vary
negative_features <- NULL  # assign as needed

# Marker genes for UMAP 
markers <- c("vas","bam","nanos","osk","orb","piwi","pgc","tud",
             "tj","Fas3","ct","br","bab1","bab2","eya",
             "grk","bcd","stau","tor",
             "Delta","N","hh","Wnt4","dpp","foxo",
             "Ilp8","Ilp6",
             "ovo","slbo","chinmo")

for (ds in datasets) {

  cat("\n=== Processing", ds, "===\n")
  path <- file.path(global_dir, ds)
  out1 <- file.path(path, "output_relaxed_new_final_may_33")
  out2 <- file.path(path, "output_stringent_new_final_may_33")
  dir.create(out1, recursive=TRUE, showWarnings=FALSE)
  dir.create(out2, recursive=TRUE, showWarnings=FALSE)

  sce0 <- load_data(path, gtf_file)

  is.mito <- grep("^MT:", rownames(sce0), value = TRUE, ignore.case = TRUE)
  is.ribo <- grep("^RPS|^RPL", rownames(sce0), value = TRUE, ignore.case = TRUE)
  is.rribo <- grep("rRNA", rownames(sce0), value = TRUE, ignore.case = FALSE)
  is.tribo <- grep("tRNA", rownames(sce0), value = TRUE, ignore.case = FALSE)
  
  is.rribo <- setdiff(is.rribo, is.mito)
  is.tribo <- setdiff(is.tribo, is.mito)
  
  all_ribo <- unique(c(is.rribo, is.ribo, is.tribo))
  
  is.prerribo <- grep("pre", is.rribo, value = TRUE, ignore.case = TRUE)
  non_pre_ribo <- setdiff(all_ribo, is.prerribo)
  
  negative_features <- unique(c(is.mito, is.tribo, is.rribo))
  negative_features <- setdiff(negative_features, is.prerribo)



  sce1 <- empty_droplet(sce0, params, negative_features, out1, ds)
  sce2 <- QC(sce1, "counts", params, out1, ds)
  
  sce3 <- gene_n_doublet(sce2, params, out1, ds)
  sce4 <- vst_norm_dim_red(sce3, params, out1, ds)
  marker_plotting(sce4, markers, out1, ds)
  saveRDS(sce4, file.path(out1, paste0(ds, "_final_no_decont.rds")))
  export_h5ad(sce4, out1, ds)

  sce_d <- decontamination(sce1, sce0, params, out2, ds)
  sce_qc2 <- QC(sce_d, "decontXcounts", params, out2, ds)
  sce_dbl2 <- gene_n_doublet(sce_qc2, params, out2, ds)
  sce_vst2 <- vst_norm_dim_red(sce_dbl2, params, out2, ds)
  marker_plotting(sce_vst2, markers, out2, ds)
  saveRDS(sce_vst2, file.path(out2, paste0(ds, "_final_with_decont.rds")))
  export_h5ad(sce_vst2, out2, ds)
}

cat("\nAll datasets processed successfully!\n")
