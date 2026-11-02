#!/usr/bin/env Rscript
# Run one of the workflow R scripts outside Snakemake, for testing.
# Usage: Rscript run_script_standalone.R <script.R> <spec.json>
# spec.json: {"input": {...}, "output": {...}, "params": {...}, "log": ["file"], "threads": 4, "config": {...}}
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
spec <- fromJSON(args[2], simplifyVector = TRUE)
setClass("Snakemake", representation(input = "list", output = "list", params = "list", log = "list",
                                     threads = "numeric", config = "list", wildcards = "list"))
as_list <- function(x) if (is.null(x)) list() else as.list(x)
snakemake <- new("Snakemake", input = as_list(spec$input), output = as_list(spec$output), params = as_list(spec$params),
                 log = as_list(spec$log), threads = if (is.null(spec$threads)) 1 else spec$threads,
                 config = as_list(spec$config), wildcards = list())
for (p in unlist(spec$output)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
for (p in unlist(spec$log)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
source(args[1])
