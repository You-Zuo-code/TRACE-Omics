# Run from the app root; babelgene is needed only to rebuild this snapshot.
stopifnot(as.character(utils::packageVersion("babelgene")) == "22.9")
d <- get("orthologs_df", asNamespace("babelgene"))
d <- d[d$taxon_id %in% c(9823L, 10090L, 10116L), ]
d$organism <- c(`9823` = "Pig", `10090` = "Mouse", `10116` = "Rat")[as.character(d$taxon_id)]
d$relationship <- "ambiguous"
for (sp in unique(d$organism)) {
  ix <- which(d$organism == sp)
  z <- d[ix, ]
  pairs <- unique(z[, c("entrez", "human_entrez")])
  pairs <- pairs[complete.cases(pairs), ]
  forward <- table(pairs$entrez); reverse <- table(pairs$human_entrez)
  ok <- !is.na(z$entrez) & !is.na(z$human_entrez) &
    forward[as.character(z$entrez)] == 1L & reverse[as.character(z$human_entrez)] == 1L
  d$relationship[ix[which(ok)]] <- "one_to_one_in_HCOP_snapshot"
}
saveRDS(d, "data/species_orthologs.rds", compress = "xz", version = 3)
meta <- list(source = "HGNC HCOP via babelgene", version = "22.9",
  url = "https://igordot.github.io/babelgene/", built = as.character(Sys.Date()),
  species = c(Human = 9606, Pig = 9823, Mouse = 10090, Rat = 10116),
  rows = nrow(d), sha256 = digest::digest(file = "data/species_orthologs.rds", algo = "sha256"),
  rule = paste("Bidirectional unique Entrez-gene pairs in the complete species snapshot,",
               "before any input or support filtering; at least three supporting resources.",
               "This is snapshot uniqueness, not an Ensembl homology classification."))
jsonlite::write_json(meta, "data/species_orthologs_provenance.json", pretty = TRUE, auto_unbox = TRUE)
print(table(d$organism, d$relationship))
