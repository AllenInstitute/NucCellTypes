load("//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Nuclei_Noisemodel/R_objects/samp.dat_MouseNuclei6960_star.Rdata")
load("C:/Users/trygveb/Documents/Projects/Cell_vs_Nuc/cache/all_layers/map_nuc_to_cells_all_layers.rda")


samp.ids <- names(max.cor.nuc[(max.cor.nuc - max.cor.all) > 0.2])
samp.dat[row.names(samp.dat) %in% samp.ids, ]

samp.excl <- as.factor(samp.dat$Class[match(names(max.cor.nuc), row.names(samp.dat))])

plot(max.cor.nuc, max.cor.all, col = as.numeric(samp.excl), cex = 0.5)
legend("bottomright", fill = 1:nlevels(samp.excl), 
       legend = levels(samp.excl))
abline(a = 0, b = 1)
