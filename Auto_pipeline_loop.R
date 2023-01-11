library(flowCore)
library(flowDensity)
library(flowViz)
library(flowWorkspace)
cwd_0 <- "/Users/mgj1343/Automatization-pipeline"
setwd(cwd_0)
samples <- read.csv(file = "/Users/mgj1343/Library/CloudStorage/OneDrive-NorthwesternUniversity/SBratchikov/FlowDensity/01_data/PASC/02_BAL_flow_cytometry/20220625_PASC_flow.csv", sep = ",")[, 1]
samples <- samples[1:39]
newdir <- "PASC_pipeline"
dir.create(newdir)
setwd("./PASC_pipeline/")
# fcs_names <- dir(path = "/Users/mgj1343/Library/CloudStorage/OneDrive-NorthwesternUniversity/SBratchikov/FlowDensity/01_data/SCRIPT/2021_06", pattern = ".fcs$", full.names = FALSE, ignore.case = TRUE)

for (i in samples) {
  newdir <- i
  dir.create(newdir)
  # setwd(paste0(cwd,'/',newdir))
  setwd(newdir)
  fcs_test <- read.FCS(paste0("/Users/mgj1343/Library/CloudStorage/OneDrive-NorthwesternUniversity/SBratchikov/FlowDensity/01_data/PASC/02_BAL_flow_cytometry/all_fcs/", i))
  ### singlets
  sngl <- flowDensity(fcs_test,
    channels = c("FSC-A", "FSC-H"), position = c(F, F),
    percentile = c(.99999, .99999), use.percentile = c(T, T),
    ellip.gate = T, scale = .8
  )
  bmp(
    file = "Singlets.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(fcs_test, c(1, 2), main = "Singlets")
  lines(sngl@filter, type = "l")
  legend("topleft", legend = c(
    paste0("count: ", sngl@cell.count),
    paste0("frequency: ", round(sngl@proportion, digits = 2)),
    paste0("total: ", nrow(sngl@flow.frame))
  ), bty = "n")
  dev.off()
  ### logicle transformation
  df_colnames <- setdiff(df_allcolnames, c("Time", "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W"))
  lgcl <- estimateLogicle(fcs_test,channels = df_colnames)
  after <- transform(fcs_test,lgcl)
  ### cd45+
  cd45 <- flowDensity(after[sngl@index], channels = c("AmCyan-A", "Qdot 655-A"), position = c(T, NA), upper = c(F, NA))
  # not_cd45 <- flowDensity(after, channels=c("AmCyan-A", "Qdot 655-A"), position=c(F,NA),upper=c(F,NA))
  bmp(
    file = "CD45+.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[sngl@index], c("AmCyan-A", "Qdot 655-A"), main = "CD45+")
  # lines(not_cd45@filter,type="l",col=2,lwd=4)
  lines(cd45@filter, type = "l")
  legend("topleft", legend = c(
    paste0("count: ", cd45@cell.count),
    paste0("frequency: ", round(cd45@proportion, digits = 2)),
    paste0("total: ", sngl@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()


  ### Live singlet cells
  live <- flowDensity(after[cd45@index], channels = c("FITC-A", "FSC-A"), position = c(F, NA), upper = c(NA, F))
  bmp(
    file = "Live.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[cd45@index], c("FITC-A", "FSC-A"), main = "Live")
  legend("topleft", legend = c(
    paste0("count: ", live@cell.count),
    paste0("frequency: ", round(live@proportion, digits = 2)),
    paste0("total: ", cd45@cell.count)
  ), bty = "n", text.col = 1)
  lines(live@filter, type = "l")
  dev.off()
  ### cd3+

  cd3 <- flowDensity(after[live@index], channels = c("PE-A", "PE-Cy7-A"), position = c(T, F))
  ### getting gates thresholds for markers
  cd3.gate <- deGate(after[live@index], channel = "PE-A")
  cd206.gate <- deGate(after[live@index],channel = "PE-Cy7-A")
  ### selecting the cd3 surrounding cells
  not_cd3 <- notSubFrame(after[live@index], channels = c("PE-A", "PE-Cy7-A"),position = c(T,F),gates=c(cd3.gate,cd206.gate))
  
  bmp(
    file = "CD3 T cells.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[live@index], c("PE-A", "PE-Cy7-A"), main = "CD3 T cells")
  lines(not_cd3@filter, type = "l", col = 2)
  lines(cd3@filter, type = "l")
  # lines(not_cd3@filter,type="l")
  legend("topleft", legend = c(
    paste0("count: ", cd3@cell.count),
    paste0("frequency: ", round(cd3@proportion, digits = 2)),
    paste0("total: ", live@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### CD8 CD4
  cd4 <- flowDensity(after[cd3@index], c("DAPI-A", "APC-A"), position = c(T, F))
  cd8 <- flowDensity(after[cd3@index], c("DAPI-A", "APC-A"), position = c(F, T))
  # not_cd45 <- flowDensity(after, channels=c("AmCyan-A", "Qdot 655-A"), position=c(F,NA),upper=c(F,NA))
  bmp(
    file = "CD8_CD4.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[cd3@index], c("DAPI-A", "APC-A"), main = "CD8/CD4")
  lines(cd8@filter, type = "l", col = 2)
  lines(cd4@filter, type = "l", col = 3)
  legend("topleft", legend = c(
    paste0("count: ", cd8@cell.count),
    paste0("frequency: ", round(cd8@proportion, digits = 2)),
    paste0("total: ", cd3@cell.count)
  ), bty = "n", text.col = 2)
  legend("bottomright", legend = c(
    paste0("count: ", cd4@cell.count),
    paste0("frequency: ", round(cd4@proportion, digits = 2)),
    paste0("total: ", cd3@cell.count)
  ), bty = "n", text.col = 3)
  dev.off()
  ### tregs
  tregs <- flowDensity(after[cd4@index], channels = c("Side Pop-A", "PE-Texas Red-A"), position = c(T, F))
  # not_cd45 <- flowDensity(after, channels=c("AmCyan-A", "Qdot 655-A"), position=c(F,NA),upper=c(F,NA))
  bmp(
    file = "Tregs.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[cd4@index], c("Side Pop-A", "PE-Texas Red-A"), main = "Tregs")
  # lines(not_cd45@filter,type="l",col=2,lwd=4)
  lines(tregs@filter, type = "l")
  legend("topleft", legend = c(
    paste0("count: ", tregs@cell.count),
    paste0("frequency: ", round(tregs@proportion, digits = 2)),
    paste0("total: ", cd4@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### Neutrophils
  neutro <- flowDensity(after[not_cd3@index], channels = c("Qdot 655-A", "Pacific Blue-A"), position = c(T, NA))
  not_neutro <- flowDensity(after[not_cd3@index], channels = c("Qdot 655-A", "Pacific Blue-A"), position = c(F, NA))
  bmp(
    file = "Neutrophils.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[not_cd3@index], c("Side Pop-A", "PE-Texas Red-A"), main = "Tregs")
  plotDens(after[not_cd3@index], c("Qdot 655-A", "Pacific Blue-A"), main = "Neutrophils")
  # lines(not_cd45@filter,type="l",col=2,lwd=4)
  lines(neutro@filter, type = "l")
  lines(not_neutro@filter, type = "l")
  legend("topleft", legend = c(
    paste0("count: ", neutro@cell.count),
    paste0("frequency: ", round(neutro@proportion, digits = 2)),
    paste0("total: ", not_cd3@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### Macrophages
  macro <- flowDensity(after[not_neutro@index], channels = c("PE-Cy7-A", "APC-A"), position = c(T, NA))
  not_macro <- flowDensity(after[not_neutro@index], channels = c("PE-Cy7-A", "APC-A"), position = c(F, NA))
  bmp(
    file = "Macrophages.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[not_neutro@index], c("PE-Cy7-A", "APC-A"), main = "Macrophages")
  lines(not_macro@filter, type = "l", col = 2)
  lines(macro@filter, type = "l")
  # lines(not_cd3@filter,type="l")
  legend("topleft", legend = c(
    paste0("count: ", macro@cell.count),
    paste0("frequency: ", round(macro@proportion, digits = 2)),
    paste0("total: ", not_neutro@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### CD206High
  cd206 <- flowDensity(after[macro@index], channels = c("PE-Cy7-A", "FSC-A"), position = c(T, NA))
  not_cd206 <- flowDensity(after[macro@index], channels = c("PE-Cy7-A", "FSC-A"), position = c(F, NA), upper = c(F, NA))
  bmp(
    file = "CD206.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[macro@index], c("PE-Cy7-A", "FSC-A"), main = "CD206")
  lines(cd206@filter, type = "l")
  lines(not_cd206@filter, type = "l", col = 2)

  legend("topleft", legend = c(
    paste0("count: ", cd206@cell.count),
    paste0("frequency: ", round(cd206@proportion, digits = 2)),
    paste0("total: ", macro@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### Eosinophils
  eosin <- flowDensity(after[not_macro@index], channels = c("APC-A", "SSC-A"), position = c(NA, T))
  not_eosin <- flowDensity(after[not_macro@index], channels = c("APC-A", "SSC-A"), position = c(NA, F))
  bmp(
    file = "Eosinophils.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[not_macro@index], c("APC-A", "SSC-A"), main = "Eosinophils")
  lines(not_eosin@filter, type = "l", col = 2)
  lines(eosin@filter, type = "l")
  # lines(not_neutro@filter,type="l")
  legend("topleft", legend = c(
    paste0("count: ", eosin@cell.count),
    paste0("frequency: ", round(eosin@proportion, digits = 2)),
    paste0("total: ", not_macro@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### NK cells
  nk <- flowDensity(after[not_eosin@index], channels = c("DAPI-A", "Side Pop-A"), position = c(F, T))
  not_nk <- flowDensity(after[not_eosin@index], channels = c("DAPI-A", "Side Pop-A"), position = c(NA, F))
  bmp(
    file = "NK cells.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[not_eosin@index], c("DAPI-A", "Side Pop-A"), main = "NK cells")
  lines(not_nk@filter, type = "l", col = 2)
  lines(nk@filter, type = "l")
  # lines(not_neutro@filter,type="l")
  legend("topleft", legend = c(
    paste0("count: ", nk@cell.count),
    paste0("frequency: ", round(nk@proportion, digits = 2)),
    paste0("total: ", not_eosin@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### Monocytes
  mono <- flowDensity(after[not_nk@index], channels = c("DAPI-A", "APC-A"), position = c(NA, T))
  not_mono <- flowDensity(after[not_nk@index], channels = c("DAPI-A", "APC-A"), position = c(NA, F))
  bmp(
    file = "Monocytes.bmp",
    width = 6, height = 4, units = "in", res = 100
  )
  plotDens(after[not_nk@index], c("DAPI-A", "APC-A"), main = "Monocytes")
  # lines(not_cd45@filter,type="l",col=2,lwd=4)
  lines(nk@filter, type = "l")
  # lines(not_neutro@filter,type="l")
  legend("topleft", legend = c(
    paste0("count: ", mono@cell.count),
    paste0("frequency: ", round(mono@proportion, digits = 2)),
    paste0("total: ", not_nk@cell.count)
  ), bty = "n", text.col = 1)
  dev.off()
  ### creating dataframe with info on percentages and counts
  Cell <- c(
    "Singlets", "CD45", "Live", "CD3+", "CD3-", "CD8 T cells", "CD4 T cells", "Tregs", "Neutrophils", "Not Neutrophils", "Marophages",
    "Not Macrophages", "CD206 high macrophages", "CD206 low macrophages", "Eosinophils", "NK cells", "Not NK", "Monocytes", "B and plasma cells"
  )
  Cell_count <- c(
    sngl@cell.count, cd45@cell.count, live@cell.count, cd3@cell.count, not_cd3@cell.count,
    cd8@cell.count, cd4@cell.count, tregs@cell.count, neutro@cell.count, not_neutro@cell.count,
    macro@cell.count, not_macro@cell.count, cd206@cell.count, not_cd206@cell.count,
    eosin@cell.count, nk@cell.count, not_nk@cell.count, mono@cell.count, not_mono@cell.count
  )
  Cell_parent_percentage <- c(
    sngl@proportion, cd45@proportion, live@proportion, cd3@proportion, not_cd3@proportion ,
    cd8@proportion, cd4@proportion, tregs@proportion, neutro@proportion, not_neutro@proportion,
    macro@proportion, not_macro@proportion, cd206@proportion, not_cd206@proportion,
    eosin@proportion, nk@proportion, not_nk@proportion, mono@proportion, not_mono@proportion
  )
  Cell_abs_percentage <- c(
    NA, NA, 1, (cd3@cell.count / live@cell.count) * 100, (not_cd3@cell.count / live@cell.count) * 100,
    (cd8@cell.count / live@cell.count) * 100, (cd4@cell.count / live@cell.count) * 100, (tregs@cell.count / live@cell.count) * 100, (neutro@cell.count / live@cell.count) * 100, (not_neutro@cell.count / live@cell.count) * 100,
    (macro@cell.count / live@cell.count) * 100, (not_macro@cell.count / live@cell.count) * 100, (cd206@cell.count / live@cell.count) * 100, (not_cd206@cell.count / live@cell.count) * 100,
    (eosin@cell.count / live@cell.count) * 100, (nk@cell.count / live@cell.count) * 100, (not_nk@cell.count / live@cell.count) * 100, (mono@cell.count / live@cell.count) * 100, (not_mono@cell.count / live@cell.count) * 100
  )
  Sample <- i
  df <- data.frame(Cell, Cell_count, Cell_parent_percentage, Cell_abs_percentage, Sample)
  write.csv(df, paste0(i, ".csv"), row.names = FALSE)
  setwd(paste0(cwd_0, "/PASC_pipeline"))
}
setwd(cwd_0)
