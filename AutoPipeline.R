library(flowCore)
library(flowDensity)
library(flowViz)
library(flowWorkspace)
fcs_test <- read.FCS("/Users/mgj1343/Library/CloudStorage/OneDrive-NorthwesternUniversity/SBratchikov/FlowDensity/01_data/PASC/02_BAL_flow_cytometry/PASC/20210416 PASC0020/20210416_LC002_001.fcs")
### singlets
sngl <- flowDensity(fcs_test,
  channels = c("FSC-A", "FSC-H"), position = c(F, F),
  percentile = c(.99999, .99999), use.percentile = c(T, T),
  ellip.gate = T, scale = .8
)
plotDens(fcs_test, c(1, 2), main = "Singlets")
lines(sngl@filter, type = "l")
legend("topleft", legend = c(
  paste0("count: ", sngl@cell.count),
  paste0("frequency: ", round(sngl@proportion, digits = 2)),
  paste0("total: ", nrow(sngl@flow.frame))
), bty = "n")

### cd45+
trans <- transformList(c("Qdot 655-A", "AmCyan-A"), lgcl)
after <- transform(sngl@flow.frame, trans)

cd45 <- flowDensity(after, channels = c("AmCyan-A", "Qdot 655-A"), position = c(T, NA), upper = c(F, NA))
# not_cd45 <- flowDensity(after, channels=c("AmCyan-A", "Qdot 655-A"), position=c(F,NA),upper=c(F,NA))
plotDens(after, c("AmCyan-A", "Qdot 655-A"), main = "CD45+")
# lines(not_cd45@filter,type="l",col=2,lwd=4)
lines(cd45@filter, type = "l")
legend("topleft", legend = c(
  paste0("count: ", cd45@cell.count),
  paste0("frequency: ", round(cd45@proportion, digits = 2)),
  paste0("total: ", sngl@cell.count)
), bty = "n", text.col = 1)
### logicle transformation
lgcl <- logicleTransform(w = 0.5, t = 262144, m = 4.5)
trans <- transformList(c("FITC-A", "FSC-A"), lgcl)
after <- transform(cd45@flow.frame, trans)
### Live singlet cells
live <- flowDensity(after, channels = c("FITC-A", "FSC-A"), position = c(F, NA), upper = c(NA, F))
plotDens(after, c("FITC-A", "FSC-A"), main = "Live")
legend("topleft", legend = c(
  paste0("count: ", live@cell.count),
  paste0("frequency: ", round(live@proportion, digits = 2)),
  paste0("total: ", cd45@cell.count)
), bty = "n", text.col = 1)
lines(live@filter, type = "l")

### cd3+
lgcl <- estimateLogicle(live@flow.frame, channels = c("PE-A", "PE-Cy7-A"), na.rm = TRUE)
# trans <- transformList(c("PE-A","PE-Cy7-A"),lgcl)
after <- transform(live@flow.frame, lgcl)

cd3 <- flowDensity(after, channels = c("PE-A", "PE-Cy7-A"), position = c(T, F))
not_cd3 <- flowDensity(after, channels = c("PE-A", "PE-Cy7-A"), position = c(F, NA))
plotDens(after, c("PE-A", "PE-Cy7-A"), main = "CD3 T cells")
lines(not_cd3@filter, type = "l", col = 2)
lines(cd3@filter, type = "l")
# lines(not_cd3@filter,type="l")
legend("topleft", legend = c(
  paste0("count: ", cd3@cell.count),
  paste0("frequency: ", round(cd3@proportion, digits = 2)),
  paste0("total: ", live@cell.count)
), bty = "n", text.col = 1)
### CD8 CD4
trans <- transformList(c("DAPI-A", "APC-A"), lgcl)
after <- transform(cd3@flow.frame, trans)
cd4 <- flowDensity(after, c("DAPI-A", "APC-A"), position = c(T, F))
cd8 <- flowDensity(after, c("DAPI-A", "APC-A"), position = c(F, T))
# not_cd45 <- flowDensity(after, channels=c("AmCyan-A", "Qdot 655-A"), position=c(F,NA),upper=c(F,NA))
plotDens(after, c("DAPI-A", "APC-A"), main = "CD8/CD4")
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

### tregs
trans <- transformList(c("Side Pop-A", "PE-Texas Red-A"), lgcl)
after <- transform(cd4@flow.frame, trans)

tregs <- flowDensity(after, channels = c("Side Pop-A", "PE-Texas Red-A"), position = c(T, F))
# not_cd45 <- flowDensity(after, channels=c("AmCyan-A", "Qdot 655-A"), position=c(F,NA),upper=c(F,NA))
plotDens(after, c("Side Pop-A", "PE-Texas Red-A"), main = "Tregs")
# lines(not_cd45@filter,type="l",col=2,lwd=4)
lines(tregs@filter, type = "l")
legend("topleft", legend = c(
  paste0("count: ", tregs@cell.count),
  paste0("frequency: ", round(tregs@proportion, digits = 2)),
  paste0("total: ", cd4@cell.count)
), bty = "n", text.col = 1)

### Neutrophils
trans <- transformList(c("Pacific Blue-A"), lgcl)
after <- transform(not_cd3@flow.frame, trans)

neutro <- flowDensity(after, channels = c("Qdot 655-A", "Pacific Blue-A"), position = c(T, NA))
not_neutro <- flowDensity(after, channels = c("Qdot 655-A", "Pacific Blue-A"), position = c(F, NA))
plotDens(after, c("Qdot 655-A", "Pacific Blue-A"), main = "Neutrophils")
# lines(not_cd45@filter,type="l",col=2,lwd=4)
lines(neutro@filter, type = "l")
lines(not_neutro@filter, type = "l")
legend("topleft", legend = c(
  paste0("count: ", neutro@cell.count),
  paste0("frequency: ", round(neutro@proportion, digits = 2)),
  paste0("total: ", not_cd3@cell.count)
), bty = "n", text.col = 1)
### Macrophages
trans <- transformList(c("APC-A"), lgcl)
after <- transform(not_neutro@flow.frame, trans)

macro <- flowDensity(after, channels = c("PE-Cy7-A", "APC-A"), position = c(T, NA))
not_macro <- flowDensity(after, channels = c("PE-Cy7-A", "APC-A"), position = c(F, NA))
plotDens(after, c("PE-Cy7-A", "APC-A"), main = "Macrophages")
lines(not_macro@filter, type = "l", col = 2)
lines(macro@filter, type = "l")
# lines(not_cd3@filter,type="l")
legend("topleft", legend = c(
  paste0("count: ", macro@cell.count),
  paste0("frequency: ", round(macro@proportion, digits = 2)),
  paste0("total: ", not_neutro@cell.count)
), bty = "n", text.col = 1)
### CD206High
cd206 <- flowDensity(macro, channels = c("PE-Cy7-A", "FSC-A"), position = c(T, NA))
not_cd206 <- flowDensity(macro, channels = c("PE-Cy7-A", "FSC-A"), position = c(F, NA))
plotDens(macro, c("PE-Cy7-A", "FSC-A"), main = "CD206")
lines(cd206@filter, type = "l")
lines(not_cd206@filter, type = "l", col = 2)

legend("topleft", legend = c(
  paste0("count: ", cd206@cell.count),
  paste0("frequency: ", round(cd206@proportion, digits = 2)),
  paste0("total: ", macro@cell.count)
), bty = "n", text.col = 1)
### Eosinophils

trans <- transformList(c("SSC-A"), lgcl)
after <- transform(not_macro@flow.frame, trans)

eosin <- flowDensity(after, channels = c("APC-A", "SSC-A"), position = c(NA, T))
not_eosin <- flowDensity(after, channels = c("APC-A", "SSC-A"), position = c(NA, F))
plotDens(after, c("APC-A", "SSC-A"), main = "Eosinophils")
lines(not_eosin@filter, type = "l", col = 2)
lines(eosin@filter, type = "l")
# lines(not_neutro@filter,type="l")
legend("topleft", legend = c(
  paste0("count: ", eosin@cell.count),
  paste0("frequency: ", round(eosin@proportion, digits = 2)),
  paste0("total: ", not_macro@cell.count)
), bty = "n", text.col = 1)

### NK cells

trans <- transformList(c("DAPI-A", "Side Pop-A"), lgcl)
after <- transform(not_eosin@flow.frame, trans)

nk <- flowDensity(after, channels = c("DAPI-A", "Side Pop-A"), position = c(F, T))
not_nk <- flowDensity(after, channels = c("DAPI-A", "Side Pop-A"), position = c(NA, F))
plotDens(after, c("DAPI-A", "Side Pop-A"), main = "NK cells")
lines(not_nk@filter, type = "l", col = 2)
lines(nk@filter, type = "l")
# lines(not_neutro@filter,type="l")
legend("topleft", legend = c(
  paste0("count: ", nk@cell.count),
  paste0("frequency: ", round(nk@proportion, digits = 2)),
  paste0("total: ", not_eosin@cell.count)
), bty = "n", text.col = 1)

### Monocytes
mono <- flowDensity(not_nk@flow.frame, channels = c("DAPI-A", "APC-A"), position = c(NA, T))
not_mono <- flowDensity(after, channels = c("DAPI-A", "APC-A"), position = c(NA, F))
plotDens(not_nk, c("DAPI-A", "APC-A"), main = "Monocytes")
# lines(not_cd45@filter,type="l",col=2,lwd=4)
lines(mono@filter, type = "l")
# lines(not_neutro@filter,type="l")
legend("topleft", legend = c(
  paste0("count: ", mono@cell.count),
  paste0("frequency: ", round(mono@proportion, digits = 2)),
  paste0("total: ", not_nk@cell.count)
), bty = "n", text.col = 1)
### creating dataframe with info on percentages and counts
Cell <- c(
  "Singlets", "CD45", "Live", "CD3+", "CD3-", "CD8 T cells", "CD4 T cells", "Tregs", "Neutrophils", "Not Neutrophils", "Marophages",
  "Not Macrophages", "CD206 high macrophages", "CD206 low macrophages", "Eosinophils", "NK cells", "Not NK", "Monocytes"
)
Cell_count <- c(
  sngl@cell.count, cd45@cell.count, live@cell.count, cd3@cell.count, not_cd3@cell.count,
  cd8@cell.count, cd4@cell.count, tregs@cell.count, neutro@cell.count, not_neutro@cell.count,
  macro@cell.count, not_macro@cell.count, cd206@cell.count, not_cd206@cell.count,
  eosin@cell.count, nk@cell.count, not_nk@cell.count, mono@cell.count
)
Cell_parent_percentage <- c(
  sngl@proportion, cd45@proportion, live@proportion, cd3@proportion, not_cd3@proportion,
  cd8@proportion, cd4@proportion, tregs@proportion, neutro@proportion, not_neutro@proportion,
  macro@proportion, not_macro@proportion, cd206@proportion, not_cd206@proportion,
  eosin@proportion, nk@proportion, not_nk@proportion, mono@proportion
)
Sample <- c("20210120_1463-BAL-00_001.fcs")

df <- data.frame(Cell, Cell_count, Cell_parent_percentage, Sample)
