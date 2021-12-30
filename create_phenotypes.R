#!/bin/R

library(data.table)
set.seed(5)

#ped <- fread("integrated_call_samples_v3.20130502.ALL.panel", data.table = F, header=T, stringsAsFactors = F, fill = T)
#ped[-1,]

ped <- fread("1000G/integrated_call_samples_v3.20130502.SAS.panel", data.table = F, header=T, stringsAsFactors = F, fill = T)
ped[-1,]

n <- nrow(ped)
ped$age <- floor(runif(n, min = 30, max = 99))
ped$bmi <- ped$t2d <- 1
colnames(ped)

t2d.prev <- data.frame(
  anc = c("EUR", "AMR", "AFR", "EAS", "SAS"),
  prev = c(7.1, 33, 12.6, 8.95, 17),
  case_bmi = c(20, 35, 30, 15, 10),
  control_bmi = c(25, 30, 30, 17, 8),
  stringsAsFactors = F
)
#t2d.prev

print("calculating new ped")
new.ped <- list()
for (anc in unique(ped$super_pop)){
  cur.ped <- ped[ped$super_pop == anc,]
  nn <- nrow(cur.ped)
  ncase <- floor(nn*t2d.prev[t2d.prev$anc == anc, "prev"]/100)
  cur.ped[sample(1:nrow(cur.ped), ncase), "t2d"] <- 2
  cur.ped[cur.ped$t2d == 2, "bmi"] <- rnorm(nrow(cur.ped[cur.ped$t2d == 2,]),
                                            mean = t2d.prev[t2d.prev$anc == anc, "case_bmi"],
                                            sd = 2)
  cur.ped[cur.ped$t2d == 1, "bmi"] <- rnorm(nrow(cur.ped[cur.ped$t2d == 1,]),
                                            mean = t2d.prev[t2d.prev$anc == anc, "control_bmi"],
                                            sd = 2)
  new.ped[[anc]] <- cur.ped
}

head(new.ped)
ped <- do.call(rbind, new.ped)

table(ped$gender)

ped$gender <- ifelse(ped$gender == "male", "M", "F")

head(ped)


#fwrite(ped,
#       file = paste0('1kg_phenotype_mock_ALL_', as.Date(Sys.time(), "%d/%m/%Y"), ".tsv"),
#       sep = "\t", row.names = F)

fwrite(ped,
       file = paste0('1000G/1kg_phenotype_mock_SAS_', as.Date(Sys.time(), "%d/%m/%Y"), ".tsv"),
       sep = "\t", row.names = F)
