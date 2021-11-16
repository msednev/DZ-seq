# this script selects entries with the highest FC among
# all cleavage sites and with the number of counts greater
# than CUTOFF and combines resulting data for all RNA 
# substrates into single dataset

library(data.table)

CUTOFF = 99



df_GC <- fread('./output/datasets/cleavage_GC.csv.gz')

dff_GC <- df_GC[order(-FC), head(.SD, 1), by = .(core, pool_mod, sub_mod, rna, rnd)][(cleaved + uncleaved) > CUTOFF]

rm(df_GC)
gc()

df_AC <- fread('./output/datasets/cleavage_AC.csv.gz')

dff_AC <- df_AC[order(-FC), head(.SD, 1), by = .(core, pool_mod, sub_mod, rna, rnd)][(cleaved + uncleaved) > CUTOFF]

rm(df_AC)
gc()

df_CC <- fread('./output/datasets/cleavage_CC.csv.gz')

dff_CC <- df_CC[order(-FC), head(.SD, 1), by = .(core, pool_mod, sub_mod, rna, rnd)][(cleaved + uncleaved) > CUTOFF]

rm(df_CC)
gc()

df_TC <- fread('./output/datasets/cleavage_TC.csv.gz')

dff_TC <- df_TC[order(-FC), head(.SD, 1), by = .(core, pool_mod, sub_mod, rna, rnd)][(cleaved + uncleaved) > CUTOFF]

rm(df_TC)
gc()

dff <- rbindlist(list(dff_GC, dff_AC, dff_CC, dff_TC))

fwrite(dff, './output/datasets/cleavage_filtered.csv.gz')
