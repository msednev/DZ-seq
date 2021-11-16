# this script collects count data into a single dataset
# with the following columns:
#     core      - sequence of deoxyribozyme core
#     substrate - residual sequence of RNA substrate
#     count     - frequency of core+substrate combination
#     pool_mod  - target modification for in-vitro selection
#     sub_mod   - modification in the RNA substrate
#     rnd       - in-vitro selection round number
#     lib_name  - library name as indicated in the SI

library(data.table)
library(stringi)



ak7u <- fread('./output/counts/AK7U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
ak7u[, pool_mod := 'C']
ak7u[, sub_mod := 'm5C']
ak7u[, rnd := 7L]
ak7u[, lib_name := "L3"]

ak7m <- fread('./output/counts/AK7M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
ak7m[, pool_mod := 'C']
ak7m[, sub_mod := 'C']
ak7m[, rnd := 7L]
ak7m[, lib_name := "L1"]

ak18u <- fread('./output/counts/AK18U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
ak18u[, pool_mod := 'C']
ak18u[, sub_mod := 'm5C']
ak18u[, rnd := 18L]
ak18u[, lib_name := "L4"]

ak18m <- fread('./output/counts/AK18M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
ak18m[, pool_mod := 'C']
ak18m[, sub_mod := 'C']
ak18m[, rnd := 18L]
ak18m[, lib_name := "L2"]

al7u <- fread('./output/counts/AL7U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
al7u[, pool_mod := 'm3C']
al7u[, sub_mod := 'C']
al7u[, rnd := 7L]
al7u[, lib_name := "L5"]

al7m <- fread('./output/counts/AL7M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
al7m[, pool_mod := 'm3C']
al7m[, sub_mod := 'm3C']
al7m[, rnd := 7L]
al7m[, lib_name := "L7"]

al18u <- fread('./output/counts/AL18U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
al18u[, pool_mod := 'm3C']
al18u[, sub_mod := 'C']
al18u[, rnd := 18L]
al18u[, lib_name := "L6"]

al18m <- fread('./output/counts/AL18M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
al18m[, pool_mod := 'm3C']
al18m[, sub_mod := 'm3C']
al18m[, rnd := 18L]
al18m[, lib_name := "L8"]

am7u <- fread('./output/counts/AM7U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
am7u[, pool_mod := 'm4C']
am7u[, sub_mod := 'C']
am7u[, rnd := 7L]
am7u[, lib_name := "L9"]

am7m <- fread('./output/counts/AM7M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
am7m[, pool_mod := 'm4C']
am7m[, sub_mod := 'm4C']
am7m[, rnd := 7L]
am7m[, lib_name := "L11"]

am18u <- fread('./output/counts/AM18U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
am18u[, pool_mod := 'm4C']
am18u[, sub_mod := 'C']
am18u[, rnd := 18L]
am18u[, lib_name := "L10"]

am18m <- fread('./output/counts/AM18M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
am18m[, pool_mod := 'm4C']
am18m[, sub_mod := 'm4C']
am18m[, rnd := 18L]
am18m[, lib_name := "L12"]

an7u <- fread('./output/counts/AN7U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
an7u[, pool_mod := 'm5C']
an7u[, sub_mod := 'C']
an7u[, rnd := 7L]
an7u[, lib_name := "L13"]

an7m <- fread('./output/counts/AN7M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
an7m[, pool_mod := 'm5C']
an7m[, sub_mod := 'm5C']
an7m[, rnd := 7L]
an7m[, lib_name := "L15"]

an18u <- fread('./output/counts/AN18U.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
an18u[, pool_mod := 'm5C']
an18u[, sub_mod := 'C']
an18u[, rnd := 18L]
an18u[, lib_name := "L14"]

an18m <- fread('./output/counts/AN18M.tsv', col.names = c('core', 'substrate', 'count'), na.strings = "")
an18m[, pool_mod := 'm5C']
an18m[, sub_mod := 'm5C']
an18m[, rnd := 18L]
an18m[, lib_name := "L16"]

DT <- rbindlist(list(ak7u, ak7m, ak18u, ak18m,
                     al7u, al7m, al18u, al18m,
                     am7u, am7m, am18u, am18m,
                     an7u, an7m, an18u, an18m))

rm(ak7u, ak7m, ak18u, ak18m,
   al7u, al7m, al18u, al18m,
   am7u, am7m, am18u, am18m,
   an7u, an7m, an18u, an18m)

# remove entries with too short or too long sequences for DNAzyme core or RNA substrate
DT <- DT[!(stri_length(core) < 12 | stri_length(core) > 25 | stri_length(substrate) > 11 | stri_length(substrate) < 4)]

fwrite(DT, './output/datasets/counts.csv.gz')