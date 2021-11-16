# this script assigns each residual RNA sequence in the count dataset
# to one of four possible RNA substrates (indicated as GC, AC, CC and TC),
# identifies its cleavage status and calculates FC for every possible cleavage
# position

# data is processed separately for each RNA substrate to reduce RAM usage

# cleavage site       1 2 3 4 5 6 7 8 
# RNA substrate: G A A N C G T A A C T

# if RNA is cleaved at position 1, the resulting residue is impossible 
# to assign to any of the RNA substrates -> average FC is calculated

library(data.table)


DT <- fread('./output/datasets/counts.csv.gz')





############## GC RNA ############## 
if (!file.exists('output/datasets/cleavage_GC.csv.gz')) {
  dt_GC <-
    DT[, .(
      uncleaved_1 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC', 'GAAGCGTAA', 'GAAGCGTA', 'GAAGCGT', 'GAAGCG', 'GAAGC', 'GAAG',
                                           'GAAGTGTAACT', 'GAAGTGTAAC', 'GAAGTGTAA', 'GAAGTGTA', 'GAAGTGT', 'GAAGTG', 'GAAGT',
                                           'GAAACGTAACT', 'GAAACGTAAC', 'GAAACGTAA', 'GAAACGTA', 'GAAACGT', 'GAAACG', 'GAAAC',
                                           'GAAATGTAACT', 'GAAATGTAAC', 'GAAATGTAA', 'GAAATGTA', 'GAAATGT', 'GAAATG', 'GAAAT',
                                           'GAACCGTAACT', 'GAACCGTAAC', 'GAACCGTAA', 'GAACCGTA', 'GAACCGT', 'GAACCG', 'GAACC', 'GAAC',
                                           'GAACTGTAACT', 'GAACTGTAAC', 'GAACTGTAA', 'GAACTGTA', 'GAACTGT', 'GAACTG', 'GAACT',
                                           'GAATCGTAACT', 'GAATCGTAAC', 'GAATCGTAA', 'GAATCGTA', 'GAATCGT', 'GAATCG', 'GAATC', 'GAAT',
                                           'GAATTGTAACT', 'GAATTGTAAC', 'GAATTGTAA', 'GAATTGTA', 'GAATTGT', 'GAATTG', 'GAATT') * count),
      cleaved_1   = sum(substrate %chin% c('GAAAAAAAAAA', 'GAAAAAAAAA', 'GAAAAAAAA', 'GAAAAAAA', 'GAAAAAA', 'GAAAAA', 'GAAAA', 'GAAA') * count),
      uncleaved_2 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC', 'GAAGCGTAA', 'GAAGCGTA', 'GAAGCGT', 'GAAGCG', 'GAAGC',
                                           'GAAGTGTAACT', 'GAAGTGTAAC', 'GAAGTGTAA', 'GAAGTGTA', 'GAAGTGT', 'GAAGTG', 'GAAGT') * count),
      cleaved_2   = sum(substrate %chin% c('GAAGAAAAAAA', 'GAAGAAAAAA', 'GAAGAAAAA', 'GAAGAAAA', 'GAAGAAA', 'GAAGAA', 'GAAGA') * count),
      uncleaved_3 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC', 'GAAGCGTAA', 'GAAGCGTA', 'GAAGCGT', 'GAAGCG', 
                                           'GAAGTGTAACT', 'GAAGTGTAAC', 'GAAGTGTAA', 'GAAGTGTA', 'GAAGTGT', 'GAAGTG') * count),
      cleaved_3   = sum(substrate %chin% c('GAAGCAAAAAA', 'GAAGCAAAAA', 'GAAGCAAAA', 'GAAGCAAA', 'GAAGCAA', 'GAAGCA',
                                           'GAAGTAAAAAA', 'GAAGTAAAAA', 'GAAGTAAAA', 'GAAGTAAA', 'GAAGTAA', 'GAAGTA') * count),
      uncleaved_4 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC', 'GAAGCGTAA', 'GAAGCGTA', 'GAAGCGT',
                                           'GAAGTGTAACT', 'GAAGTGTAAC', 'GAAGTGTAA', 'GAAGTGTA', 'GAAGTGT') * count),
      cleaved_4   = sum(substrate %chin% c('GAAGCGAAAAA', 'GAAGCGAAAA', 'GAAGCGAAA', 'GAAGCGAA', 'GAAGCGA',
                                           'GAAGTGAAAAA', 'GAAGTGAAAA', 'GAAGTGAAA', 'GAAGTGAA', 'GAAGTGA') * count),
      uncleaved_5 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC',
                                           'GAAGTGTAACT', 'GAAGTGTAAC') * count),
      cleaved_5   = sum(substrate %chin% c('GAAGCGTAAAA', 'GAAGCGTAAA',
                                           'GAAGTGTAAAA', 'GAAGTGTAAA') * count),
      uncleaved_8 = sum(substrate %chin% c('GAAGCGTAACT',
                                           'GAAGTGTAACT') * count),
      cleaved_8   = sum(substrate %chin% c('GAAGCGTAACA',
                                           'GAAGTGTAACA') * count)
    ),
    by = .(core, pool_mod, sub_mod, rnd, lib_name)]
  
  dt_GC <-
    melt(
      dt_GC,
      id.vars = c('core', 'pool_mod', 'sub_mod', 'rnd', 'lib_name'),
      measure.vars = patterns(cleaved = '^cleaved_', uncleaved = '^uncleaved_'),
      variable.name = 'pos'
    )
  
  levels(dt_GC$pos) <- c(1, 2, 3, 4, 5, 8)
  
  dt_GC[, FC := cleaved / (uncleaved + cleaved)]
  dt_GC[, rna := 'GC']
  
  fwrite(dt_GC, './output/datasets/cleavage_GC.csv.gz')
  rm(dt_GC)
  gc()
  }

############## AC RNA ##############
if (!file.exists('output/datasets/cleavage_AC.csv.gz')) {
  dt_AC <-
    DT[, .(
      uncleaved_1 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC', 'GAAGCGTAA', 'GAAGCGTA', 'GAAGCGT', 'GAAGCG', 'GAAGC', 'GAAG',
                                           'GAAGTGTAACT', 'GAAGTGTAAC', 'GAAGTGTAA', 'GAAGTGTA', 'GAAGTGT', 'GAAGTG', 'GAAGT',
                                           'GAAACGTAACT', 'GAAACGTAAC', 'GAAACGTAA', 'GAAACGTA', 'GAAACGT', 'GAAACG', 'GAAAC',
                                           'GAAATGTAACT', 'GAAATGTAAC', 'GAAATGTAA', 'GAAATGTA', 'GAAATGT', 'GAAATG', 'GAAAT',
                                           'GAACCGTAACT', 'GAACCGTAAC', 'GAACCGTAA', 'GAACCGTA', 'GAACCGT', 'GAACCG', 'GAACC', 'GAAC',
                                           'GAACTGTAACT', 'GAACTGTAAC', 'GAACTGTAA', 'GAACTGTA', 'GAACTGT', 'GAACTG', 'GAACT',
                                           'GAATCGTAACT', 'GAATCGTAAC', 'GAATCGTAA', 'GAATCGTA', 'GAATCGT', 'GAATCG', 'GAATC', 'GAAT',
                                           'GAATTGTAACT', 'GAATTGTAAC', 'GAATTGTAA', 'GAATTGTA', 'GAATTGT', 'GAATTG', 'GAATT') * count),
      cleaved_1   = sum(substrate %chin% c('GAAAAAAAAAA', 'GAAAAAAAAA', 'GAAAAAAAA', 'GAAAAAAA', 'GAAAAAA', 'GAAAAA', 'GAAAA') * count),
      uncleaved_2 = sum(substrate %chin% c('GAAACGTAACT', 'GAAACGTAAC', 'GAAACGTAA', 'GAAACGTA', 'GAAACGT', 'GAAACG', 'GAAAC',
                                           'GAAATGTAACT', 'GAAATGTAAC', 'GAAATGTAA', 'GAAATGTA', 'GAAATGT', 'GAAATG', 'GAAAT') * count),
      cleaved_2   = sum(substrate %chin% c('GAAAAAAAAAA', 'GAAAAAAAAA', 'GAAAAAAAA', 'GAAAAAAA', 'GAAAAAA', 'GAAAAA', 'GAAAA') * count),
      uncleaved_3 = sum(substrate %chin% c('GAAACGTAACT', 'GAAACGTAAC', 'GAAACGTAA', 'GAAACGTA', 'GAAACGT', 'GAAACG', 
                                           'GAAATGTAACT', 'GAAATGTAAC', 'GAAATGTAA', 'GAAATGTA', 'GAAATGT', 'GAAATG') * count),
      cleaved_3   = sum(substrate %chin% c('GAAACAAAAAA', 'GAAACAAAAA', 'GAAACAAAA', 'GAAACAAA', 'GAAACAA', 'GAAACA',
                                           'GAAATAAAAAA', 'GAAATAAAAA', 'GAAATAAAA', 'GAAATAAA', 'GAAATAA', 'GAAATA') * count),
      uncleaved_4 = sum(substrate %chin% c('GAAACGTAACT', 'GAAACGTAAC', 'GAAACGTAA', 'GAAACGTA', 'GAAACGT',
                                           'GAAATGTAACT', 'GAAATGTAAC', 'GAAATGTAA', 'GAAATGTA', 'GAAATGT') * count),
      cleaved_4   = sum(substrate %chin% c('GAAACGAAAAA', 'GAAACGAAAA', 'GAAACGAAA', 'GAAACGAA', 'GAAACGA',
                                           'GAAATGAAAAA', 'GAAATGAAAA', 'GAAATGAAA', 'GAAATGAA', 'GAAATGA') * count),
      uncleaved_5 = sum(substrate %chin% c('GAAACGTAACT', 'GAAACGTAAC',
                                           'GAAATGTAACT', 'GAAATGTAAC') * count),
      cleaved_5   = sum(substrate %chin% c('GAAACGTAAAA', 'GAAACGTAAA',
                                           'GAAATGTAAAA', 'GAAATGTAAA') * count),
      uncleaved_8 = sum(substrate %chin% c('GAAACGTAACT',
                                           'GAAATGTAACT') * count),
      cleaved_8   = sum(substrate %chin% c('GAAGCGTAACA',
                                           'GAAGTGTAACA') * count)
    ),
    by = .(core, pool_mod, sub_mod, rnd, lib_name)]
  
  dt_AC <-
    melt(
      dt_AC,
      id.vars = c('core', 'pool_mod', 'sub_mod', 'rnd', 'lib_name'),
      measure.vars = patterns(cleaved = '^cleaved_', uncleaved = '^uncleaved_'),
      variable.name = 'pos'
    )
  
  levels(dt_AC$pos) <- c(1, 2, 3, 4, 5, 8)
  
  dt_AC[, FC := cleaved / (uncleaved + cleaved)]
  dt_AC[, rna := 'AC']
  
  fwrite(dt_AC, './output/datasets/cleavage_AC.csv.gz')
  rm(dt_AC)
  gc()
  }

############## CC RNA ##############
if (!file.exists('output/datasets/cleavage_CC.csv.gz')) {
  dt_CC <-
    DT[, .(
      uncleaved_1 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC', 'GAAGCGTAA', 'GAAGCGTA', 'GAAGCGT', 'GAAGCG', 'GAAGC', 'GAAG',
                                           'GAAGTGTAACT', 'GAAGTGTAAC', 'GAAGTGTAA', 'GAAGTGTA', 'GAAGTGT', 'GAAGTG', 'GAAGT',
                                           'GAAACGTAACT', 'GAAACGTAAC', 'GAAACGTAA', 'GAAACGTA', 'GAAACGT', 'GAAACG', 'GAAAC',
                                           'GAAATGTAACT', 'GAAATGTAAC', 'GAAATGTAA', 'GAAATGTA', 'GAAATGT', 'GAAATG', 'GAAAT',
                                           'GAACCGTAACT', 'GAACCGTAAC', 'GAACCGTAA', 'GAACCGTA', 'GAACCGT', 'GAACCG', 'GAACC', 'GAAC',
                                           'GAACTGTAACT', 'GAACTGTAAC', 'GAACTGTAA', 'GAACTGTA', 'GAACTGT', 'GAACTG', 'GAACT',
                                           'GAATCGTAACT', 'GAATCGTAAC', 'GAATCGTAA', 'GAATCGTA', 'GAATCGT', 'GAATCG', 'GAATC', 'GAAT',
                                           'GAATTGTAACT', 'GAATTGTAAC', 'GAATTGTAA', 'GAATTGTA', 'GAATTGT', 'GAATTG', 'GAATT') * count),
      cleaved_1   = sum(substrate %chin% c('GAAAAAAAAAA', 'GAAAAAAAAA', 'GAAAAAAAA', 'GAAAAAAA', 'GAAAAAA', 'GAAAAA', 'GAAAA', 'GAAA') * count),
      uncleaved_2 = sum(substrate %chin% c('GAACCGTAACT', 'GAACCGTAAC', 'GAACCGTAA', 'GAACCGTA', 'GAACCGT', 'GAACCG', 'GAACC',
                                           'GAACTGTAACT', 'GAACTGTAAC', 'GAACTGTAA', 'GAACTGTA', 'GAACTGT', 'GAACTG', 'GAACT') * count),
      cleaved_2   = sum(substrate %chin% c('GAACAAAAAAA', 'GAACAAAAAA', 'GAACAAAAA', 'GAACAAAA', 'GAACAAA', 'GAACAA', 'GAACA') * count),
      uncleaved_3 = sum(substrate %chin% c('GAACCGTAACT', 'GAACCGTAAC', 'GAACCGTAA', 'GAACCGTA', 'GAACCGT', 'GAACCG', 
                                           'GAACTGTAACT', 'GAACTGTAAC', 'GAACTGTAA', 'GAACTGTA', 'GAACTGT', 'GAACTG') * count),
      cleaved_3   = sum(substrate %chin% c('GAACCAAAAAA', 'GAACCAAAAA', 'GAACCAAAA', 'GAACCAAA', 'GAACCAA', 'GAACCA',
                                           'GAACTAAAAAA', 'GAACTAAAAA', 'GAACTAAAA', 'GAACTAAA', 'GAACTAA', 'GAACTA') * count),
      uncleaved_4 = sum(substrate %chin% c('GAACCGTAACT', 'GAACCGTAAC', 'GAACCGTAA', 'GAACCGTA', 'GAACCGT',
                                           'GAACTGTAACT', 'GAACTGTAAC', 'GAACTGTAA', 'GAACTGTA', 'GAACTGT') * count),
      cleaved_4   = sum(substrate %chin% c('GAACCGAAAAA', 'GAACCGAAAA', 'GAACCGAAA', 'GAACCGAA', 'GAACCGA',
                                           'GAACTGAAAAA', 'GAACTGAAAA', 'GAACTGAAA', 'GAACTGAA', 'GAACTGA') * count),
      uncleaved_5 = sum(substrate %chin% c('GAACCGTAACT', 'GAACCGTAAC',
                                           'GAACTGTAACT', 'GAACTGTAAC') * count),
      cleaved_5   = sum(substrate %chin% c('GAACCGTAAAA', 'GAACCGTAAA',
                                           'GAACTGTAAAA', 'GAACTGTAAA') * count),
      uncleaved_8 = sum(substrate %chin% c('GAACCGTAACT',
                                           'GAACTGTAACT') * count),
      cleaved_8   = sum(substrate %chin% c('GAACCGTAACA',
                                           'GAACTGTAACA') * count)
    ),
    by = .(core, pool_mod, sub_mod, rnd, lib_name)]
  
  dt_CC <-
    melt(
      dt_CC,
      id.vars = c('core', 'pool_mod', 'sub_mod', 'rnd', 'lib_name'),
      measure.vars = patterns(cleaved = '^cleaved_', uncleaved = '^uncleaved_'),
      variable.name = 'pos'
    )
  
  levels(dt_CC$pos) <- c(1, 2, 3, 4, 5, 8)
  
  dt_CC[, FC := cleaved / (uncleaved + cleaved)]
  dt_CC[, rna := 'CC']
  
  fwrite(dt_CC, './output/datasets/cleavage_CC.csv.gz')
  rm(dt_CC)
  gc()
  }


############## TC RNA ##############
if (!file.exists('output/datasets/cleavage_TC.csv.gz')) {
  dt_TC <-
    DT[, .(
      uncleaved_1 = sum(substrate %chin% c('GAAGCGTAACT', 'GAAGCGTAAC', 'GAAGCGTAA', 'GAAGCGTA', 'GAAGCGT', 'GAAGCG', 'GAAGC', 'GAAG',
                                           'GAAGTGTAACT', 'GAAGTGTAAC', 'GAAGTGTAA', 'GAAGTGTA', 'GAAGTGT', 'GAAGTG', 'GAAGT',
                                           'GAAACGTAACT', 'GAAACGTAAC', 'GAAACGTAA', 'GAAACGTA', 'GAAACGT', 'GAAACG', 'GAAAC',
                                           'GAAATGTAACT', 'GAAATGTAAC', 'GAAATGTAA', 'GAAATGTA', 'GAAATGT', 'GAAATG', 'GAAAT',
                                           'GAACCGTAACT', 'GAACCGTAAC', 'GAACCGTAA', 'GAACCGTA', 'GAACCGT', 'GAACCG', 'GAACC', 'GAAC',
                                           'GAACTGTAACT', 'GAACTGTAAC', 'GAACTGTAA', 'GAACTGTA', 'GAACTGT', 'GAACTG', 'GAACT',
                                           'GAATCGTAACT', 'GAATCGTAAC', 'GAATCGTAA', 'GAATCGTA', 'GAATCGT', 'GAATCG', 'GAATC', 'GAAT',
                                           'GAATTGTAACT', 'GAATTGTAAC', 'GAATTGTAA', 'GAATTGTA', 'GAATTGT', 'GAATTG', 'GAATT') * count),
      cleaved_1   = sum(substrate %chin% c('GAAAAAAAAAA', 'GAAAAAAAAA', 'GAAAAAAAA', 'GAAAAAAA', 'GAAAAAA', 'GAAAAA', 'GAAAA', 'GAAA') * count),
      uncleaved_2 = sum(substrate %chin% c('GAATCGTAACT', 'GAATCGTAAC', 'GAATCGTAA', 'GAATCGTA', 'GAATCGT', 'GAATCG', 'GAATC',
                                           'GAATTGTAACT', 'GAATTGTAAC', 'GAATTGTAA', 'GAATTGTA', 'GAATTGT', 'GAATTG', 'GAATT') * count),
      cleaved_2   = sum(substrate %chin% c('GAATAAAAAAA', 'GAATAAAAAA', 'GAATAAAAA', 'GAATAAAA', 'GAATAAA', 'GAATAA', 'GAATA') * count),
      uncleaved_3 = sum(substrate %chin% c('GAATCGTAACT', 'GAATCGTAAC', 'GAATCGTAA', 'GAATCGTA', 'GAATCGT', 'GAATCG', 
                                           'GAATTGTAACT', 'GAATTGTAAC', 'GAATTGTAA', 'GAATTGTA', 'GAATTGT', 'GAATTG') * count),
      cleaved_3   = sum(substrate %chin% c('GAATCAAAAAA', 'GAATCAAAAA', 'GAATCAAAA', 'GAATCAAA', 'GAATCAA', 'GAATCA',
                                           'GAATTAAAAAA', 'GAATTAAAAA', 'GAATTAAAA', 'GAATTAAA', 'GAATTAA', 'GAATTA') * count),
      uncleaved_4 = sum(substrate %chin% c('GAATCGTAACT', 'GAATCGTAAC', 'GAATCGTAA', 'GAATCGTA', 'GAATCGT',
                                           'GAATTGTAACT', 'GAATTGTAAC', 'GAATTGTAA', 'GAATTGTA', 'GAATTGT') * count),
      cleaved_4   = sum(substrate %chin% c('GAATCGAAAAA', 'GAATCGAAAA', 'GAATCGAAA', 'GAATCGAA', 'GAATCGA',
                                           'GAATTGAAAAA', 'GAATTGAAAA', 'GAATTGAAA', 'GAATTGAA', 'GAATTGA') * count),
      uncleaved_5 = sum(substrate %chin% c('GAATCGTAACT', 'GAATCGTAAC',
                                           'GAATTGTAACT', 'GAATTGTAAC') * count),
      cleaved_5   = sum(substrate %chin% c('GAATCGTAAAA', 'GAATCGTAAA',
                                           'GAATTGTAAAA', 'GAATTGTAAA') * count),
      uncleaved_8 = sum(substrate %chin% c('GAATCGTAACT',
                                           'GAATTGTAACT') * count),
      cleaved_8   = sum(substrate %chin% c('GAATCGTAACA',
                                           'GAATTGTAACA') * count)
    ),
    by = .(core, pool_mod, sub_mod, rnd, lib_name)]
  
  dt_TC <-
    melt(
      dt_TC,
      id.vars = c('core', 'pool_mod', 'sub_mod', 'rnd', 'lib_name'),
      measure.vars = patterns(cleaved = '^cleaved_', uncleaved = '^uncleaved_'),
      variable.name = 'pos'
    )
  
  levels(dt_TC$pos) <- c(1, 2, 3, 4, 5, 8)
  
  dt_TC[, FC := cleaved / (uncleaved + cleaved)]
  dt_TC[, rna := 'UC']
  
  fwrite(dt_TC, './output/datasets/cleavage_TC.csv.gz')
  rm(dt_TC)
  gc()
  }