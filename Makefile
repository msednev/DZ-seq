SHELL := /bin/bash
RAW_DIR = ./input/raw_reads
CORE_DIR = ./output/processed_reads/core
SUBSTRATE_DIR = ./output/processed_reads/substrate
COUNTS_DIR = ./output/counts
DATASETS_DIR = ./output/datasets
LEFT_BINDING_ARM=GTGACGCGACTAGTTAC
RIGHT_BINDING_ARM=TTCATTCAGTTGGCGACTCC
CONNECTING_LOOP=TTCATTCAGTTGGCGACTCCRATAGACTGAAT


usage: Makefile
	@echo "usage: "
	@sed -n 's/^##//p' $<

RAW_READS := $(wildcard  $(RAW_DIR)/*.fq.gz)
CORES := $(patsubst $(RAW_DIR)/%.fq.gz, $(CORE_DIR)/%.fq, $(RAW_READS))

$(CORE_DIR)/%.fq: $(RAW_DIR)/%.fq.gz
	@echo "extracting deoxyribozyme cores from $*..."
	@mkdir -p $(CORE_DIR)/logs
		@cutadapt -g "$(LEFT_BINDING_ARM);e=0.2...$(RIGHT_BINDING_ARM);e=0.3"\
		-o $@ $< > $(CORE_DIR)/logs/$*.log 2> /dev/null

## extract-core:	extract deoxyribozyme cores from raw reads with cutadapt
extract-core: $(CORES)

SUBSTRATES := $(patsubst $(RAW_DIR)/%.fq.gz, $(SUBSTRATE_DIR)/%.fq, $(RAW_READS))

$(SUBSTRATE_DIR)/%.fq: $(RAW_DIR)/%.fq.gz
	@echo "extracting substrate sequences from $*..."
	@mkdir -p $(SUBSTRATE_DIR)/logs
	@cutadapt -g "$(CONNECTING_LOOP);e=0.3" -o $@  $<\
	       	> $(SUBSTRATE_DIR)/logs/$*.log 2> /dev/null

## extract-substrate:	extract sequences of RNA substrates from raw reads with cutadapt
extract-substrate: $(SUBSTRATES)

COUNTS := $(patsubst $(RAW_DIR)/%.fq.gz, $(COUNTS_DIR)/%.tsv, $(RAW_READS))

$(COUNTS_DIR)/%.tsv: $(CORE_DIR)/%.fq $(SUBSTRATE_DIR)/%.fq
	@echo "counting reads for $*..."
	@mkdir -p $(COUNTS_DIR)
	@paste $^ | awk 'NR%4==2 {print}' | datamash -s -g 1,2 count 1 | sort -nr -k3,3 > $@

$(DATASETS_DIR)/counts.csv.gz: $(COUNTS) scripts/process_counts.R
	@mkdir -p $(DATASETS_DIR)
	@echo "compiling counts dataset..."
	@Rscript scripts/process_counts.R



## counts:	compile dataset with count data
counts: $(DATASETS_DIR)/counts.csv.gz


$(DATASETS_DIR)/cleavage.csv.gz: $(DATASETS_DIR)/counts.csv.gz scripts/process_cleavage.R
	@echo "compiling cleavage dataset..."
	@Rscript scripts/process_cleavage.R
	@Rscript scripts/filter_cleavage.R

## cleavage:	compile dataset with cleavage data
cleavage: $(DATASETS_DIR)/cleavage.csv.gz


## clean:	delete all generated files
clean:
	rm -rf ./output


.PHONY: clean extract-core extract-substrate counts cleavage usage
