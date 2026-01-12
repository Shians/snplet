library(tidyverse)

human_assembly_report <- read_lines("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_assembly_report.txt")[-(1:61)] %>%
    paste(collapse = "\n") %>%
    read_tsv() %>%
    filter(`Sequence-Role` == "assembled-molecule")

mouse_assembly_report <- read_lines("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.20_GRCm38/GCF_000001635.20_GRCm38_assembly_report.txt")[-(1:41)] %>%
    paste(collapse = "\n") %>%
    read_tsv() %>%
    filter(`Sequence-Role` == "assembled-molecule")

human_assembly_report <- human_assembly_report %>%
    set_names(make.names(names(.), unique = TRUE)) %>%
    set_names(tolower(names(.)) %>%
        str_replace_all("\\.", "_") %>%
        str_remove_all("^x__")
    ) %>%
    select(sequence_name, genbank_accn, refseq_accn, ucsc_style_name) %>%
    rename(
        genbank_human = genbank_accn,
        refseq_human = refseq_accn
    )

mouse_assembly_report <- mouse_assembly_report %>%
    set_names(make.names(names(.), unique = TRUE)) %>%
        set_names(tolower(names(.)) %>%
        str_replace_all("\\.", "_") %>%
        str_remove_all("^x__")
    ) %>%
    select(sequence_name, genbank_accn, refseq_accn, ucsc_style_name) %>%
    rename(
        genbank_mouse = genbank_accn,
        refseq_mouse = refseq_accn
    )

full_join(human_assembly_report, mouse_assembly_report) %>%
    rename(
        numeric = sequence_name,
        ucsc = ucsc_style_name
    ) %>%
    relocate(ucsc, .before = numeric) %>%
    write_tsv("inst/extdata/chr_name_table.tsv")
