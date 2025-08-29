check_awk_available <- function() {
    if (Sys.which("awk") == "") {
        stop("awk is not available. Please install awk to use this function.")
    }
}

check_bcftools_available <- function() {
    if (Sys.which("bcftools") == "") {
        stop("bcftools is not available. Please install bcftools to use this function.")
    }
}

check_samtools_available <- function() {
    if (Sys.which("samtools") == "") {
        stop("samtools is not available. Please install samtools to use this function.")
    }
}

check_gunzip_available <- function() {
    if (Sys.which("gunzip") == "") {
        stop("gunzip is not available. Please install gunzip to use this function.")
    }
}

gff_exons_to_bed <- function(gff_file, output_bed) {
    check_awk_available()

    gff_file <- normalizePath(gff_file, mustWork = TRUE)
    output_bed <- normalizePath(output_bed, mustWork = FALSE)

    if (!file.exists(gff_file)) {
        stop("GFF file does not exist: ", gff_file)
    }
    
    cmd <- sprintf(
        "awk '$3 == \"exon\" {print $1 \"\\t\" ($4-1) \"\\t\" $5}' %s > %s",
        shQuote(gff_file),
        shQuote(output_bed)
    )
    
    exit_code <- system(cmd)
    if (exit_code != 0) {
        stop("Failed to extract exons from GFF file")
    }
    
    if (!file.exists(output_bed)) {
        stop("Output BED file was not created")
    }
    
    return(invisible(output_bed))
}

filter_vcf_by_bed <- function(input_vcf, bed, output_vcf) {
    check_bcftools_available()
    
    input_vcf <- normalizePath(input_vcf, mustWork = TRUE)
    bed <- normalizePath(bed, mustWork = TRUE)
    output_vcf <- normalizePath(output_vcf, mustWork = FALSE)
    
    if (!file.exists(input_vcf)) {
        stop("Input VCF file does not exist: ", input_vcf)
    }
    
    if (!file.exists(bed)) {
        stop("BED file does not exist: ", bed)
    }
    
    # Determine if output should be compressed and use system2 to capture stderr directly
    is_compressed_output <- grepl("\\.gz$", output_vcf)
    if (is_compressed_output) {
        result <- system2(
            "bcftools", 
            args = c("view", "-R", bed, input_vcf, "-Oz", "-o", output_vcf),
            stderr = TRUE
        )
    } else {
        result <- system2(
            "bcftools", 
            args = c("view", "-R", bed, input_vcf, "-o", output_vcf),
            stderr = TRUE
        )
    }
    
    # Check if command failed
    if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
        error_output <- paste(result, collapse = "\n")
        
        # Check for common indexing errors
        if (grepl("Could not retrieve index file|could not load index", error_output)) {
            stop("bcftools failed: VCF file is not indexed. ",
                 "Please create an index with: bcftools index ", basename(input_vcf), "\n",
                 "Original error: ", error_output)
        } else {
            stop("bcftools view failed to filter VCF by BED regions\n",
                 "Error: ", error_output)
        }
    }
    
    if (!file.exists(output_vcf)) {
        stop("Output VCF file was not created")
    }
    
    return(invisible(output_vcf))
}

bam_depth_at_positions <- function(bam_file, vcf_file, output_file) {
    check_samtools_available()
    check_awk_available()
    check_gunzip_available()
    
    bam_file <- normalizePath(bam_file, mustWork = TRUE)
    vcf_file <- normalizePath(vcf_file, mustWork = TRUE)
    output_file <- normalizePath(output_file, mustWork = FALSE)
    
    if (!file.exists(bam_file)) {
        stop("BAM file does not exist: ", bam_file)
    }
    
    if (!file.exists(vcf_file)) {
        stop("VCF file does not exist: ", vcf_file)
    }
    
    temp_bed <- tempfile(fileext = ".bed")
    on.exit(unlink(temp_bed), add = TRUE)
    
    # Determine if VCF is compressed
    is_compressed <- grepl("\\.gz$", vcf_file)
    
    if (is_compressed) {
        bed_cmd <- sprintf(
            "gunzip -c %s | awk '!/^#/ {print $1 \"\\t\" ($2-1) \"\\t\" $2}' > %s",
            shQuote(vcf_file),
            shQuote(temp_bed)
        )
    } else {
        bed_cmd <- sprintf(
            "cat %s | awk '!/^#/ {print $1 \"\\t\" ($2-1) \"\\t\" $2}' > %s",
            shQuote(vcf_file),
            shQuote(temp_bed)
        )
    }
    
    bed_exit_code <- system(bed_cmd)
    if (bed_exit_code != 0) {
        stop("Failed to extract positions from VCF file")
    }
    
    depth_cmd <- sprintf(
        "samtools depth -b %s %s > %s",
        shQuote(temp_bed),
        shQuote(bam_file),
        shQuote(output_file)
    )
    
    depth_exit_code <- system(depth_cmd)
    if (depth_exit_code != 0) {
        stop("samtools depth failed to calculate depth at positions")
    }
    
    if (!file.exists(output_file)) {
        stop("Output depth file was not created")
    }
    
    return(invisible(output_file))
}