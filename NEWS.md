# snplet 1.0.0

* Changed terminology from sample to barcode throughout codebase
* Changed get_sample_info to get_barcode_info with backwards compatibility alias
* Changed vireo_file parameter to be optional in import functions
* Changed vdj_file parameter to be optional in import_cellsnp
* Changed argument order in import_cellsnp (gene_annotation now before vdj_file)
* Added aggregate_count_df function
* Added samples file to example data
* Added comprehensive test suite
* Added annotation to tests
* Added example data files
* Added OTH count import functionality
* Added expression matrix export functionality
* Added metadata update utilities (add_barcode_metadata and add_snp_metadata)
* Added setters for barcode_info and snp_info slots
* Added handling for missing VDJ information when exporting SNPData
* Added improved error messages for missing clonotype data in to_expr_matri
* Fixed exporting barcodes functionality
* Fixed namespacing issues
* Fixed importing without vireo files
