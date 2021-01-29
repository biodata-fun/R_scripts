#
# SCRIPT to generate the smoothed wig file with the coverage values
# For smoothing the signal, the FFT from nucleR will be used
# AUTHOR: ernesto lowy (ernestolowy@gmail.com)
#
# USAGE: Rscript --vanilla generate_wig_files.R <BAM>
library(nucleR)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

bam_file<-c(args[1])

prefix<-tools::file_path_sans_ext(basename(bam_file))

reads <- readBAM(bam_file, type="paired")
cat("Number of reads fetched from file: ", length(reads),"\n")

# Process the reads, but now trim each read to 40bp around the dyad
reads_filt_trim <- processReads(reads, type="paired", fragmentLen=200, trim=40)
cat("Number of reads after filtering fragms with len>200: ", length(reads_filt_trim),"\n")

# Calculate the coverage, directly in reads per million (r.p.m.)
cover_trim <- coverage.rpm(reads_filt_trim)

# Per chromosome
# now , let's check the Power spectrum for chr I,II,III
#cover_clean_I <- filterFFT(cover_trim$I, pcKeepComp=0.02, showPowerSpec=TRUE)
#cover_clean_II <- filterFFT(cover_trim$II, pcKeepComp=0.02, showPowerSpec=TRUE)
#cover_clean_III <- filterFFT(cover_trim$III, pcKeepComp=0.02, showPowerSpec=TRUE)

# export wig files
#export.wig(cover_clean_I, name="FFT coverages (972-t0-20U-n1_S1) for chrI",filepath="972-t0-20U-n1_S1",chrom="I")
#export.wig(cover_clean_II, name="FFT coverages (972-t0-20U-n1_S1) for chrII",filepath="972-t0-20U-n1_S1",chrom="II")
#export.wig(cover_clean_III, name="FFT coverages (972-t0-20U-n1_S1) for chrIII",filepath="972-t0-20U-n1_S1",chrom="III")

# All chromosomes
cover_clean <- filterFFT(cover_trim, pcKeepComp=0.02, showPowerSpec=FALSE)
track_name<-paste("FFT coverages", prefix)
export.wig(cover_clean, name=track_name, filepath=prefix)