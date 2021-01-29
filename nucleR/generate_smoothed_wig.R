library(nucleR)
library(ggplot2)
library(GenomicRanges)

bam_file<-c("972-t0-20U-n1_S1.merged.bam")

reads <- readBAM(bam_file, type="paired")

cat("Number of reads fetched from file: ", length(reads),"\n")

# Process the paired end reads, but discard those with length > 200
reads_filt <- processReads(reads, type="paired", fragmentLen=200)

cat("Number of reads after filtering fragms with len>200: ", length(reads_filt),"\n")

# Process the reads, but now trim each read to 40bp around the dyad
reads_filt_trim <- processReads(reads, type="paired", fragmentLen=200, trim=40)

# Calculate the coverage, directly in reads per million (r.p.m.)
cover_orig <- coverage.rpm(reads_filt)
cover_trim <- coverage.rpm(reads_filt_trim)

# Compare both coverages, the dyad is much more clear in trimmed version
t1 <- as.vector(cover_orig[[1]])[1:5000]
t2 <- as.vector(cover_trim[[1]])[1:5000]

plot_data <- rbind(
    data.frame(
        x=seq_along(t1),
        y=t1,
    coverage="original"
    ),
    data.frame(
        x=seq_along(t2),
        y=t2,
        coverage="trimmed"
    )
)

qplot(x=x, y=y, color=coverage, data=plot_data, geom="line",
xlab="position", ylab="coverage")

# Let's try to call nucleosomes from the trimmed version
# First of all, let's remove some noise with FFT
# Power spectrum will be not plotted as this is not allowed with `SimpleRleList`
# The calculation will be done for all chromsomes here
cover_clean <- filterFFT(cover_trim, pcKeepComp=0.02, showPowerSpec=FALSE)

# now , let's check the Power spectrum for chr I,II,III
cover_clean_I <- filterFFT(cover_trim$I, pcKeepComp=0.02, showPowerSpec=TRUE)
cover_clean_II <- filterFFT(cover_trim$II, pcKeepComp=0.02, showPowerSpec=TRUE)
cover_clean_III <- filterFFT(cover_trim$III, pcKeepComp=0.02, showPowerSpec=TRUE)

# export wig files
export.wig(cover_clean_I, name="FFT coverages (972-t0-20U-n1_S1) for chrI",filepath="972-t0-20U-n1_S1.chrI.wig")

plot_data <- rbind(
    data.frame(
        x=1:4000,
        y=as.vector(cover_trim[[1]])[1:4000],
        coverage="noisy"
    ),
    data.frame(
        x=1:4000,
        y=as.vector(cover_clean[[1]])[1:4000],
        coverage="filtered"
    )
)

qplot(x=x, y=y, color=coverage, data=plot_data, geom="line", xlab="position", ylab="coverage")

# And how similar? Let's see the correlation
cor(cover_clean[[1]], as.vector(cover_trim[[1]]))

# Now it's time to call for peaks, first just as points
# See that the score is only a measure of the height of the peak
# Here, we only plot some peaks for chr I
peaks_I <- peakDetection(cover_clean_I, threshold="25%", score=TRUE)
plotPeaks(peaks_I[[1]][1:5000], cover_clean_I[1:5000], threshold="25%")

# Do the same as previously, but now we will create the nucleosome calls:
peaks_I <- peakDetection(cover_clean_I, width=147, threshold="25%", score=TRUE)
plotPeaks(peaks_I[[1]][1:5000], cover_clean_I[1:5000], threshold="25%")