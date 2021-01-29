library(nucleR)

# load example data
# GRanges object
data(nucleosome_htseq)

# Filter reads and remove noise: discard the reads longer than 200bp 
# (threshold given to only keep mononucleosomes), remove noise due to MNase 
# efficiency by trimming reads to use only its central part (50bp around the dyad)
# Another GRanges object
reads_trim <- processReads(nucleosome_htseq, type="paired", fragmentLen=200, trim=50)

# Obtain the normalized coverage (the count of how many reads are mapped to each position, 
# divided by the total number of reads and multiplied by one milion)
# RleList object
cover_trim <- coverage.rpm(reads_trim)

# Smooth the coverage signal using the Fast Fourier Transformation
fft_ta <- filterFFT(cover_trim, pcKeepComp=0.01, showPowerSpec=TRUE)

# Detect peaks in the smoothed coverage which correspond to nucleosome 
# dyads and score them according to their fuzziness level
# GRanges object
peaks <- peakDetection(fft_ta, threshold="25%", score=TRUE, width=147)

