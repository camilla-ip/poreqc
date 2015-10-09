# Usage:
# Rscript ont_runqc_readstats.R readstats.txt callstats.txt activechannels.txt filldelay.txt rundesc outdir outprefix imgszconst
# Input file is a tab-separated file with one header and columns as produced by extract_readstats.py and extract_callstats.py

# =========================================================================== #
# Constants                                                                   #
# =========================================================================== #

# Palettes from http://blog.mollietaylor.com/2012/10/color-blindness-and-palette-choice.html
# 1=green, 2=orange, 3=purple, 4=pink, 5=green, 6=yellow, 7=brown, 8=grey
alpha <- 0.5
cbPaletteLt <- c(rgb(0.40,0.76,0.65, alpha), rgb(0.99,0.55,0.38, alpha), rgb(0.55,0.63,0.80, alpha), rgb(0.91,0.54,0.76, alpha),
                 rgb(0.65,0.85,0.33, alpha), rgb(1.00,0.85,0.18, alpha), rgb(0.90,0.77,0.58, alpha), rgb(0.70,0.70,0.70, alpha))
cbPaletteDk <- c(rgb(0.11,0.62,0.47, alpha), rgb(0.85,0.37,0.01, alpha), rgb(0.46,0.44,0.70, alpha), rgb(0.91,0.16,0.54, alpha),
                 rgb(0.40,0.65,0.12, alpha), rgb(0.90,0.67,0.01, alpha), rgb(0.65,0.46,0.11, alpha), rgb(0.26,0.26,0.26, alpha))

gridmajcol <- rgb(0.78,0.82,0.80, alpha)
gridmincol <- rgb(0.78,0.82,0.80, alpha)

ticksz <- 0.6
labelsz <- 0.5
titlesz <- 0.5
legendsz <- 0.3
plotdpi <- 300

kb_number <- 1e+3
mb_number <- 1e+6
gb_number <- 1e+9
tb_number <- 1e+12
read_length_bucketsz <- 1000.0
read_length_bucketsz_small <- 100.0
channel_range <- 1:512

# =========================================================================== #
# Command-line arguments                                                      #
# =========================================================================== #

args <- commandArgs(trailingOnly = TRUE)
inreadstatspath <- args[1]
incallstatspath <- args[2]
inactivepath <- args[3]
indelaypath <- args[4]
rundesc <- args[5]
outdir <- args[6]
outprefix <- args[7]
imgszconst <- as.numeric(args[8])

debugmode <- FALSE
if (debugmode) {
inreadstatspath <- '/Users/ip/work/ont/20141019-qc-v0.2.3/debug_cDNA_50MP_SSII_freshPSM/readstats.txt'
incallstatspath <- '/Users/ip/work/ont/20141019-qc-v0.2.3/debug_cDNA_50MP_SSII_freshPSM/callstats.txt'
inactivepath <- '/Users/ip/work/ont/20141019-qc-v0.2.3/debug_cDNA_50MP_SSII_freshPSM/activechannels.txt'
indelaypath <- '/Users/ip/work/ont/20141019-qc-v0.2.3/debug_cDNA_50MP_SSII_freshPSM/filldelay.txt'
rundesc <- 'htgsaureus2'
outdir <- '/Users/ip/work/ont/20141019-qc-v0.2.3/debug_cDNA_50MP_SSII_freshPSM'
outprefix <- 'htgsaureus2'
imgszconst <- 50
}
pngpdfdir <- sprintf("%s/pngpdf", outdir)
pnghtmldir <- sprintf("%s/pnghtml", outdir)

imgrectwidth <- as.integer(imgszconst * 2)
imgrectheight <- as.integer(imgszconst * 0.75)
imgsquarewidth <- as.integer(imgszconst)
imgsquareheight <- as.integer(imgszconst * 0.8)

library(ggplot2)
library(plyr)
library(reshape2)
library(grid)

readdata <- read.table(inreadstatspath, header=TRUE, sep="\t")
calldata <- read.table(incallstatspath, header=TRUE, sep="\t")
activechannels <- read.table(inactivepath, header=TRUE, sep="\t")
filldelay <- read.table(indelaypath, header=TRUE, sep="\t")

if (! file.exists(pngpdfdir)){
    dir.create(pngpdfdir)
}

# =========================================================================== #
# runtime plots                                                               #
# =========================================================================== #

pngpath <- sprintf('%s/%s_runtime_readcompletionrate.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
D <- data.frame(aggregate(rep(1, length(readdata$read_end_hours_since_expt_start_qtrhrbucketname)),
    by=list(readdata$read_end_hours_since_expt_start_qtrhrbucketname), FUN=sum))
timebucket <- D$Group.1
readcount <- D$x
totalreadcount <- formatC(max(D$x), format="d", big.mark=',')
plottitle <- sprintf("%s : Read completion rate (max reads completed in a 15m interval = %s)", rundesc, totalreadcount)
xlabel <- 'Time (h)'
ylabel <- 'Completed reads (K)'
p <- ggplot(D, aes(x=timebucket, y=readcount/1000.0)) +
    geom_line(stat="identity", colour=cbPaletteDk[3], fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8]) }
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_runtime_eventcompletionrate.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
D <- data.frame(aggregate(readdata$event_count, by=list(readdata$read_end_hours_since_expt_start_qtrhrbucketname), FUN=sum))
timebucket <- D$Group.1
eventcount <- D$x/mb_number
maxeventcount <- formatC(max(D$x), format="d", big.mark=',')
plottitle <- sprintf("%s : Event completion rate (max events completed in a 15m interval = %s)", rundesc, maxeventcount)
xlabel <- 'Time (h)'
ylabel <- 'Completed events (M)'
p <- ggplot(D, aes(x=timebucket, y=eventcount)) +
    geom_line(stat="identity", colour=cbPaletteDk[3], fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8])
}
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_runtime_readcompletioncount.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
D <- data.frame(aggregate(rep(1, length(readdata$read_end_hours_since_expt_start_qtrhrbucketname)),
    by=list(readdata$read_end_hours_since_expt_start_qtrhrbucketname), FUN=sum))
timebucket <- D$Group.1
readcount <- cumsum(D$x)/1000.0
totalreadscompleted <- formatC(max(readcount), format="d", big.mark=',')
plottitle <- sprintf("%s : Read yield (total = %s)", rundesc, totalreadscompleted)
xlabel <- 'Time (h)'
ylabel <- 'Completed reads (K)'
p <- ggplot(D, aes(x=timebucket, y=readcount)) +
    geom_line(colour=cbPaletteDk[3], fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
#    xlim(0, max(timebucket)) +
#    ylim(0, ceiling(max(readcount) / 100.0) * 100.0) +
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8])
}
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_runtime_eventcompletioncount.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
G <- data.frame(aggregate(readdata$event_count,
    by=list(readdata$read_end_hours_since_expt_start_qtrhrbucketname), FUN=sum))
timebucket <- G$Group.1
eventcount <- cumsum(G$x) / mb_number
totalyield = formatC(max(cumsum(G$x)), format="d", big.mark=',')
plottitle <- sprintf("%s : Event yield (total = %s events)", rundesc, totalyield)
xlabel <- 'Time (h)'
ylabel <- 'Events (M)'
p <- ggplot(G, aes(x=timebucket, y=eventcount)) +
    geom_line(colour=cbPaletteDk[3], fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8]) }
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_runtime_activechannels.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
timebucket <- activechannels$time_bucket_since_expt_start_hrs[activechannels$time_bucket_since_expt_start_hrs > 0]
numactive <- activechannels$num_active_channels[activechannels$time_bucket_since_expt_start_hrs > 0]
active_channels <- data.frame(timebucket, numactive)
maxactive <- formatC(max(active_channels), format="d", big.mark=',')
numuniquechannels <- formatC(length(unique(readdata$channel_number)), format="d", big.mark=",")
plottitle <- sprintf("%s : Active channels (max active = %s, max concurrently = %s)", rundesc, numuniquechannels, maxactive)
xlabel <- 'Time (h)'
ylabel <- 'Active channels'
p <- ggplot(active_channels, aes(x=timebucket, y=numactive)) +
    geom_line(colour=cbPaletteDk[3], fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8]) }
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_runtime_channelspeedovertime.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
min_eventsS = sprintf("%d", min_events)
max_events <- max(100000, calldata$max_events[1])
max_eventsS = formatC(max_events, format="d", big.mark=',')
event_count <- readdata$event_count[readdata$event_count >= min_events & readdata$event_count <= max_events]
read_end_time <- readdata$read_end_hours_since_expt_start[readdata$event_count >= min_events & readdata$event_count <= max_events]
read_duration_time <- readdata$read_duration_seconds[readdata$event_count >= min_events & readdata$event_count <= max_events]
porespeed <- event_count / read_duration_time
plotdata <- data.frame(read_end_time=read_end_time, porespeed=porespeed)
plotdata$read_end_time[which(is.nan(plotdata$read_end_time))] = NA
plotdata$read_end_time[which(plotdata$read_end_time==Inf)] = NA
plotdata$porespeed[which(is.nan(plotdata$porespeed))] = NA
plotdata$porespeed[which(plotdata$porespeed==Inf)] = NA
H <- lm(plotdata$porespeed ~ plotdata$read_end_time)
slope <- H$coefficients[2]
intercept <- H$coefficients[1]
smallnum_excluded <- length(readdata$event_count[readdata$event_count < min_events])
bignum_excluded <- length(readdata$event_count[readdata$event_count > max_events])
max_eventsS = formatC(max_events, format="d", big.mark=',')
plottitle <- sprintf("%s : Channel speed (excl. %d reads len<%s and %d reads with len>%s)", rundesc, smallnum_excluded, min_eventsS, bignum_excluded, max_eventsS)
xlabel <- 'Time (h)'
ylabel <- 'Channel speed (events/s)'
p <- ggplot(plotdata, aes(x=read_end_time, y=porespeed)) +
    geom_point(shape=16, size=0.25, colour=cbPaletteLt[2], fill=cbPaletteLt[2]) +
    #geom_smooth(method=lm, colour=cbPaletteDk[4], fill=cbPaletteLt[4]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8])
}
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_runtime_chiptemp.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
time <- calldata$read_start_hours_since_expt_start
temperature <- c(calldata$asic_temp, calldata$heatsink_temp)
temptype <- c(rep("chip", length(calldata$asic_temp)), rep("heat-sink", length(calldata$heatsink_temp)))
plotdata <- data.frame(time=time, temp=temperature, temptype=temptype)
plottitle <- sprintf("%s : Flow-cell temperature", rundesc)
xlabel <- 'Time (h)'
ylabel <- 'Temperature (degrees C)'
p <- ggplot(plotdata, aes(x=time, y=temperature, group=temptype, colour=temptype)) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    theme(legend.title=element_blank()) +
    scale_fill_manual(
        values=c(cbPaletteDk[2], cbPaletteLt[7]),
        name="",
        breaks=c("chip", "heat-sink")) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(0.1, "cm"), legend.key.height = unit(0.15, "cm")) +
    theme(legend.position=c(0.85,0.6)) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_runtime_2Dmeanbqovertime.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
time <- calldata$read_start_hours_since_expt_start[calldata$twod_fastq_seqlen > 0]
meanbq <- calldata$twod_fastq_bqmean[calldata$twod_fastq_seqlen > 0]
if (length(time) == 0) {
    time <- calldata$read_start_hours_since_expt_start
    meanbq <- rep(0, length(time))
}
plotdata <- data.frame(time=time, meanbq=meanbq)
H <- lm(plotdata$meanbq ~ plotdata$time)
slope <- H$coefficients[2]
intercept <- H$coefficients[1]
plottitle <- sprintf("%s : 2D-read base quality over time (linear regression with 95%% confidence interval, bq = %.3f + %.3f*time)", rundesc, intercept, slope)
xlabel <- 'Time (h)'
ylabel <- 'Mean base quality of read'
p <- ggplot(plotdata, aes(x=time, y=meanbq)) +
    geom_point(shape=16, size=0.25, colour=cbPaletteLt[5], fill=cbPaletteLt[5]) +
    geom_smooth(method=lm, colour=cbPaletteDk[4], fill=cbPaletteLt[4]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8]) }
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

print(sprintf("Generating %s", pngpath))
tmtwod <- calldata$read_start_hours_since_expt_start[calldata$twod_fastq_seqlen > 0]
tmtemp <- calldata$read_start_hours_since_expt_start[calldata$template_fastq_seqlen > 0]
tmcomp <- calldata$read_start_hours_since_expt_start[calldata$complement_fastq_seqlen > 0]
bqtwod <- calldata$twod_fastq_bqmean[calldata$twod_fastq_seqlen > 0]
bqtemp <- calldata$template_fastq_bqmean[calldata$template_fastq_seqlen > 0]
bqcomp <- calldata$complement_fastq_bqmean[calldata$complement_fastq_seqlen > 0]
type_all <- c(rep("2D", length(tmtwod)), rep("temp", length(tmtemp)), rep("comp", length(tmcomp)))
time_all <- c(tmtwod, tmtemp, tmcomp)
bq_all <- c(bqtwod, bqtemp, bqcomp)
twod_median <- median(bqtwod)
temp_median <- median(bqtemp)
comp_median <- median(bqcomp)
plotdata <- data.frame(type=type_all, time=time_all, bq=bq_all)
plotdata$type <- factor(plotdata$type, levels=c("2D", "temp", "comp"), labels=c("2D", "temp", "comp"))
plottitle <- sprintf("%s : Read base qualities over time (median 2D=%.1f, temp=%.1f, comp=%.1f)",
    rundesc, twod_median, temp_median, comp_median)
xlabel <- 'Time (h)'
ylabel <- 'Mean base quality'
p <- ggplot(plotdata, aes(x = time, y = bq, color = type, shape = type)) +
    geom_point(shape=16, size=0.25) +
    scale_colour_manual(values=c(cbPaletteDk[5], cbPaletteDk[7], cbPaletteDk[6])) +
    guides(colour = guide_legend(override.aes = list(size=1))) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid")) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(0.1, "cm"), legend.key.height = unit(0.15, "cm")) +
    theme(legend.position=c(0.8,0.8))
if (max(timebucket) > 24) {
    p + geom_vline(aes(xintercept = 24), linetype = "dashed", size=0.25, colour=cbPaletteDk[8]) }
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

# WARNING: uses plotdata variable from previous plot
pngpath <- sprintf('%s/%s_runtime_meanbqboxplots.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
plottitle <- sprintf("%s : Mean base qualities", rundesc)
xlabel <- ''
ylabel <- 'Mean base quality'
colourlist <- c()
filllist <- c()
if ("2D" %in% unique(plotdata$type)) {
    colourlist <- c(colourlist, cbPaletteDk[5])
    filllist <- c(filllist, cbPaletteLt[5])
}
if ("temp" %in% unique(plotdata$type)) {
    colourlist <- c(colourlist, cbPaletteDk[7])
    filllist <- c(filllist, cbPaletteLt[7])
}
if ("comp" %in% unique(plotdata$type)) {
    colourlist <- c(colourlist, cbPaletteDk[6])
    filllist <- c(filllist, cbPaletteLt[6])
}
#if (length(intersect(c("2D"), unique(plotdata$type))) == 0) {
#    plotdata <- insertRow(plotdata, c("2D", 0, 0), 1)
#}
#if (length(intersect(c("temp"), unique(plotdata$type))) == 0) {
#    plotdata <- insertRow(plotdata, c("temp", 0, 0), 1)
#}
#if (length(intersect(c("comp"), unique(plotdata$type))) == 0) {
#    plotdata <- insertRow(plotdata, c("comp", 0, 0), 1)
#}
#    geom_boxplot(colour=c(cbPaletteDk[5], cbPaletteDk[7], cbPaletteDk[6]),
#        fill=c(cbPaletteLt[5], cbPaletteLt[7], cbPaletteLt[6]),
#        outlier.colour=cbPaletteLt[8], outlier.size=1.0, lwd=0.2) +
p <- ggplot(plotdata, aes(factor(type), bq)) +
    geom_boxplot(colour=colourlist,
        fill=filllist,
        outlier.colour=cbPaletteLt[8], outlier.size=1.0, lwd=0.2) +
    scale_fill_manual(values=c(cbPaletteLt[5], cbPaletteLt[7], cbPaletteLt[6])) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title.y=element_text(size=rel(labelsz))) +
    theme(axis.title.x=element_blank()) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgsquarewidth, height=imgsquareheight, units="mm", dpi=plotdpi)

# =========================================================================== #
# readlen plots                                                               #
# =========================================================================== #

pngpath <- sprintf('%s/%s_readlen_histall.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS = sprintf("%.1f", min_events/kb_number)
max_eventsS = formatC(max_events/1000.0, format="d", big.mark=',')
eventmask <- readdata$event_count >= min_events & readdata$event_count <= max_events
maxusedreadlen <- formatC(max(readdata$event_count[eventmask]), format="d", big.mark=',')
maxreadlen <- formatC(max(readdata$event_count), format="d", big.mark=',')
A <- data.frame(event_count = readdata$event_count[eventmask]/1000.0)
plottitle <- sprintf("%s : Event-read length distribution (max = %s)", rundesc, maxusedreadlen)
xlabel <- sprintf('Callable read length (K, range=[%sK , %sK])', min_eventsS, max_eventsS)
ylabel <- 'Frequency'
p <- ggplot(A, aes(x=A$event_count)) +
    geom_histogram(fill=cbPaletteDk[3]) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgsquarewidth, height=imgsquareheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_readlen_histlower.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS = sprintf("%.1f", min_events/1000.0)
max_eventsS = formatC(max_events, format="d", big.mark=',')
max_thisplot <- median(readdata$event_count) * 3.0
max_thisplotS <- sprintf("%.3f", max_thisplot/kb_number)
eventmask <- readdata$event_count >= min_events & readdata$event_count <= max_events & readdata$event_count < max_thisplot
maxusedreadlen <- formatC(max(readdata$event_count[eventmask]), format="d", big.mark=',')
maxreadlen <- formatC(max(readdata$event_count), format="d", big.mark=',')
A <- data.frame(event_count = readdata$event_count[eventmask]/kb_number)
plottitle <- sprintf("%s : Event-read lengths (lower values only)", rundesc, maxusedreadlen)
xlabel <- sprintf('Callable read length (K, median is dotted)', min_eventsS, max_thisplotS)
ylabel <- 'Frequency'
p <- ggplot(A, aes(x=A$event_count)) +
    geom_histogram(fill=cbPaletteDk[3]) +
    geom_vline(aes(xintercept = as.integer(median(readdata$event_count)/kb_number)),
        linetype = "dashed", size=0.25, colour=cbPaletteDk[8]) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_readlen_histupper.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS <- sprintf("%.1f", min_events/1000.0)
max_eventsS <- formatC(max_events/1000.0, format="d", big.mark=',')
min_thisplot <- median(readdata$event_count) * 3.0
min_thisplotS <- sprintf("%.3f", min_thisplot/1000.0)
eventmask <- readdata$event_count >= min_events & readdata$event_count <= max_events & readdata$event_count > min_thisplot
maxusedreadlen <- formatC(max(readdata$event_count[eventmask]), format="d", big.mark=',')
maxreadlen <- formatC(max(readdata$event_count), format="d", big.mark=',')
A <- data.frame(event_count = readdata$event_count[eventmask]/1000.0)
plottitle <- sprintf("%s : Event-read lengths (upper values only)", rundesc, maxusedreadlen)
xlabel <- sprintf('Callable read length (K, median is dotted)', min_eventsS, min_thisplotS)
ylabel <- 'Frequency'
p <- ggplot(A, aes(x=A$event_count)) +
    geom_histogram(fill=cbPaletteDk[3], binwidth=1) +
    geom_vline(aes(xintercept = as.integer(median(readdata$event_count)/1000.0)),
        linetype = "dashed", size=0.25, colour=cbPaletteDk[8]) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgsquarewidth, height=imgsquareheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_readlen_yield.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS <- sprintf("%.1f", min_events/1000.0)
max_eventsS <- formatC(max_events/1000.0, format="d", big.mark=',')
plottitle <- sprintf("%s : Callable event yield by read length", rundesc)
yield <- data.frame(aggregate(readdata$event_count[readdata$event_count >= min_events & readdata$event_count <= max_events],
    by=list(ceiling(readdata$event_count[readdata$event_count >= min_events & readdata$event_count <= max_events] / read_length_bucketsz) * read_length_bucketsz), FUN=sum))
readlengroup <- yield$Group.1 / read_length_bucketsz
totalevents <- yield$x
conversion <- 1
units <- ''
if (max(totalevents) > kb_number) {
    conversion = kb_number
    units = ' (K)' }
if (max(totalevents) > mb_number) {
    conversion = mb_number
    units = ' (M)' }
if (max(totalevents) > gb_number) {
    conversion = gb_number
    units = ' (G)' }
totalevents <- yield$x / conversion
xlabel <- sprintf('Read length increment (K) (range = [%s, %s])', min_eventsS, max_eventsS)
ylabel <- sprintf('Event yield%s', units)
p <- ggplot(yield, aes(x=readlengroup, y=totalevents)) +
    geom_bar(stat="identity", fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_readlen_yieldlower.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS <- sprintf("%.1f", min_events/kb_number)
max_eventsS <- formatC(max_events/1000.0, format="d", big.mark=',')
max_thisplot <- median(readdata$event_count[readdata$event_count >= min_events & readdata$event_count <= max_events]) * 10.0
max_thisplotS <- sprintf("%.3f", max_thisplot/kb_number)
plottitle <- sprintf("%s : Callable event yield by read length (lower values only)", rundesc)
eventmask <- readdata$event_count >= min_events & readdata$event_count <= max_events & readdata$event_count < max_thisplot
yield <- data.frame(aggregate(readdata$event_count[eventmask], by=list(ceiling(readdata$event_count[eventmask] / read_length_bucketsz_small) * read_length_bucketsz_small), FUN=sum))
readlengroup <- yield$Group.1 / read_length_bucketsz_small / 10.0
totalevents <- yield$x
conversion <- 1
units <- ''
if (max(totalevents) > kb_number) {
    conversion = kb_number
    units = ' (K)' }
if (max(totalevents) > mb_number) {
    conversion = mb_number
    units = ' (M)' }
if (max(totalevents) > gb_number) {
    conversion = gb_number
    units = ' (G)' }
totalevents <- yield$x / conversion
xlabel <- sprintf('Read length increment (K) (range = [%s , %s]', min_eventsS, max_thisplotS)
ylabel <- sprintf('Event yield%s', units)
p <- ggplot(yield, aes(x=readlengroup, y=totalevents)) +
    geom_bar(stat="identity", fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

# =========================================================================== #
# basecalled plots                                                            #
# =========================================================================== #

pngpath <- sprintf('%s/%s_basecall_readyield.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
object <-c(
    rep("Event files", length(readdata$fast5_filename[readdata$event_count >= min_events & readdata$event_count <= max_events])),
    rep("Call files", length(calldata$fast5_filename[calldata$event_count >= min_events & calldata$event_count <= max_events])),
    rep("Template", length(calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0])),
    rep("Complement", length(calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0])),
    rep("2D", length(calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0])))
count <- c(
    rep(1, length(readdata$event_count[readdata$event_count >= min_events & readdata$event_count <= max_events])),
    rep(1, length(calldata$fast5_filename[calldata$event_count >= min_events & calldata$event_count <= max_events])),
    rep(1, length(calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0])),
    rep(1, length(calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0])),
    rep(1, length(calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0])))
lab <- c(
    rep("Event files", length(readdata$fast5_filename[readdata$event_count >= min_events & readdata$event_count <= max_events])),
    rep("Call files", length(calldata$fast5_filename[calldata$event_count >= min_events & calldata$event_count <= max_events])),
    rep("Template", length(calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0])),
    rep("Complement", length(calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0])),
    rep("2D", length(calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0])))
plottitle <- sprintf("%s : Called read yield from callable event-reads", rundesc)
xlabel <- ''
ylabel <- 'Read count (K)'
B <- data.frame(aggregate(count/kb_number, by=list(object), FUN=sum))
B$Group.1 <- factor(B$Group.1, levels=c("Event files", "Call files", "Template", "Complement", "2D"), labels=c("Event files", "Call files", "Template", "Complement", "2D"))
object_list <- B$Group.1
count_total <- B$x
p <- ggplot(B, aes(x=object_list, y=count_total, label=object_list, fill=object_list)) +
    geom_bar(stat="identity") +
    theme(legend.position = "none") +
    scale_fill_manual(
        values=c(cbPaletteDk[3], cbPaletteLt[3], cbPaletteDk[7], cbPaletteDk[6], cbPaletteDk[5]),
        name="",
        breaks=c("Event files", "Call files", "Template", "Complement", "2D")) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.1, "cm"), legend.key.height = unit(.15, "cm")) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgsquarewidth, height=imgsquareheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_baseyield.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
object <-c(
    rep("Event files", length(readdata$fast5_filename[readdata$event_count >= min_events & readdata$event_count <= max_events])),
    rep("Call files", length(calldata$fast5_filename[calldata$event_count >= min_events & calldata$event_count <= max_events])),
    rep("Template", length(calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0])),
    rep("Complement", length(calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0])),
    rep("2D", length(calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0])))
count <- c(
    readdata$event_count[readdata$event_count >= min_events & readdata$event_count <= max_events],
    calldata$event_count[calldata$event_count >= min_events & calldata$event_count <= max_events],
    calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0],
    calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0],
    calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0])
B <- data.frame(aggregate(count/mb_number, by=list(object), FUN=sum))
B$Group.1 <- factor(B$Group.1, levels=c("Event files", "Call files", "Template", "Complement", "2D"),
    labels=c("Event files", "Call files", "Template", "Complement", "2D"))
object_list <- B$Group.1
count_total <- B$x
plottitle <- sprintf("%s : Called base recovery from callable event-reads", rundesc)
xlabel <- ''
ylabel <- 'Event or base count (M)'
p <- ggplot(B, aes(x=object_list, y=count_total, label=object_list, fill=object_list)) +
    geom_bar(stat="identity") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    scale_fill_manual(
        values=c(cbPaletteDk[3], cbPaletteLt[3], cbPaletteDk[7], cbPaletteDk[6], cbPaletteDk[5]),
        name="",
        breaks=c("Event files", "Call files", "Template", "Complement", "2D")) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.1, "cm"), legend.key.height = unit(.15, "cm")) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgsquarewidth, height=imgsquareheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_eventto2Dreadlen.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
#D <- data.frame(eventcount=calldata$event_count[calldata$event_count >= min_events & calldata$event_count <= max_events]/kb_number,
#    seqlen=calldata$twod_fastq_seqlen[calldata$event_count >= min_events & calldata$event_count <= max_events]/kb_number)
if (length(calldata$event_count[calldata$event_count >= min_events & calldata$event_count <= max_events & calldata$twod_fastq_seqlen]) == 0) {
D <- data.frame(eventcount=c(0), seqlen=c(0))
} else {
D <- data.frame(eventcount=calldata$event_count[calldata$event_count >= min_events & calldata$event_count <= max_events & calldata$twod_fastq_seqlen]/kb_number,
    seqlen=calldata$twod_fastq_seqlen[calldata$event_count >= min_events & calldata$event_count <= max_events & calldata$twod_fastq_seqlen > 0]/kb_number)
}
H <- lm(D$seqlen ~ D$eventcount)
slope <- H$coefficients[2]
intercept <- H$coefficients[1]
plottitle <- sprintf("%s : Event and 2D read lengths\n(linear regression with 95%% confidence interval, eventlen = %.3f*2Dreadlen + %.3f)", rundesc, slope, intercept)
xlabel <- 'Event read length (K)'
ylabel <- '2D read length (K)'
p <- ggplot(D, aes(x=eventcount, y=seqlen)) +
    geom_point(shape=16, size=0.25, colour=cbPaletteLt[5], fill=cbPaletteLt[5]) +
    geom_smooth(method=lm, colour=cbPaletteDk[4], fill=cbPaletteLt[4]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgsquarewidth, height=imgsquareheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_readlendistall.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
tempreadlen <- calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0]/kb_number
compreadlen <- calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0]/kb_number
twodreadlen <- calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0]/kb_number
readtypelist <- c(rep("Temp", length(tempreadlen)), rep("Comp", length(compreadlen)), rep("2D", length(twodreadlen)))
readlengthlist <- c(tempreadlen, compreadlen, twodreadlen)
if (length(calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0]) == 0) {
    max2DlenS <- formatC(0, format="d", big.mark=',')
} else {
    max2DlenS <- formatC(max(calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0]), format="d", big.mark=',')
}
if (length(calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0]) == 0) {
    maxTemplenS <- formatC(0, format="d", big.mark=',')
} else {
    maxTemplenS <- formatC(max(calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0]), format="d", big.mark=',')
}
if (length(calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0]) == 0) {
    maxComplenS <- formatC(0, format="d", big.mark=',')
} else {
    maxComplenS <- formatC(max(calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0]), format="d", big.mark=',')
}
if (length(readtypelist) == 0 || length(readlengthlist) == 0) {
    D <- data.frame(readtype = c(0), readlen = c(0))
    D$readtype <- factor(D$readtype, levels=c("2D", "Temp", "Comp"), labels=c("2D", "Temp", "Comp"))
} else {
    D <- data.frame(readtype = readtypelist, readlen = readlengthlist)
    D$readtype <- factor(D$readtype, levels=c("2D", "Temp", "Comp"), labels=c("2D", "Temp", "Comp"))
}
plottitle <- sprintf("%s : Called read lengths (max 2D=%s, T=%s, C=%s)", rundesc, max2DlenS, maxTemplenS, maxComplenS)
xlabel <- 'Read length (K)'
ylabel <- 'Count'
p <- ggplot(D, aes(x=D$readlen, fill=readtype)) +
    geom_histogram(position="dodge") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    scale_fill_manual(
        values=c(cbPaletteDk[5], cbPaletteDk[7], cbPaletteDk[6]), 
        name="",
        breaks=c("2D", "Temp", "Comp"),
        labels=c("2D", "Temp", "Comp")) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text = element_text(size=rel(ticksz))) +
    theme(axis.title = element_text(size=rel(labelsz))) +
    theme(plot.title = element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.1, "cm"), legend.key.height = unit(.15, "cm")) +
    theme(legend.position = c(0.75,0.6)) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.25, "cm"), legend.key.height = unit(.25, "cm")) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_readlendist2D.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
read_length <- calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0]/kb_number
if (length(read_length) == 0) {
    D <- data.frame(readlen = c(0))
    max2DlenS <- formatC(0, format="d", big.mark=',')
} else {
    D <- data.frame(readlen = read_length)
    max2DlenS <- formatC(max(calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0]), format="d", big.mark=',')
}
plottitle <- sprintf("%s : Called 2D read lengths (max=%s)", rundesc, max2DlenS)
xlabel <- 'Read length (K)'
ylabel <- 'Count'
p <- ggplot(D, aes(x=D$readlen)) +
    geom_histogram(fill=cbPaletteDk[5], binwidth=1) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text = element_text(size=rel(ticksz))) +
    theme(axis.title = element_text(size=rel(labelsz))) +
    theme(plot.title = element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.1, "cm"), legend.key.height = unit(.15, "cm")) +
    theme(legend.position = c(0.75,0.6)) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.25, "cm"), legend.key.height = unit(.25, "cm")) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid")) 
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_readlendisttemp.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
read_length <- calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0]/kb_number
D <- data.frame(readlen = read_length)
maxTemplenS <- formatC(max(calldata$template_fastq_seqlen[calldata$template_fastq_seqlen > 0]), format="d", big.mark=',')
plottitle <- sprintf("%s : Called template read lengths (max=%s)", rundesc, maxTemplenS)
xlabel <- 'Read length (K)'
ylabel <- 'Count'
p <- ggplot(D, aes(x=D$readlen)) +
    geom_histogram(fill=cbPaletteDk[7], binwidth=1) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text = element_text(size=rel(ticksz))) +
    theme(axis.title = element_text(size=rel(labelsz))) +
    theme(plot.title = element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.1, "cm"), legend.key.height = unit(.15, "cm")) +
    theme(legend.position = c(0.75,0.6)) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.25, "cm"), legend.key.height = unit(.25, "cm")) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_readlendistcomp.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
read_length <- calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0]/kb_number
if (length(read_length) == 0) {
    D <- data.frame(readlen = c(0))
    maxComplenS <- formatC(0, format="d", big.mark=',')
} else {
    D <- data.frame(readlen = read_length)
    maxComplenS <- formatC(max(calldata$complement_fastq_seqlen[calldata$complement_fastq_seqlen > 0]), format="d", big.mark=',')
}
plottitle <- sprintf("%s : Called complement read lengths (max=%s)", rundesc, maxTemplenS)
xlabel <- 'Read length (K)'
ylabel <- 'Count'
p <- ggplot(D, aes(x=D$readlen)) +
    geom_histogram(fill=cbPaletteDk[6], binwidth=1) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text = element_text(size=rel(ticksz))) +
    theme(axis.title = element_text(size=rel(labelsz))) +
    theme(plot.title = element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.1, "cm"), legend.key.height = unit(.15, "cm")) +
    theme(legend.position = c(0.75,0.6)) +
    theme(legend.text = element_text(size=rel(legendsz))) +
    theme(legend.key.width = unit(.25, "cm"), legend.key.height = unit(.25, "cm")) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_channelspeedtemp2comp.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
x <- calldata$template_num_events[calldata$template_num_events > 0]
y <- calldata$template_events_duration[calldata$template_num_events > 0]
template_pore_speed <- x / y
x <- calldata$complement_num_events[calldata$template_num_events > 0]
y <- calldata$complement_events_duration[calldata$template_num_events > 0]
complement_pore_speed <- x / y
if (length(complement_pore_speed) > 0 && sum(!is.nan(complement_pore_speed)) == 0) {
    plotdata <- data.frame(tps=c(0), cps=c(0))
    slope <- 0
    intercept <- 0
} else {
    plotdata <- data.frame(tps=template_pore_speed, cps=complement_pore_speed)
    H <- lm(plotdata$cps ~ plotdata$tps)
    slope <- H$coefficients[2]
    intercept <- H$coefficients[1]
}
plottitle <- sprintf("%s : Channel speed of template and complement\n(linear regression with 95%% confidence interval, compspeed = %.3f*tempspeed + %.3f)", rundesc, slope, intercept)
xlabel <- 'Template channel speed (events/s)'
ylabel <- 'Complement channel speed (events/s)'
p <- ggplot(plotdata, aes(x=tps, y=cps)) +
    geom_point(shape=16, size=0.25, colour=cbPaletteLt[2], fill=cbPaletteLt[2]) +
    geom_smooth(method=lm, colour=cbPaletteDk[4], fill=cbPaletteLt[4]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgsquarewidth, height=imgsquareheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_basecall_2Dmeanbqvs2Dreadlen.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
readlen <- calldata$twod_fastq_seqlen[calldata$twod_fastq_seqlen > 0]
meanbq <- calldata$twod_fastq_bqmean[calldata$twod_fastq_seqlen > 0]
if (length(readlen) == 0 || length(meanbq) == 0) {
    plotdata <- data.frame(readlen=c(0), meanbq=c(0))
    slope <- 0
    intercept <- 0
} else {
    plotdata <- data.frame(readlen=readlen, meanbq=meanbq)
    H <- lm(plotdata$meanbq ~ plotdata$readlen)
    slope <- H$coefficients[2]
    intercept <- H$coefficients[1]
}
plottitle <- sprintf("%s : 2D read length and mean base quality\n(linear regression with 95%% confidence interval, meanbq = %.3f*readlen + %.3f)", rundesc, slope, intercept)
xlabel <- 'Length of 2D read'
ylabel <- 'Mean base quality of 2D read'
p <- ggplot(plotdata, aes(x=readlen, y=meanbq)) +
    geom_point(shape=16, size=0.25, colour=cbPaletteLt[5], fill=cbPaletteLt[5]) +
    geom_smooth(method=lm, colour=cbPaletteDk[4], fill=cbPaletteLt[4]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

# =========================================================================== #
# channel plots                                                               #
# =========================================================================== #

pngpath <- sprintf('%s/%s_channel_readcountcallable.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS = sprintf("%d", min_events)
max_eventsS = formatC(max_events, format="d", big.mark=',')
plottitle <- sprintf("%s : Callable event-reads (length=[%s, %s])", rundesc, min_eventsS, max_eventsS)
xlabel <- 'Channel number'
ylabel <- 'Read count'
B <- data.frame(aggregate(rep(1, length(readdata$channel_number[readdata$event_count >= min_events & readdata$event_count <= max_events])),
    by=list(readdata$channel_number[readdata$event_count >= min_events & readdata$event_count <= max_events]), FUN=sum))
channel_num <- B$Group.1
read_freq <- B$x
p <- ggplot(B, aes(x=channel_num, y=read_freq)) +
    geom_bar(stat="identity", fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    xlim(min(channel_range), max(channel_range)) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_channel_readcounttooshort.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
min_eventsS = sprintf("%d", min_events)
plottitle <- sprintf("%s : Uncallable event-reads (length < %s)", rundesc, min_eventsS)
xlabel <- 'Channel number'
ylabel <- 'Read count'
B <- data.frame(aggregate(rep(1, length(readdata$channel_number[readdata$event_count < min_events])),
    by=list(readdata$channel_number[readdata$event_count < min_events]), FUN=sum))
channel_num <- B$Group.1
read_freq <- B$x
p <- ggplot(B, aes(x=channel_num, y=read_freq)) +
    geom_bar(stat="identity", fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    xlim(min(channel_range), max(channel_range)) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_channel_readcounttoolong.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
max_events <- max(100000, calldata$max_events[1])
max_eventsS = formatC(max_events, format="d", big.mark=',')
plottitle <- sprintf("%s : Uncallable event-reads (length > %s)", rundesc, max_eventsS)
xlabel <- 'Channel number'
ylabel <- 'Read count'
if (sum(readdata$event_count > max_events) > 0) {
    B <- data.frame(aggregate(rep(1, length(readdata$channel_number[readdata$event_count > max_events])),
        by=list(readdata$channel_number[readdata$event_count > max_events]), FUN=sum))
    channel_num <- B$Group.1
    read_freq <- B$x
} else {
    B <- data.frame(channel_num=1:512, read_freq=rep(0, 512))
    channel_num <- channel_range
    read_freq <- rep(0, length(channel_num))
}
p <- ggplot(B, aes(x=channel_num, y=read_freq)) +
    geom_bar(stat="identity", fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    xlim(min(channel_range), max(channel_range)) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_channel_readlenmean.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS = sprintf("%d", min_events)
max_eventsS = formatC(max_events, format="d", big.mark=',')
plottitle <- sprintf("%s : Callable event-read lengths per channel (length=[%s, %s])", rundesc, min_eventsS, max_eventsS)
xlabel <- 'Channel number'
ylabel <- 'Read length (K), with smoothed mean and sd'
C <- data.frame(x=readdata$channel_number[readdata$event_count >= min_events & readdata$event_count <= max_events],
    y=readdata$event_count[readdata$event_count >= min_events & readdata$event_count <= max_events]/kb_number)
p <- ggplot(C, aes(x=x,y=y)) +
    geom_point(size=0.25, colour=cbPaletteLt[3]) +
    geom_smooth(colour=cbPaletteDk[4], fill=cbPaletteLt[4], lwd=0.3) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    xlim(min(channel_range), max(channel_range)) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor = element_line(colour = gridmincol, size=0.25, linetype="solid"))
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_channel_initfilldelay.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
plottitle <- sprintf("%s : Initial channel fill delay", rundesc)
xlabel <- 'Channel number (active only)'
ylabel <- 'Fill delay (min)'
p <- ggplot(filldelay, aes(factor(filldelay$channel_number), filldelay$initfilldelay_seconds/60.0)) +
    geom_bar(stat="identity", fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz*0.7))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major.y = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor.y = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.x = element_blank())
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_channel_refilldelay.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
plottitle <- sprintf("%s : Channel re-fill delay", rundesc)
xlabel <- 'Channel number'
ylabel <- 'Fill delay (min)'
D <- aggregate(x = filldelay$refilldelay_seconds[filldelay$refilldelay_seconds > 0]/60.0,
    by=list(channelnumber = filldelay$channel_number[filldelay$refilldelay_seconds > 0]), FUN=mean, na.rm=TRUE)
channel_number <- D$channelnumber
refill_delay <- D$x
p <- ggplot(D, aes(x=channel_number, y=refill_delay)) +
    geom_bar(stat="identity", fill=cbPaletteDk[3]) +
    theme(legend.position = "none") +
    labs(x = xlabel, y = ylabel) +
    ggtitle(plottitle) +
    theme(axis.text=element_text(size=rel(ticksz*0.7))) +
    theme(axis.title=element_text(size=rel(labelsz))) +
    theme(plot.title=element_text(size=rel(titlesz))) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.x = element_text(angle=90, vjust=1)) +
    theme(plot.background = element_rect(fill = "white")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major.y = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.minor.y = element_line(colour = gridmajcol, size=0.35, linetype="solid")) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.x = element_blank())
ggsave(pngpath, width=imgrectwidth, height=imgrectheight, units="mm", dpi=plotdpi)

pngpath <- sprintf('%s/%s_channel_activitytime.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
min_events <- max(0, calldata$min_events[1])
max_events <- max(100000, calldata$max_events[1])
min_eventsS = sprintf("%d", min_events)
max_eventsS = formatC(max_events/1000.0, format="g", big.mark=',')
plottitle <- sprintf("%s : channel activity over time", rundesc, slope, intercept)
xlabel <- 'Time (h)'
ylabel <- 'Channel number'
png(pngpath, width=12*imgszconst, height=16*imgszconst, units="px")
par(mfrow=c(1,1), oma=c(0,0,0,0))
plot(c(), main=plottitle, xlab=xlabel, ylab=ylabel, xlim=c(0, max(readdata$read_end_hours_since_expt_start)), ylim=c(min(channel_range), max(channel_range)))
max_length <- max(readdata$read_end_hours_since_expt_start)
for (i in channel_range) {
    lines(c(0, max_length), c(i, i), col=cbPaletteLt[8], lwd=1)
}
for(i in 1:length(readdata$fast5_filename)) {
    x0 <- readdata$read_start_hours_since_expt_start[i]
    x1 <- readdata$read_end_hours_since_expt_start[i]
    event_count <- readdata$event_count[i]
    y0 <- readdata$channel_number[i]
    y1 <- y0
    if (event_count < min_events) { colour <- cbPaletteDk[1] } else
    if (event_count <= 1000)      { colour <- cbPaletteDk[2] } else
    if (event_count <= 10000)     { colour <- cbPaletteDk[3] } else
    if (event_count <= 25000)     { colour <- cbPaletteDk[4] } else
    if (event_count <= 50000)     { colour <- cbPaletteDk[5] } else
    if (event_count <= max_events){ colour <- cbPaletteDk[6] } else
                                  { colour <- cbPaletteDk[7] }
    lines(c(x0, x1), c(y0, y1), col=colour, lwd=2)
}
legend(
    max(readdata$read_end_hours_since_expt_start)*0.80,
    max(channel_range),
    pch=15,
    c(sprintf("<%s", min_eventsS), "<=1K", "<=10K", "<=25K", "<=50K",
        sprintf("<=%sK", max_eventsS), sprintf(">%sK", max_eventsS)),
    col=c(cbPaletteDk[1], cbPaletteDk[2], cbPaletteDk[3],
        cbPaletteDk[4], cbPaletteDk[5], cbPaletteDk[6], cbPaletteDk[7]),
    cex=1.5)

# =========================================================================== #
# targets plots                                                               #
# =========================================================================== #

dev.off()
