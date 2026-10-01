# Reusable genome plotting functions. Coordinates are 1-based, closed intervals.
# Synthetic chromosomes have no relationship to a reference assembly.
COL <- c(loss = '#2878A0', neutral = '#D9DFE5', gain = '#D65A47', signal = '#267D73')

gallery_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = 'bold', color = '#182B3A'),
      plot.subtitle = ggplot2::element_text(color = '#596773', size = 10),
      strip.text = ggplot2::element_text(face = 'bold'),
      legend.position = 'bottom', plot.margin = ggplot2::margin(10, 15, 10, 10))
}

make_example_data <- function(path = 'data', seed = 42) {
  set.seed(seed)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  chromosomes <- data.frame(chr = paste0('chr', 1:6), length = c(120, 105, 95, 85, 75, 65) * 1e6)
  bins <- do.call(rbind, lapply(seq_len(nrow(chromosomes)), function(i) {
    n <- chromosomes$length[i] / 1e6
    x <- seq_len(n)
    cn <- rep(0, n)
    cn[x >= floor(n * .24) & x <= floor(n * .43)] <- if (i %% 2) .8 else -.8
    cn[x >= floor(n * .68) & x <= floor(n * .8)] <- if (i %% 2) -.6 else .65
    data.frame(chr = chromosomes$chr[i], start = (x - 1) * 1e6 + 1,
      end = x * 1e6,
      signal = pmax(0, 12 + 7 * sin(x / 7 + i) + 14 * (cn > .5) + rnorm(n, 0, 2)))
  }))
  links <- data.frame(chr1 = c('chr1','chr1','chr2','chr2','chr3','chr4','chr5','chr6'),
    start1 = c(32,45,28,75,34,22,53,18)*1e6+1,
    chr2 = c('chr3','chr5','chr4','chr6','chr5','chr6','chr1','chr3'),
    start2 = c(35,24,30,20,28,47,82,68)*1e6+1)
  links$end1 <- links$start1 + 2e6 - 1
  links$end2 <- links$start2 + 2e6 - 1
  genes <- do.call(rbind, lapply(seq_len(nrow(chromosomes)), function(i) {
    n <- 48
    starts <- sort(sample(seq_len(chromosomes$length[i] - 6e5), n))
    data.frame(chr = chromosomes$chr[i], start = starts, end = starts + 5e5,
      gene = paste0('DEMO', i, '_', seq_len(n)), strand = sample(c('+','-'), n, TRUE),
      target = seq_len(n) %in% c(12, 36))
  }))
  bins$gc <- pmin(.7, pmax(.3, .48 + .06*sin(seq_len(nrow(bins))/9) + rnorm(nrow(bins), 0, .025)))
  bins$gene_density <- vapply(seq_len(nrow(bins)), function(i) {
    mid <- (genes$start + genes$end)/2
    sum(genes$chr == bins$chr[i] & mid >= bins$start[i] & mid <= bins$end[i]) /
      ((bins$end[i] - bins$start[i] + 1)/1e6)
  }, numeric(1))
  regions <- data.frame(chr = chromosomes$chr, start = floor(chromosomes$length*.24)+1,
    end = floor(chromosomes$length*.43), label = paste0('R', 1:6))
  for (nm in c('chromosomes','bins','links','genes','regions'))
    write.csv(get(nm), file.path(path, paste0(nm, '.csv')), row.names = FALSE)
  invisible(NULL)
}

read_genome_data <- function(path = 'data') {
  d <- setNames(lapply(c('chromosomes','bins','links','genes','regions'), function(nm)
    read.csv(file.path(path, paste0(nm, '.csv')), stringsAsFactors = FALSE)),
    c('chromosomes','bins','links','genes','regions'))
  validate_genome_data(d)
  d$bins <- d$bins[order(match(d$bins$chr, d$chromosomes$chr), d$bins$start), ]
  d
}

validate_genome_data <- function(d) {
  required <- list(chromosomes = c('chr','length'),
    bins = c('chr','start','end','signal','gc','gene_density'),
    links = c('chr1','start1','end1','chr2','start2','end2'),
    genes = c('chr','start','end','gene','strand','target'),
    regions = c('chr','start','end','label'))
  for (nm in names(required)) {
    if (!is.data.frame(d[[nm]]) || !all(required[[nm]] %in% names(d[[nm]])))
      stop('Missing table or columns: ', nm)
    if (anyNA(d[[nm]][required[[nm]]])) stop('Missing values in ', nm)
  }
  g <- d$chromosomes
  if (!nrow(g) || anyDuplicated(g$chr) || any(!nzchar(g$chr)) ||
      !is.numeric(g$length) || any(!is.finite(g$length)) || any(g$length < 1 | g$length %% 1 != 0))
    stop('Chromosomes require unique names and positive integer lengths.')
  intervals <- function(x, chr = 'chr', start = 'start', end = 'end') {
    a <- x[[start]]; b <- x[[end]]
    if (!is.numeric(a) || !is.numeric(b) || any(!is.finite(a)) || any(!is.finite(b)))
      stop('Coordinates must be finite numbers.')
    if (any(!x[[chr]] %in% g$chr) || any(a < 1 | b < a | a %% 1 != 0 | b %% 1 != 0) ||
        any(b > g$length[match(x[[chr]], g$chr)])) stop('Invalid or out-of-bounds genomic interval.')
  }
  if (!nrow(d$bins)) stop('At least one bin is required.')
  intervals(d$bins); intervals(d$genes); intervals(d$regions)
  intervals(d$links, 'chr1','start1','end1'); intervals(d$links, 'chr2','start2','end2')
  for (nm in c('signal','gc','gene_density'))
    if (!is.numeric(d$bins[[nm]]) || any(!is.finite(d$bins[[nm]]))) stop('Invalid ', nm)
  if (any(d$bins$gc < 0 | d$bins$gc > 1) || any(d$bins$gene_density < 0)) stop('Invalid GC fraction or gene density.')
  if (!is.logical(d$genes$target)) stop('Gene target flags must be TRUE or FALSE.')
  if (any(d$bins$signal < 0)) stop('Signal must be nonnegative.')
  for (chr in unique(d$bins$chr)) {
    x <- d$bins[d$bins$chr == chr, ]; x <- x[order(x$start), ]
    if (nrow(x) > 1 && any(x$start[-1] <= head(x$end, -1))) stop('Bins must not overlap.')
  }
  if (any(!d$genes$strand %in% c('+','-')) || anyDuplicated(d$genes$gene) || any(!nzchar(d$genes$gene)))
    stop('Genes require unique names and strands + or -.')
  invisible(TRUE)
}

plot_circular_genome <- function(d, links_only = FALSE) {
  validate_genome_data(d)
  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  graphics::par(mar = c(4, 1, 3, 1))
  circlize::circos.par(start.degree = 90, gap.degree = 5,
    cell.padding = c(0, 0, 0, 0), track.margin = c(.012, .012))
  g <- d$chromosomes
  circlize::circos.initialize(g$chr, xlim = cbind(0, g$length / 1e6))
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = .075,
    bg.col = '#EDF1F5', bg.border = 'white', panel.fun = function(x, y) {
      chr <- circlize::get.cell.meta.data('sector.index')
      lim <- circlize::get.cell.meta.data('xlim')
      circlize::circos.text(mean(lim), .5, chr, facing = 'bending.inside', cex = .8, col = '#182B3A')
      ticks <- seq(0, max(lim), by = 20)
      circlize::circos.axis(h = 'top', major.at = ticks, labels = ticks,
        labels.cex = .48, major.tick.length = .12)
    })
  if (!links_only) {
    density_max <- max(1, d$bins$gene_density)
    signal_max <- max(1, d$bins$signal)
    for (metric in c('gene_density','gc','signal')) {
      ymax <- switch(metric, gene_density = density_max, gc = 1, signal = signal_max)
      color <- switch(metric, gene_density = '#586B9C', gc = '#DCA957', signal = '#267D73')
      circlize::circos.trackPlotRegion(ylim = c(0, ymax), track.height = .14,
        bg.col = '#F3F5F7', bg.border = NA, panel.fun = function(x, y) {
          b <- d$bins[d$bins$chr == circlize::get.cell.meta.data('sector.index'), ]
          circlize::circos.rect((b$start-1)/1e6, 0, b$end/1e6, b[[metric]], col = color, border = NA)
        })
    }
    graphics::legend('bottom', inset = -.14, xpd = NA, bty = 'n', cex = .7, horiz = TRUE,
      fill = c('#586B9C','#DCA957','#267D73'),
      legend = c(sprintf('Genes/Mb (0–%.0f)', density_max), 'GC (0–1)', sprintf('Signal (0–%.0f)', signal_max)))
  }
  for (i in seq_len(nrow(d$links))) {
    l <- d$links[i, ]
    circlize::circos.link(l$chr1, c(l$start1-1, l$end1)/1e6,
      l$chr2, c(l$start2-1, l$end2)/1e6, col = '#6876AD66', border = NA)
  }
  graphics::title(if (links_only) 'Genomic interval connections (Mb)' else 'Circular genome overview (Mb)',
    cex.main = 1.15, col.main = '#182B3A')
  if (links_only) graphics::mtext('Ribbons connect supplied intervals; width represents interval span',
    side = 1, line = 1, cex = .7, col = '#596773')
}

plot_positions <- function(d) {
  validate_genome_data(d)
  g <- d$chromosomes; r <- d$regions; t <- d$genes[d$genes$target, ]
  for (nm in c('g','r','t')) {
    x <- get(nm); x$y <- match(x$chr, rev(g$chr)); assign(nm, x)
  }
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = g, ggplot2::aes(xmin = 0, xmax = length/1e6, ymin = y-.12, ymax = y+.12),
      fill = '#E6EBEF', color = '#A3AFB9', linewidth = .4) +
    ggplot2::geom_rect(data = r, ggplot2::aes(xmin = (start-1)/1e6, xmax = end/1e6, ymin = y-.12, ymax = y+.12),
      fill = '#DCA957', alpha = .9) +
    ggplot2::geom_text(data = r, ggplot2::aes(x = (start+end-1)/2e6, y = y-.32, label = label),
      color = '#987123', size = 3) +
    ggplot2::geom_point(data = t, ggplot2::aes(x = (start+end-1)/2e6, y = y), color = '#365B80', size = 2) +
    ggrepel::geom_text_repel(data = t, ggplot2::aes(x = (start+end-1)/2e6, y = y, label = gene),
      nudge_y = .32, direction = 'x', size = 3, seed = 42, min.segment.length = 0, max.overlaps = Inf) +
    ggplot2::scale_y_continuous(breaks = seq_len(nrow(g)), labels = rev(g$chr), expand = ggplot2::expansion(add = .6)) +
    ggplot2::labs(title = 'Chromosome positions', subtitle = 'Target gene midpoints (blue) and highlighted intervals (gold)',
      x = 'Position (Mb)', y = NULL) + gallery_theme()
}

plot_linear_genome <- function(d) {
  validate_genome_data(d)
  b <- d$bins; b$chr <- factor(b$chr, levels = d$chromosomes$chr)
  long <- do.call(rbind, lapply(c('gene_density','gc','signal'), function(nm) {
    x <- b; x$value <- x[[nm]]; x$metric <- nm; x
  }))
  long$metric <- factor(long$metric, levels = c('gene_density','gc','signal'),
    labels = c('Genes/Mb','GC fraction','Signal (a.u.)'))
  ggplot2::ggplot(long) +
    ggplot2::geom_rect(ggplot2::aes(xmin = (start-1)/1e6, xmax = end/1e6, ymin = 0, ymax = value, fill = metric)) +
    ggplot2::facet_grid(metric ~ chr, scales = 'free', space = 'free_x') +
    ggplot2::scale_fill_manual(values = c('#586B9C','#DCA957','#267D73'), guide = 'none') +
    ggplot2::scale_x_continuous(expand = c(0,0), breaks = function(lim) seq(0, max(0, max(lim)-1), by = 40)) +
    ggplot2::labs(title = 'Linear genome tracks', subtitle = 'Identical bins and metrics to the circular view',
      x = 'Position within chromosome (Mb)', y = NULL) + gallery_theme() +
    ggplot2::theme(panel.spacing.x = grid::unit(.65, 'lines'))
}

plot_region <- function(d, chromosome = 'chr1', start = 20e6, end = 60e6) {
  validate_genome_data(d)
  if (length(chromosome) != 1 || !chromosome %in% d$chromosomes$chr ||
      length(start) != 1 || length(end) != 1 || !is.finite(start) || !is.finite(end) ||
      start %% 1 != 0 || end %% 1 != 0 || start < 1 || end <= start || end > d$chromosomes$length[match(chromosome, d$chromosomes$chr)])
    stop('Invalid zoom region.')
  b <- d$bins[d$bins$chr == chromosome & d$bins$end >= start & d$bins$start <= end, ]
  if (!nrow(b)) stop('No bins overlap the zoom region.')
  genes <- d$genes[d$genes$chr == chromosome & d$genes$end >= start & d$genes$start <= end, ]
  axis <- ggplot2::scale_x_continuous(limits = c((start-1)/1e6, end/1e6), expand = c(0, 0))
  b$left <- pmax(b$start-1, start-1)/1e6; b$right <- pmin(b$end, end)/1e6
  p1 <- ggplot2::ggplot(b) +
    ggplot2::geom_rect(ggplot2::aes(xmin = left, xmax = right, ymin = 0, ymax = signal), fill = COL['signal']) +
    axis + ggplot2::labs(y = 'Signal (a.u.)', x = NULL) + gallery_theme()
  p2 <- ggplot2::ggplot(b, ggplot2::aes(x = (left+right)/2, y = gc)) +
    ggplot2::scale_y_continuous(limits = c(0,1)) +
    ggplot2::geom_point(color = '#B88735', size = 1.8) + axis +
    ggplot2::labs(y = 'GC fraction', x = NULL) + gallery_theme()
  pd <- ggplot2::ggplot(b) +
    ggplot2::geom_rect(ggplot2::aes(xmin = left, xmax = right, ymin = 0, ymax = gene_density), fill = '#586B9C') +
    axis + ggplot2::labs(y = 'Genes/Mb', x = NULL) + gallery_theme()
  genes$lane <- seq_len(nrow(genes)) %% 3 + 1
  genes$left <- pmax(genes$start-1, start-1)/1e6; genes$right <- pmin(genes$end, end)/1e6
  genes$from <- ifelse(genes$strand == '+', genes$left, genes$right)
  genes$to <- ifelse(genes$strand == '+', genes$right, genes$left)
  p3 <- ggplot2::ggplot(genes) +
    ggplot2::geom_segment(ggplot2::aes(x = from, xend = to, y = lane, yend = lane),
      arrow = grid::arrow(length = grid::unit(2, 'mm')), linewidth = 1, color = '#43546A') +
    ggrepel::geom_text_repel(ggplot2::aes(x = (left+right)/2, y = lane, label = gene),
      nudge_y = .3, direction = 'x', seed = 42, size = 3, box.padding = .25, max.overlaps = Inf, min.segment.length = 0) + axis +
    ggplot2::scale_y_continuous(limits = c(.6, 3.8), breaks = NULL) +
    ggplot2::labs(x = paste0(chromosome, ' position (Mb)'), y = 'Gene spans') + gallery_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  patchwork::wrap_plots(pd, p2, p1, p3, ncol = 1, heights = c(1, 1, 1, 1)) +
    patchwork::plot_annotation(title = paste0('Regional detail / ', chromosome),
      subtitle = 'Aligned gene density, GC, signal and strand-aware gene spans',
      theme = gallery_theme() + ggplot2::theme(plot.margin = ggplot2::margin(22, 15, 10, 30)))
}

# Convert base graphics to vector grid output for patchwork composition.
# A null device avoids creating Rplots.pdf during capture.
circular_grob <- function(d, links_only = FALSE) {
  grDevices::pdf(file = NULL, width = 8, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.grabExpr(grid::grid.draw(gridGraphics::echoGrob(
    function() plot_circular_genome(d, links_only))), width = 8, height = 8)
}

regional_grob <- function(p) {
  grDevices::pdf(file = NULL, width = 11, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  patchwork::patchworkGrob(p)
}
