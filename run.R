# Run from this directory: Rscript run.R or source('run.R').
# Install once:
# install.packages(c('ggplot2', 'patchwork', 'circlize', 'gridGraphics',
#                    'ggrepel', 'ragg', 'svglite'), repos = 'https://cloud.r-project.org')
# Tested: R 4.6.1; ggplot2 4.0.3; patchwork 1.3.2; circlize 0.4.18;
#         ragg 1.5.2; svglite 2.2.2; gridGraphics 0.5.1; ggrepel 0.9.8.
required <- c('ggplot2','patchwork','circlize','gridGraphics','ggrepel','ragg','svglite')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Install missing packages: ', paste(missing, collapse = ', '))
source('plots.R')

regenerate_example <- FALSE # TRUE overwrites the five example CSVs.
zoom_chr <- 'chr1'
zoom_start <- 20e6
zoom_end <- 60e6
image_dpi <- 180
if (regenerate_example) make_example_data()
d <- read_genome_data()

circular <- patchwork::wrap_elements(full = circular_grob(d))
linear <- plot_linear_genome(d)
regional <- plot_region(d, zoom_chr, zoom_start, zoom_end)
positions <- plot_positions(d)
links <- patchwork::wrap_elements(full = circular_grob(d, links_only = TRUE))
combined <- patchwork::wrap_plots(circular, positions, linear,
  patchwork::wrap_elements(full = regional_grob(regional)),
  design = 'AB\nCC\nDD', heights = c(1, .72, 1.05), widths = c(1, 1.1)) +
  patchwork::plot_annotation(title = 'Genome architecture across scales',
    subtitle = 'Synthetic genome / gene density, GC, interval signal and regional detail',
    tag_levels = 'A', theme = gallery_theme())
plots <- list(positions = positions, linear = linear, circular = circular, links = links, regional = regional, combined = combined)
dir.create('figures', showWarnings = FALSE)
for (nm in names(plots)) {
  size <- switch(nm, circular = c(8,8), positions = c(10,6), linear = c(13,6), links = c(8,8), regional = c(11,8), combined = c(15,19))
  for (ext in c('png','svg')) ggplot2::ggsave(file.path('figures', paste0(nm,'.',ext)),
    plots[[nm]], width = size[1], height = size[2], dpi = image_dpi, bg = 'white',
    device = if (ext == 'png') ragg::agg_png else svglite::svglite)
}
message('Done: 5 individual plots + 1 composite, each in PNG and SVG.')
