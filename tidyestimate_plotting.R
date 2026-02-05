suppressPackageStartupMessages({
  library(tidyverse)
})

# assumes final_scores exists from the script above
# (or read it back in)
 final_scores <- read.delim("Output_from_tidyestrimate/tidyestimate_scores.tsv", header=T, check.names=F, sep="\t")

# Output paths
out_purity_pdf   <- "/path/tidyestimate_purity_by_MLIVER.pdf"
out_purity_png   <- "/path/tidyestimate_purity_by_MLIVER.png"

out_scatter_pdf  <- "/path/tidyestimate_immune_vs_stromal.pdf"
out_scatter_png  <- "/path/tidyestimate_immune_vs_stromal.png"

out_est_pdf      <- "/path/tidyestimate_ESTIMATE_by_MLIVER.pdf"
out_est_png      <- "/path/tidyestimate_ESTIMATE_by_MLIVER.png"

# 1) Purity by MLIVER
p_purity <- ggplot(final_scores, aes(x = MLIVER, y = purity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.8, size = 1.6) +
  theme_classic() +
  labs(x = "MLIVER", y = "Tumor purity (tidyestimate)", title = "Tumor purity by Liver metastasis status")

ggsave(out_purity_pdf, p_purity, width = 5, height = 4)
ggsave(out_purity_png, p_purity, width = 5, height = 4, dpi = 300)

# 2) Immune vs Stromal scatter
p_scatter <- ggplot(final_scores, aes(x = stromal, y = immune, shape = MLIVER)) +
  geom_point(alpha = 0.85, size = 2) +
  theme_classic() +
  labs(x = "Stromal score", y = "Immune score", title = "Immune vs Stromal (tidyestimate)")

ggsave(out_scatter_pdf, p_scatter, width = 5, height = 4)
ggsave(out_scatter_png, p_scatter, width = 5, height = 4, dpi = 300)

# 3) ESTIMATE score by MLIVER
p_est <- ggplot(final_scores, aes(x = MLIVER, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.8, size = 1.6) +
  theme_classic() +
  labs(x = "MLIVER", y = "ESTIMATE score", title = "ESTIMATE score by Liver metastasis status")

ggsave(out_est_pdf, p_est, width = 5, height = 4)
ggsave(out_est_png, p_est, width = 5, height = 4, dpi = 300)

cat("\nSaved plots:\n",
    out_purity_pdf, "\n", out_purity_png, "\n",
    out_scatter_pdf, "\n", out_scatter_png, "\n",
    out_est_pdf, "\n", out_est_png, "\n")


############## Plotting for Tumor Purity and Estimate Score ###########
library(tidyverse)

infile <- "/path/tidyestimate_scores.tsv"

df <- read.delim(infile, sep = "\t", check.names = FALSE)

# Ensure numeric
df <- df %>%
  mutate(
    estimate = as.numeric(estimate),
    purity   = as.numeric(purity),
    MLIVER   = as.character(MLIVER)
  )

# Color mapping: Y (LiverMet) red, N black
df <- df %>%
  mutate(color = ifelse(MLIVER == "Y", "LiverMet (Y)", "NonLiverMet (N)"))

p <- ggplot(df, aes(x = estimate, y = purity, color = color)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "grey40") +
  theme_classic() +
  labs(
    x = "ESTIMATE score",
    y = "Tumor purity",
    color = ""
  ) +
  scale_color_manual(values = c("NonLiverMet (N)" = "black", "LiverMet (Y)" = "red"))

print(p)

# Save (optional)
ggsave("/path/purity_vs_ESTIMATE_MLIVER_color.pdf",
       p, width = 6, height = 5)
