runs <- c(
  "s1_LE479_NM_T", "LE489KG_Rep", "LE493BR_Rep", "LE497NA_Rep", "LE497ST_Rep", "LE-501-DM_Rep", "LE-517-LD",
  "LE-511-MW", "LE_569_HH", "LE_577_WR", "LE_579_DE", "LE-583-KM", "LE-585-HM", "LE-597-EG_GEX", "LE-595-SV_GEX",
  "LE-593-KP_RNAseq", "LE-589-BE_RNAseq"
)
library(here)
rmarkdown::render(
  "analysis-report-samples.Rmd",
  output_file = here("figs", "my_markdown"),
  params = list(run = runs[1])
)


filenames <- paste(runs, "-report.html")
params <- map(runs, ~list(runs = .))

map2(
  filenames,
  params,
  ~rmarkdown::render("analysis-report-samples.Rmd", output_file = here::here("figs", x.), params = y. )
)
