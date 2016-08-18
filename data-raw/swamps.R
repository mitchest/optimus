library(tidyr)
library(dplyr)

# access .csv from:
# http://www.ltern.org.au/knb/metacat?action=read&qformat=html&docid=ltern.84.15
# but it sits behind a password, so it's in /data-raw
ltern_long <- read.csv("data-raw/ltern.85.2-kuhs_vegetation_data_floristics_2014_p110t165.csv", stringsAsFactors = F)

ltern_long %>%
  select(-comment) %>%
  spread(species, frequency) -> swamps

devtools::use_data(swamps, overwrite = TRUE)
