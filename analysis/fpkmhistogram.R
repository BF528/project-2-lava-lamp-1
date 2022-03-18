library(tidyverse)

path = '/projectnb/bf528/users/lava-lamp/project-2-lava-lamp-1/samples/P0_1_tophat/P0_1_cufflinks/genes.fpkm_tracking'
values <- read.table(path, header = TRUE)

filtered <- values %>% filter(FPKM >= 1 ) %>% mutate(logFPKM = log10(FPKM)) 


ggplot(filtered, aes(logFPKM)) +
  geom_histogram(bins = 13, fill = 'Grey', color = 'Black') +
  labs(title = "Distribution of Log10 FPKM Values for FPKM > 1", y="Count", x="Log10 FPKM")




