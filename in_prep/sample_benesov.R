library(tidyverse)
DataBenesov <- read_delim("~/Descargas/LepsetalEcography2019Data/DataBenesov.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

d <- DataBenesov %>% 
  select(-c(Code, Year_fact, `Year_fact*Treatment`,
                 Block, Treatment))
devtools::load_all()

a = clean_community(x = d, time_col = "Year_quant", community_col = "Treatment_Block", filter = FALSE)

b <- split(a, f = a$Treatment_Block)
c <- lapply(b, function(x) {
  logvar_ratio(x = x[,!colnames(x) %in% "Treatment_Block"],
            time_col = "Year_quant", term = "three", log = FALSE)
  }
  )
data.frame(Treatment_Block = names(c),
           t3 = unlist(c)) %>% 
  separate(Treatment_Block, sep ="_", into = c("Block", "Treatment")) %>% 
  ggplot(aes(x = Treatment, y = t3-1))+
  geom_boxplot()
