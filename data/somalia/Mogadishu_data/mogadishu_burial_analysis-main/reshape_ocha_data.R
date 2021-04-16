  ocha <- read_excel("banadir_covid_data.xlsx", sheet = "ocha")
    # remove tibble nonsense
    ocha <- as.data.frame(ocha)

    
a <- melt(ocha, id = c("date", "variable"))    
colnames(a) <- c("date", "variable", "cemetery", "deaths")
write.csv(a, "a.csv", row.names = F)
