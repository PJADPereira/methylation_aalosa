library(tidyr)
library(car)

options(digits = 10)

## anova_data is provided as an example as it was constructed from other analysis
data <- read.table("./methyl_analysis/anova_data.txt",header=TRUE,stringsAsFactors = FALSE, na.strings = "NA")



AC <- data[which(data$Mutation=="A>G"),c(2:5)]
longAC <- AC %>% pivot_longer(cols = everything(), names_to = "Group",values_to = "Value")
names(longAC) <- c("Group","Freq")
anova_model <- aov(Freq ~ Group, data = longAC)







tukey_results <- TukeyHSD(anova_model,conf.level=0.95)

# View Tukey's HSD results
print(tukey_results)

# Plot Tukey's HSD results
plot(tukey_results)

