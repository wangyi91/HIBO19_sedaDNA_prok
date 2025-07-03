library(decontam)
library(tidyverse)
library(scico)

# metadata
mt = read.csv(file.path("../metadata/HIBO_library_metadata.csv"))

input_dir = "../deContamination/output"
rank="species"

taxa_names <- read.csv(file.path(input_dir, paste("ra_matrix.", rank, ".csv", sep = "")), header=F, nrows = 1) %>% 
    t %>% as.data.frame %>% filter(V1!="Label") %>% `colnames<-`("taxa_name")

taxa <- read.csv(file.path(input_dir, paste("ra_matrix.", rank, ".csv", sep = ""))) 
taxa_matrix <- taxa %>% select(-Label) %>% as.matrix

dmg_df <- read.csv(file.path(input_dir, paste("dmg_matrix.", rank, ".csv", sep = "")))
dmg_matrix <- dmg_df %>% select(-Label) %>% as.matrix

sample_data <- read.csv(file.path(input_dir, paste("sample_data.", rank, ".csv", sep = "")))
sample_data$iscontrol <- sample_data$sample_type != "sample"

method = "frequency"# "either",# "frequency", "prevalence",
method = "prevalence"
method = "frequency"# "either",# "frequency", "prevalence",
contamdf <- isContaminant(taxa_matrix, 
                               conc = sample_data$molarity, 
                               neg = sample_data$iscontrol, 
                               method = method, 
                               threshold=0.1)
table(contamdf$contaminant)

tax_unchanged = cbind(contamdf, taxa_names) %>% filter(contaminant==TRUE) %>% select(taxa_name) %>% as.list
tax = contamdf %>% filter(contaminant==TRUE) %>% rownames

tax_dict <- setNames(tax_unchanged$taxa_name, tax)

# Export taxa list for downstream processing in julia
#write.csv(tax_unchanged$taxa_name, "output/contam_taxa_name.csv", row.names = F, quote = TRUE)
writeLines(tax_unchanged$taxa_name, "output/contam_taxa_name.txt")



# Plotting:

# Flatten the data
all_values <- as.vector(as.matrix(dmg_matrix))
subset_values <- as.vector(as.matrix(dmg_matrix[, colnames(dmg_matrix) %in% tax]))

# Combine into a data frame for plotting
plot_data <- data.frame(
  value = c(all_values, subset_values),
  group = c(rep("All taxa", length(all_values)),
            rep("Contaminants", length(subset_values))
           )) %>% na.omit

# Plot overlaid histograms
ggplot(plot_data, aes(x = value, fill = group)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity", bins = 50) +
  scale_fill_manual(values = c("All taxa" = "darkgrey", "Contaminants" = "darkblue")) +
  labs(title = paste('"',method,'" method', sep=""),
       x = "DNA damage",
       y = "Density") +
  theme_linedraw()

ggsave(paste("output/density_plot_", method, ".pdf", sep=""))


# plot dmg against age
dmg_contam <- dmg_df %>% select(Label,all_of(tax)) %>% left_join(select(mt,yrBP,Label), by="Label")
ra_contam <- taxa %>% select(Label,all_of(tax))
  
dmg_contam_long <- pivot_longer(dmg_contam, cols = tax, names_to = "tax_name", values_to = "damage")# %>% na.omit()
ra_contam_long <- pivot_longer(ra_contam, cols = tax, names_to = "tax_name", values_to = "relative_abd")

dt <- left_join(dmg_contam_long, ra_contam_long,by=c("Label","tax_name"))

dt$tax_name <- tax_dict[dt$tax_name]

ggplot(dt, aes(x = damage, 
               y = factor(Label, levels = dmg_contam$Label), 
               color = tax_name, 
               size = relative_abd)) +
  geom_point() +
  guides(color = "none") +
  scale_size_continuous(range = c(0, 5), trans = "sqrt", 
                        breaks = c(0.001, 0.01, 0.05, 0.2, 0.5), 
                        name = "Relative abundance") +
  labs(x = "DNA damage", y = "", title = "") +
  #scale_color_scico_d(palette = "batlow", name = "Taxa name") +
  theme_bw() + 
  theme(#legend.text = element_text(face = "italic"),
        legend.position = "right")  # or "bottom" / "top" / "none"

ggsave("output/dotplot_alltaxa_inctrl.pdf", width = 6, height = 6, units = "in")
