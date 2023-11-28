library(FSA)

#nonparametric tests 
indice<-bio_indices$SR
region<-bio_indices$Bioregions
kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")

pal <- c(
  "#002939",
  "#003946",
  "#005f73ff", 
  "#0e9396",
  "#0e9396",
  "#94d2bdff",
  "#94d2bdff",
  "#94d2bdff",
  "#f1f1de",
  "#f1f1de")
  "#e9d8a6",
  "#ebc97d",
  "#ecba53",
  "#ee9b00ff", 
  "#dc8101", 
  "#ca6702",
  "#c35303", 
  "#bb3e03ff",
  "#ae2012ff",
  "#9b2226",
  "#340A00")


bio_indices <- bio_indices %>%
  mutate(code = case_when(
    str_detect(Bioregions, "Turkey_bs1") ~ "TR (West)",
    str_detect(Bioregions, "Turkey_bs2") ~ "TR (East)",
    str_detect(Bioregions, "Turkey_marmara") ~ "TR (Marmara)",
    str_detect(Bioregions, "Ukraine_azov") ~ "UA (Azov)",
    str_detect(Bioregions, "Ukraine") ~ "UA",
    str_detect(Bioregions, "Russia_azov") ~ "RU (Azov)",
    str_detect(Bioregions, "Russia") ~ "RU",
    str_detect(Bioregions, "Georgia") ~ "GE",
    str_detect(Bioregions, "Bulgaria") ~ "BG",
    str_detect(Bioregions, "Romania") ~ "RO"
  ))


library(ggstatsplot)
library(statsExpressions)
library(rlang)
library(glue)

#use this for countries = box plots = 6
pal <- c(
  "#002939",
  "#003946",
  "#005f73ff", 
  "#0e9396",
  "#94d2bdff",
  "#f1f1de")
countries<-ggbetweenstats(data = bio_indices,
               x = Countries,
               y = SR,
               fill = bio_indices$Countries,
               plot.type = "box",
               type = "non-parametric",
               #mean.color = 'black',
              ggtheme = ggpubr::theme_pubclean(),
              boxplot.args = list(color="black", alpha=1, fill=pal, width=0.7),
              pairwise.comparisons = T,
              pairwise.display = "ns",
              p.adjust.method ="bonferroni",
              xlab = "Bioregions",
              ylab= "Species Richness (SR)",
              point.args = list(alpha=0),
              violin.args = list(alpha=0, width=0),
            centrality.type = "parametric",
            centrality.label.args = list(
              size = 3,
              fill="white",
                          label.r=.5,
              nudge_x = 0.2,
              nudge_y = c(2, 2, 0, 0, 0, -2),
              segment.linetype = 4,
              min.segment.length = 0),
            centrality.point.args = list(size = 3, color = "darkred"),
              ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
              title = "",
messages = FALSE)
countries
ggsave("./Maps/countriesboxplot_pubclean.png", countries, width = 10, height = 6.5, dpi = 300, bg = "transparent")



pal <- c(
  "#002939",
  "#003946",
  "#005f73ff", 
  "#0e9396",
  "#0e9396",
  "#94d2bdff",
  "#94d2bdff",
  "#94d2bdff",
  "#f1f1de",
  "#f1f1de")

bioregions<-ggbetweenstats(
  data = bio_indices,
  x = code,
  y = SR,
  fill = bio_indices$Countries,
  plot.type = "box",
  type = "non-parametric",
  #mean.color = 'black',
  ggtheme = ggpubr::theme_pubclean(),
  boxplot.args = list(color = "black", alpha = 1, fill = pal, width = 0.7),
  pairwise.comparisons = T,
  pairwise.display = "ns",
  p.adjust.method = "bonferroni",
  xlab = "Bioregions",
  ylab = "Species Richness (SR)",
  point.args = list(alpha = 0),
  violin.args = list(alpha = 0, width = 0),
  centrality.type = "parametric",
  centrality.label.args = list(
    size = 2.5,
    fill = "white",
    label.r = .5,
    nudge_x = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,-.3, -.3),
    nudge_y = c(0, 2, 0, 0, -1, 0, 2, 0, -2, -2),
    segment.linetype = 4,
    min.segment.length = 0),
  centrality.point.args = list(size = 3, color = "darkred"),
  ggsignif.args = list(textsize = 3,tip_length = 0.01, na.rm = TRUE),
  messages = FALSE)
bioregions

#+theme(text=element_text(family="Jost*", size=14))
#+coord_flip()+theme(text=element_text(family="Jost*", size=14))
ggsave("./Maps/bioregionsboxplot_pubclean2.png", bioregions, width = 10, height = 3.25, dpi = 300, bg = "transparent")

?ggsave

box<-ggplot(bio_indices, aes(x = indice, y = region, fill = Countries)) +
  geom_boxplot(color="black")+
  scale_fill_manual(values=pal)+
  theme_ggstatsplot()+
  labs(title = "Taxonomic Diversity (SR) among Bioregions",
       x = "Species Richness (SR)",
       y = "Bioregions",
       fill = "Bioregions") +
  theme(legend.position = "none")
box

# Create a column for significance level based on the adjusted p-value
posthoc_dunn$res$significance <- ifelse(posthoc_dunn$res$`P.adj` < 0.05, "***", 
                                        ifelse(posthoc_dunn$res$`P.adj` < 0.01, "**", 
                                               ifelse(posthoc_dunn$res$`P.adj` < 0.1, "*", "ns")))

# Create a bar plot for the Z-values with significance level as labels
ggplot(posthoc_dunn$res, aes(x =Z, y =  reorder(Comparison, -Z), fill = significance)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = significance), vjust = -0.5, size = 3.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Z-value") +
  ylab("Comparison") +
  labs(fill = "Significance") +
  scale_fill_brewer()+
  theme_classic2()


indice<-bio_indices$fric
region<-bio_indices$Bioregions
kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")


ggplot(bio_indices, aes(x = region, y = indice, fill = Countries)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Greens")+
  theme_classic2()+
  labs(title = "Functional diversity (fric) among Bioregions",
       x = "Bioregions",
       y = "Functional diversity (fric)",
       fill = "Bioregions") +
  theme(legend.position = "none")

# Create a column for significance level based on the adjusted p-value
posthoc_dunn$res$significance <- ifelse(posthoc_dunn$res$`P.adj` < 0.05, "***", 
                                        ifelse(posthoc_dunn$res$`P.adj` < 0.01, "**", 
                                               ifelse(posthoc_dunn$res$`P.adj` < 0.1, "*", "ns")))

# Create a bar plot for the Z-values with significance level as labels
ggplot(posthoc_dunn$res, aes(x =Z, y =  reorder(Comparison, -Z), fill = significance)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = significance), vjust = -0.5, size = 3.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Z-value") +
  ylab("Comparison") +
  labs(fill = "Significance") +
  scale_fill_brewer(palette="Greens")+
  theme_classic2()



indice<-bio_indices$PD
region<-bio_indices$Countries
kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")


ggplot(bio_indices, aes(x = region, y = indice, fill = Countries)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Reds")+
  theme_classic2()+
  labs(title = "Phylogenetic Diversity (PD) among Countries",
       x = "Bioregions",
       y = "Phylogenetic Diversity (PD)",
       fill = "Bioregions") +
  theme(legend.position = "none")

# Create a column for significance level based on the adjusted p-value
posthoc_dunn$res$significance <- ifelse(posthoc_dunn$res$`P.adj` < 0.05, "***", 
                                        ifelse(posthoc_dunn$res$`P.adj` < 0.01, "**", 
                                               ifelse(posthoc_dunn$res$`P.adj` < 0.1, "*", "ns")))

# Create a bar plot for the Z-values with significance level as labels
ggplot(posthoc_dunn$res, aes(x =Z, y =  reorder(Comparison, -Z), fill = significance)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = significance), vjust = -0.5, size = 3.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Z-value") +
  ylab("Comparison") +
  labs(fill = "Significance") +
  scale_fill_brewer(palette="Reds")+
  theme_classic2()


