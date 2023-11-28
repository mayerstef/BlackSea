# Create decision tree
graph <- "
digraph structs {
  node [shape=box, style=filled, fillcolor=lightblue];
  edge [color=gray50];
  
  struct1 [label='All Species\n Generate raster maps from available data sources'];
  
  struct2 [label='Review each map for sensibility and accuracy \n Is one map clearly more accurate than the others?'];
  
  struct3 [label='Use the most accurate map', fillcolor=forestgreen];
  
  struct4 [label='Are there maps that could be combined?'];
  
  struct5 [label='Combine maps to create more accurate representation', fillcolor=forestgreen];
  
  struct6 [label='Use the maps available, acknowledging their limitations', fillcolor=lightcoral];
  
  struct1 -> struct2;
  
  struct2:s -> struct3:n [label='Yes', fontsize=10];
  
  struct2:s -> struct4:n [label='No', fontsize=10];
  
  struct4:s -> struct5:n [label='Yes', fontsize=10];
  
  struct4:s -> struct6:n [label='No', fontsize=10];
}
"

# Render the graph
DiagrammeR::grViz(graph)

