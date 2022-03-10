#Build your newick tree for bpp with ape package
#Build your tree in reverse, from the bottom to the top 

library(ape)

text.string<-
  "((((((((((((Psammophis_sudanensis, Psammophis_orientalis), Psammophis_subtaeniatus), (Psammophis_afroccidentalis, Psammophis_rukwae)), Psammophis_sibilans), (((Psammophis_philipsii_occidentalis, Psammophis_brevirostris), (Psammophis_mossambicus, (Psammophis_phillipsi, Psammophis_phillipsii)), Psammophis_leopardinus)), Psammophis_praeornatus), Psammophis_lineatus), (Psammophis_biseriatus, Psammophis_tanganicus)), ((Psammophis_punctulatus_trivirgatus, Psammophis_elegans), (Psammophis_aegyptius, Psammophis_schokari))), Psammophis_angolensis), ((((((Psammophis_leightoni, Psammophis_namibensis), Psammophis_trinasalis), Psammophis_notostictus), Psammophis_jallae), Psammophis_trigrammus), Psammophis_lineolatus)), Psammophis_crucifer));"
vert.tree<-read.tree(text=text.string)

plot(vert.tree,no.margin=TRUE,edge.width=2)

is.binary.tree(vert.tree)
