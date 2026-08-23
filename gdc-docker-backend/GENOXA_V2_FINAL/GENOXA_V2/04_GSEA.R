# Define selectInput options for organism database
organismDbOptions <- selectInput("organismDb", "Select Organism Database:",
                                 choices = c(
                                   "Human (org.Hs.eg.db)" = "org.Hs.eg.db", 
                                   "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db"
                                 ), selected = "org.Hs.eg.db"
                          )

# Define selectInput options for key type
keytypeOptions <- selectInput("keytype", "Select Key Type:",
                              choices = c(
                                "ACCNUM" = "ACCNUM",
                                "ALIAS" = "ALIAS",
                                "ENSEMBL" = "ENSEMBL",
                                "ENSEMBLPROT" = "ENSEMBLPROT",
                                "ENSEMBLTRANS" = "ENSEMBLTRANS",
                                "ENTREZID" = "ENTREZID",
                                "GENENAME" = "GENENAME",
                                "REFSEQ" = "REFSEQ",
                                "SYMBOL" = "SYMBOL",
                                "UNIGENE" = "UNIGENE",
                                "UNIPROT" = "UNIPROT"
                              ),
                              selected = "ENSEMBL"
                        )
