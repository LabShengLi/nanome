function(input, output, session) {
  library(shiny)
  dataInput <- reactive({
    #dataframe creation
    df %>%
      filter(Dataset %in% input$Dataset, 
             Chrom %in% input$Chrom,
             Strand %in% input$Strand,
             Singleton %in% input$Singleton,
             Genomic_location %in% input$Genomic_location, 
             CpG_location %in% input$CpG_location,
             CG_density %in% input$CG_density,
             Repetitive_regions %in% input$Repetitive_regions)
  })
  output$table1 <- renderDT({
    DT::datatable(dataInput())
  })
  
  
  # Reference: https://github.com/romanhaa/cerebroApp/blob/7fb55f08d1cb9611039e9b73f2859f3a61f4fd42/inst/shiny/v1.3/about/server.R
  output[["about"]] <- renderText({
    paste0(
      '<b>Overview</b>
      <br>
      This application is an online DNA methylation database to display the DNA methylation levels detected by nanopore sequencing
      and bisulfite sequencing data across different genomic contexts.
      <br/>
      Seven methylation calling tools are included in the dataset: Nanopolish, Megalodon, DeepSignal, Guppy, Tombo, METEORE, DeepMod
      <br/>
      Note: 
      <br/>
      1) All start coordinates in the application are 1-based, not 0-based.
      <br/>
      2) Only overlapped sites covered by >= 5 reads in BS-seq and >=3 reads in methylation calling tools are considered
      <br/>
      3) Genomic regions is annotated with promoter > exon > intron > intergenic regions(refer to intergenic) precedence
      <br/>
      4) Promoter is define as the 2000bp up an down of the transcriptional start codon
      <br/>
      <br>
      <b>Funding</b>
      <br>This study was funded and supported by the Jackson Laboratory.<br>
      <br>
      <b>Contact</b><br>
      Ziwei Pan: <a href="mailto:ziwei.pan@jax.org">ziwei.pan@jax.org</a><br>
      <br>
      <b>Citation</b><br>
      If you used the dataset for your research, please cite the following <a href=https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02510-z target="Link to publication">publication</a>:
      <br>
      Liu, Y., Rosikiewicz, W., Pan, Z. et al. 
      DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation. 
      Genome Biol 22, 295 (2021).,
      <a href=https://doi.org/10.1186/s13059-021-02510-z title="DOI number" target="_blank">https://doi.org/10.1186/s13059-021-02510-z.</a><br>
      <br>
      <b>License</b><br>
      The data presented here associated with the publication in the Citation section.
      <br/>
      The application is distributed under the terms of the <a href=https://opensource.org/licenses/MIT title="MIT license" target="_blank">MIT license.</a><br>
      <br>'
    )
  })
  
}
