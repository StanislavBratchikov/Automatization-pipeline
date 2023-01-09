# Automatization-pipeline

This project aims to provide an automatized pipeline for gating cell populations
in fcs files. 
## Intro

Pipeline was inspired by  [this paper](https://pubmed.ncbi.nlm.nih.gov/25378466/)
that has introduced the flowDensity algorithm that 
> facilitates reproducible, high-throughput analysis of flow cytometry data by automating a predefined manual gating approach

It [performed](./Papers%20and%20refs/20220627%20PASC%20BAL%20flow.pdf.pdf) best in `flowCAP-III` competition.


Algorithm is based on estimation of regions surrounding cell populations via
using such characteristics of the marker density distribution as number, height,
width of peaks and the slope of a distribution curve. 
![Alt text](./Papers%20and%20refs/Screen%20Shot%202023-01-09%20at%201.56.02%20PM.png?raw=true "flowDensity perfomance sneak peak")

[AutoPipeline.R](./AutoPipeline.R) provides full gating of 1 given fcs file 
based on [this gating strategy](./Papers%20and%20refs/20220627%20PASC%20BAL%20flow.pdf).
[Auto_pipeline_loop.R](./Auto_pipeline_loop.R) provides
same gating but for a folder of given fcs files.

```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.