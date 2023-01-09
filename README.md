# Automatization-pipeline

This project aims to provide an automatized pipeline for gating cell populations
in fcs files. 
## Intro

Pipeline was inspired by  [this paper](https://pubmed.ncbi.nlm.nih.gov/25378466/)
that has introduced the flowDensity algorithm that 
> facilitates reproducible, high-throughput analysis of flow cytometry data by automating a predefined manual gating approach

It [performed]('Papers and refs'/'20220627 PASC BAL flow'.pdf) best in `flowCAP-III` competition.


Algorithm is based on estimation of regions surrounding cell populations via
using such characteristics of the marker density distribution as number, height,
width of peaks and the slope of a distribution curve. 
![Alt text]('Papers and refs'/'Screen Shot 2023-01-09 at 1.56.02 PM.png'?raw=true "flowDensity perfomance sneak peak")

[AutoPipeline.R](../AutoPipeline.R) provides full gating of 1 given fcs file 
based on [this gating strategy](../'Papers and refs'/'20220627 PASC BAL flow.pdf').
[Auto_pipeline_loop.R](../Automatization-pipeline/Auto_pipeline_loop.R) provides
same gating but for a folder of given fcs files.

```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.