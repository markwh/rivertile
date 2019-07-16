# rivertile

Functions to work with SWOT rivertile data as produced by [RiverObs](https://github.com/SWOTAlgorithms/RiverObs). Includes uncertainty quantification and validation for SWOT river processing algorithms. 

## Use cases

- Read SWOT river products into convenient tabular format -- **read-products.R**
  - pixel cloud netcdf
  - rivertile netcdf
    - return either nodes or reaches
  - pixcvec netcdf
- Read prior database node and reaches -- **read-priordb.R**
- Transformations of SWOT river products
  - convert to simple features / spatial frames using **sf** package
  - Reprocess aggregation from nodes to reaches (and pixels to nodes?)
- Visualizations
  - maps
  - scatterplots
  - validation--residual distributions

- Validation
  - might belong in a different package--but don't worry about that for now.
  
## Specifications

Functions should reflect specifications in the PDD for these products. 

## Design questions

- Define an S3 class for rivertile data? 
  - I'm currently returning data.frames with netcdf attributes. This may be awkward when processing--attributes may not carry through. 
- Do I try to run RiverObs from within R? 

## TODO

- Add test data
- Add unit tests
- Provide instructions for handling netcdf dependency
- Make default file names more elegant, or remove completely. 
  - default names are only in validation, right? This should probably be a separate package. 

