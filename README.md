# phylodyn2

This package originated as a copy of [phylodyn](https://github.com/mdkarcher/phylodyn), an R package with the purpose of facilitating phylodynamic inference and analyses.
`phylodyn2` has reduced, specific functionality compared to phylodyn, and it implements methodology developed for accounting for reporting delays, by incorporating reporting probabilities in the preferential sampling model for real-time phylodynamic analyses.

The peer-reviewed manuscript describing the methodolody is available at <https://doi.org/10.1371/journal.pcbi.1012970> and the repository with code to reproduce the manuscript is available at <https://github.com/CatalinaMedina/reporting-delays-in-phylodynamics-paper>.

## Installation

```{r}
# install.packages("devtools")
devtools::install_github("CatalinaMedina/phylodyn2")
```
