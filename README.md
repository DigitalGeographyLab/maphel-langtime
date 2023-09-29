### Information

This repo is for Python scripts for *"Capturing urban diversity through languages: long-term changes in multilingual residential neighbourhoods in the Helsinki Metropolitan Area"* published in [Population, Space and Place](https://doi.org/10.1002/psp.2717).

### Pre-analysis steps

* Apply for access to individual-level register data of [Statistics Finland](https://www.stat.fi/tup/tutkijapalvelut/fiona-etakayttojarjestelma_en.html)
* Calculate linguistic diversity of register data per grid cell with [scikit-bio]() similar to how it was done in [a previous article](https://github.com/DigitalGeographyLab/maphel-urbanlingdiv/blob/main/preprocessing/neighborhood_diversities.py) or with language family information in this article's appendices with [calculate_langfam_diversity.py](calculate_langfam_diversity.py).

### Suggested running order of scripts

| Step | Script | Description | Input | Output |
| ---- | :----- | :---------- | :---- | :----- |
| 1 | [add_langfam_genus.py](add_langfam_genus.py) | Adds information on the family and genus per each language | First language information per individual in a CSV file | Language family and genus information per individual |
| 2 | [calculate_langfam_diversity.py](calculate_langfam_diversity.py) | Calculates diversity of language families per grid cell | Output from Step 1 | Geopackages |
| 3 | [calculate_speakers.py](calculate_speakers.py) | Calculates the number of speakers per a grid cell across the language groups used in the article | First language information per individual in a CSV file | Geopackage |
| 4 | [calculate_language_changes.py](calculate_language_changes.py) | Calculates how many times individuals have changed their language | First language information per individual in a CSV file | Pickled dataframe on changes |
| 5 | [HMA_overall_diversityplot.py](HMA_overall_diversityplot.py) | Plots subplots in Figure 1 and full Figure 2 | Individual-level data from Statistics Finland | PDF Figures |
| 6 | [get_gridcell_histories.py](get_gridcell_histories.py) | Gets annual linguistic diversity values per grid cell | Geopackage of diversity information per grid cell | Geopackage of metric-specific grid cell histories |
| 7 | [calculate_livingspace_commutes.py](calculate_livingspace_commutes.py) | Aggregates the records on living space and commute distances to grid cells | FOLK commuting and employment statistics | Pickled data frame |
| 8 | [add_livingspace_commute_to_HMAgrid.py](add_livingspace_commute_to_HMAgrid.py) | Adds the calculated values from Step 5 to the spatial grid data | Step 5 | Geopackage |
| 9 | [plot_diversities_in_new_old_stable_grids.py](plot_diversities_in_new_old_stable_grids.py) | Provides non-normalized plots for Figure 5 | Output from step 4 | PDF Graph and pickled dataframes |
| 10 | [norm_grid_trajectories.py](norm_grid_trajectories.py) | Calculates the normalized values of grid trajectories | Geopackage of diversity information per grid cell | PDF Graphs and pickled dataframes |
| 11 | [popweigh_grid_trajectories.py](popweigh_grid_trajectories.py) | Calculates the population-weighted average values of grid trajectories | Geopackage of diversity information per grid cell | PDF Graphs and pickled dataframes |
| 12 | [plot_grid_diversitiyes_figure5.py](plot_grid_diversitiyes_figure5.py) | Plots final Figure 5 | Outputs from steps 8 and 9 | PDF Graph |
| 13 | [plot_divs_est_som_87-19.py](plot_divs_est_som_87-19.py) | Plots figures 9 and 10, yields outputs | Output from Step 1 | Figures 9 and 10, dataframes |
| 14 | [plot_weighted_divs_est_som_87-19.py](plot_weighted_divs_est_som_87-19.py) | Plots Figure 4 | Output from Step 2| Geopackage with stability classficiations |
| 15 | [calculate_spatial_markov.py](calculate_spatial_markov.py) | Calculate and plot Markov chain matrices | Output from Step 5 | PDF matrices (Figure 6) and matrices as pickled dataframes |
| 16 | [compare_markov.py](compare_markov.py) | Compares Markov probability matrices with Jensen-Shannon distances | Outputs Figure 7 | PNG file |
