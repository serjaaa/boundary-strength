
## Overview

The repository `boundary-strength` contains the Python code designed to extract boundaries and their associated strengths from a set of partitions of transport matrices. When using this code, please acknowledge the authors by citing  [Ser-Giacomi et al. (2020)](#references).



## Index
This documentation is organized as follows:

- [Boundary strength code](#boundary-strength-code)
	- [Description](#description)
	- [Inputs](#inputs)
	- [Outputs](#outputs)
	
- [References](#references)



## Boundary strength code

#### Description
The code called `interpolated_boundary_strength.py` extracts boundaries from a set of partitions of transport matrices. First, for each partition all the boundaries are extracted and their strengths calculated. Then an average of such boundary maps is taken across several partitions.  This code has been used to obtain boundary maps of two sets of partitions (associated to historical and scenario climate simulations) and their significant difference.

#### Inputs
The input files are:

- A grid file `.grid` with the latitude, longitude and percentage of land of each node (land_vec). The format for a row is: *(latitude longitude land_vec)* and the line number is the node label.

- A set of membership files `.memb` containing the community to which each node belongs to for each partition considered. The format for a row is: *(node community)*

- A set of network of community files `.netcomm` containing the flows between each pair of communities for each partition considered. The format for a row is:  *(community community flow)*

Each `.memb` and `.netcomm` file name start with `simulation_startdate` and after includes the date associated to the partition in the format *(ddmmyyyy)*.

The input parameters need to be inserted in this first part of the code:
```
## ################################  INPUT PARAMETERS

## grid file name with node positions
gridname = '' + '.grid'

## number of nodes
num_nodes = 

## number of partitions 
contr_num_sim = 

## list of day and month dates
vec_daymonths = []

## list of year dates 
vec_years = range()
```


#### Outputs
The mean and standard deviation of boundary occurrence and strength across partitions are saved in the following part of the code:
```
## calculating mean and sd   
for i in range(0,num_nodes):
    cbound_occ_mean[i] = np.mean(cbound_occ[i,])
    cbound_occ_sd[i] = np.std(cbound_occ[i,])
    cbound_str_mean[i] = np.mean(cbound_str[i,])
    cbound_str_sd[i] = np.std(cbound_str[i,])
```




















## References

[[Ser-Giacomi et al. 2020]](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL089941) Ser-Giacomi, E.,  (2020). Impact of Climate Change on Surface Stirring and Transport in the Mediterranean Sea. *Geophysical Research Letters*, 47, https://doi.org/10.1029/2020GL089941
















~~
