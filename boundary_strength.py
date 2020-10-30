
import math

import numpy as np
import scipy
from scipy import spatial
from scipy import stats


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


## ################################  READING GRID                                                                                                                         

## lon, lat, landing vector (indicates the percentage of land)                                                                                                                                                 
lon=np.zeros((num_nodes))
lat=np.zeros((num_nodes))
land_vec=np.zeros((num_nodes))

## reading grid, note that the coordinates on the file refer to the centers of the nodes
i=0                                                                                      
f=open(gridname)
for line in f:
    lat[i] = float(line.split()[0])
    lon[i] = float(line.split()[1])
    land_vec[i] = float(line.split()[2])
    i=i+1


## ################################  CREATING TREE OBJECTS TO FIND NEIGHBOURS OF EACH NODE             

## reshaping lat and lon in a single array
latlon = np.c_[lat,lon]

## creating the tree file from which we can extract neighbours label and distances
tree=scipy.spatial.KDTree(latlon)

## list of 5 nearest neighbours
dist4,neig4=tree.query(latlon, k=5)

## list of 9 nearest neighbours
dist8,neig8=tree.query(latlon, k=9)
    

# ## ################################  DEFINING BOUNDARY VECTORS, STARTING READING LOOP

## vector of boundaies occurrence
cbound_occ = np.zeros((num_nodes,contr_num_sim))
cbound_occ_mean = np.zeros(num_nodes)
cbound_occ_sd = np.zeros(num_nodes)

## vector of boundary strength
cbound_str = np.zeros((num_nodes,contr_num_sim))
cbound_str_mean = np.zeros(num_nodes)
cbound_str_sd = np.zeros(num_nodes)

## starting loop   
ccounter=0                                                                                                                                                                                               
for year in vec_years:
    for daymonth in vec_daymonths:

        ## simulation name file                                                                                                                      
        simul_name = ('simul' + '_startdate' + str(daymonth).zfill(4) + str(year).zfill(4))
                      

        ## ################################  READING MEMBERSHIP                                                                                                                      
                     
        ## name file
        membname= simul_name + '.memb'
        memb = np.zeros(num_nodes)
                      
        f=open(membname)
        for line in f:
            memb[int(line.split()[0])] = int(line.split()[1])

        num_comm = max(memb) + 1

       
	## ################################  READING NETWORK OF COMMUNITIES FILE                                                                                                   

        ## name file                                                                                                               
        netcomname=  simul_name + '.netcomm'

        ## netcomm matrix
        netcomm_mat = np.zeros((num_comm,num_comm))

        ## opening and reading the first line header                                                                                      
        f=open(netcomname)
        f.readline()

        ## reading loop
        for line in f:
            netcomm_mat[int(line.split()[0]), int(line.split()[1])] = float(line.split()[2])


        ## ################################  FINDING ISOLATED NODES                                                                                                                  

        ## vector of isolated nodes
        iso_occ = np.zeros(num_nodes)

        ## loop over all nodes                                                                                                                                     
        for node in range(0,num_nodes):
        ## condition to exclude isolated clusters of nodes                                                                                                                  
            count_similar = 0
            ## loop on neigs asvoiding the node itself                                                                                                                      
            for i in range(1,len(neig8[node,:])):
                if ( memb[neig8[node,i]] == memb[node] ):
                    count_similar = count_similar + 1
            ## isolated node counting                                                                                                                                          
            if (count_similar<3):
                iso_occ[node] = 1


        ## ################################  FINDING BOUNDARIES                                                                                                                         

        ## loop over all nodes
        for node in range(0,num_nodes):
        ## condition to exclude isolated clusters of nodes
            count_similar = 0
            ## loop on neigs avoiding the node itself
            for i in range(1,len(neig8[node,:])):
                if ( memb[neig8[node,i]] == memb[node] ):
                    count_similar = count_similar + 1
            ## condition to avoid isolated nodes           
            if (count_similar>2):
                ## condition to avoid calculation on land
                if (land_vec[neig4[node,0]]<1 and land_vec[neig4[node,1]]<1 and land_vec[neig4[node,2]]<1 and land_vec[neig4[node,3]]<1 and land_vec[neig4[node,4]]<1):
                    ## finding indexes of east, south and south-east nodes
                    node_s = int(neig4[node,][ lat[neig4[node,]] == min(lat[neig4[node,]])][0])
                    node_e = int(neig4[node,][ lon[neig4[node,]] == max(lon[neig4[node,]])][0])
                    node_se = int(neig8[node,][ (lon[neig8[node,]] > lon[node])*(lat[neig8[node,]] < lat[node]) ][0])
                    ## condition of finding a different memb in east, south and south-east node and they are not isolated
                    if ( ((memb[node_s] != memb[node]) and (iso_occ[node_s] == 0)) or 
                         ((memb[node_e] != memb[node]) and (iso_occ[node_e] == 0)) or 
                         ((memb[node_se] != memb[node]) and (iso_occ[node_se] == 0)) ):
                         cbound_occ[node,ccounter] = 1
                        ## calculation of boundary strength using netcomm matrix
                         if ((memb[node_s] != memb[node]) and (iso_occ[node_s] == 0)):
                             p_ij = netcomm_mat[memb[node],memb[node_s]] / np.sum( netcomm_mat[memb[node],:] )
                             p_ji = netcomm_mat[memb[node_s],memb[node]] / np.sum( netcomm_mat[memb[node_s],:] )
                             cbound_str[node,ccounter] = ( 1 - p_ij)*(1 - p_ji)
                         if ((memb[node_e] != memb[node]) and (iso_occ[node_e] == 0)):
                             p_ij = netcomm_mat[memb[node],memb[node_e]] / np.sum( netcomm_mat[memb[node],:] )
                             p_ji = netcomm_mat[memb[node_e],memb[node]] / np.sum( netcomm_mat[memb[node_e],:] )
                             cbound_str[node,ccounter] = ( 1 - p_ij)*(1 - p_ji)
                         if ((memb[node_se] != memb[node]) and (iso_occ[node_se] == 0)):
                             p_ij = netcomm_mat[memb[node],memb[node_se]] / np.sum( netcomm_mat[memb[node],:] )
                             p_ji = netcomm_mat[memb[node_se],memb[node]] / np.sum( netcomm_mat[memb[node_se],:] )
                             cbound_str[node,ccounter] = ( 1 - p_ij)*(1 - p_ji)
        ccounter = ccounter + 1
        print simul_name

## calculating mean and sd                                                                                                                                          
for i in range(0,num_nodes):
    cbound_occ_mean[i] = np.mean(cbound_occ[i,])
    cbound_occ_sd[i] = np.std(cbound_occ[i,])
    cbound_str_mean[i] = np.mean(cbound_str[i,])
    cbound_str_sd[i] = np.std(cbound_str[i,])


