# 3D intercellular assembly of decellularized matrix recapitulates in vivo tumor heterogeneity

This code accompanies the paper entitled: "3D intercellular assembly of decellularized matrix recapitulates in vivo tumor heterogeneity." (CCBR-1223)

To reproduce these results, follow these steps on a Linux OS (tested on Oracle Linux Server 7.9):

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/3D-intercellular-assembly-of-decellularized-matrix-recapitulates-in-vivo-tumor-heterogeneity.git```

2.  The input files for this pipeline will be available upon request. Please reach out to the authors before continue to following steps

3.  Install docker and build the docker container:
    * Navigate to the cloned repository directory. 
    * Move to the ./Docker_file/ directory of this repo

4.  Build the container:
    * ```docker build --tag nidap-analysis-r3 .```

5.  Navidate to the cloned repository directory, Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container:
    * ```docker run -ti -v $(pwd)/src:/mnt nidap-analysis-r3```
    
6.  Run the following code.
    * ```cd /mnt```
    * ```bash run_pipeline.sh```

