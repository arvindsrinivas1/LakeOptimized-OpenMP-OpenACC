NOTE: This repo has been moved from my personal university repository to this personal public repository


# LakeOptimized-OpenMP-OpenACC
Optimized a serial lake simulation 10X faster using OpenMP, OpenACC.

The lake simulation logic generates the state of the lake after pebbles are dropped on it. This is written using Centralized Finite Difference. 
I use OpenMP and OpenACC here to run it 10X faster. 

It was run on a Opetron Core with 16 threads allocated:
  srun -n16 -popteron --pty /bin/bash

