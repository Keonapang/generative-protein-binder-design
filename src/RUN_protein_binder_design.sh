# Keona Pang
# Apr 18, 2025

# Hardware requirements: 4 x A100 CRUSOE instance
# 4 x GPU, 47 GiB GPU memory
# 128 GB SSD drive space, 60GiB RAM,24 CPU

#INPUT FILES TO UPLOAD TO JUPYTER NOTEBOOK
#   - `protein-binder-design_v3.ipynb` uploaded to a [public Github repo](https://github.com/Keonapang/generative-protein-binder-design/blob/main/src/protein-binder-design.ipynb)
#  - `docker-compose.yaml` (3MB) from the original [BioNeMo repo](https://github.com/NVIDIA-BioNeMo-blueprints/generative-protein-binder-design/blob/main/deploy/docker-compose.yaml)
#  - `cycle1_alphafold2_output.pdb` (80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)

# Get from brev.nvidia
brew install brevdev/homebrew-/brev && brev login --token <****> # Install the CLI
brev shell <instance-name> # find instance-name on brev.nvidias # Open a terminal locally

# Run on Terminal
sudo apt-get install -y docker-compose
sudo apt install python3.11

# Log in to Docker
export NGC_CLI_API_KEY=<enter-key> # enter personal API key
docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"

## Create the nim cache directory to download any model data to your local/server disk 
mkdir -p ~/.cache/nim
chmod -R 777 ~/.cache/nim    ## Make it writable by the NIM
export HOST_NIM_CACHE=~/.cache/nim

# Run Docker compose
docker compose 

# check status of NIMs
curl localhost:8082/v1/health/ready # RFdiffusion
curl localhost:8083/v1/health/ready # Protein MPNN

pip install requests

# EXECUTION
cycle = 1 # cycle 1 or cycle 2
num_seq = 4 # number of sequences to generate per target
diffusion = 30 # number of diffusion steps (15-30 recommended)
temp = 0.2 # sampling temperature (range: 0-1) to adjust the probability values for the 20 amino acids at each position, controls the diversity of the design outcomes
python SCRIPT_protein_binder_design.py --cycle "1" --num_seq 4 --diffusion 30 --temp 0.2

