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
export NGC_CLI_API_KEY=<enter-key> #  export NGC_CLI_API_KEY=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
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
curl localhost:8084/v1/health/ready # Protein MPNN

pip install requests

# EXECUTION
cycle = 1 # cycle 1 or cycle 2
num_seq = 4 # number of sequences to generate per target
diffusion = 30 # number of diffusion steps (15-30 recommended)
temp = 0.2 # sampling temperature (range: 0-1) to adjust the probability values for the 20 amino acids at each position, controls the diversity of the design outcomes

for cycle in "1" "2"; do
    python3.11 1_protein_binder_design.py --cycle "$cycle" --num_seq 5 --diffusion 25 --temp 0.4
done

for cycle in "1A" "1B" "1C" "1D"; do #  "2A" "2B" "2C" "2D"
    python3.11 2_protein_binder_design.py --cycle "$cycle" --num_seq 1 --diffusion 25 --temp 0.3
done

# ----------------------------------------------------
for cycle in "1A" "1B" "1C" "1D"; do
    python3.11 2_protein_binder_design.py --cycle "$cycle" --num_seq 2 --diffusion 25 --temp 0.3
done

python3.11 2_protein_binder_design.py --cycle "1" --num_seq 2 --diffusion 25 --temp 0.2
python3.11 2_protein_binder_design.py --cycle "2" --num_seq 2 --diffusion 25 --temp 0.1
python3.11 2_protein_binder_design.py --cycle "2" --num_seq 2 --diffusion 25 --temp 0.2

# ----------------------------------------------------

python3.11 -m pip install biopython
python3.11 -m pip install prodigy-prot


# ----------------------------------------------------
# Convert JSON to PDB
# ----------------------------------------------------

json="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/EVO2/OUTPUT_RF/2_cycle1A_A400_440.json"
Rscript "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Keona_scripts/generative-protein-binder-design/src/convert_json_to_pdb.r" $json


