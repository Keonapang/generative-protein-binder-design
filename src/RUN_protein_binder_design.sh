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

# ----------------------------------------------------
for cycle in "1A" "1B" "1C" "1D"; do #  "2A" "2B" "2C" "2D"
    python3.11 2_protein_binder_design.py --cycle "$cycle" --num_seq 1 --diffusion 25 --temp 0.3
done

# ----------------------------------------------------
for cycle in "1A" "1B" "1C" "1D" "2A" "2B" "2C" "2D"; do
    python3.11 /home/ubuntu/2_protein_binder_design.py --cycle "$cycle" --num_seq 1 --diffusion 30 --temp 0.35
done
# ----------------------------------------------------

cycle="1" #  cycle="2"
# python3.11 3_protein_binder_design.py --cycle "$cycle" --num_seq 4 --diffusion 25 --temp 0.4
python /home/ubuntu/3_protein_binder_design.py --cycle "$cycle" --num_seq 4 --diffusion 50 --temp 0.5
# python3.11 3_protein_binder_design.py --cycle "$cycle" --num_seq 5 --diffusion 25 --temp 0.4

python 3_protein_binder_design.py --cycle "1" --num_seq 4 --diffusion 50 --temp 0.5

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






docker run --rm --gpus all nvidia/cuda nvidia-smi


export NGC_CLI_API_KEY=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
docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"

mkdir -p ~/.cache/nim
chmod -R 777 ~/.cache/nim    
export HOST_NIM_CACHE=~/.cache/nim

docker compose up
docker run --gpus all -it nvcr.io/nim/deepmind/alphafold2-multimer:2.1
docker container ls


# check health of container
curl localhost:8084/v1/health/ready # alphafold2-multimer
docker run -it --rm --runtime=nvidia --gpus all ubuntu nvidia-smi

# host environment
nvidia-smi --query-gpu=index,name,utilization.gpu,utilization.memory,memory.used,memory.total,temperature.gpu --format=csv
docker exec -it 8ff714cb790d bash

export NGC_CLI_API_KEY=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
export CUDA_VISIBLE_DEVICES=0
export NIM_PARALLEL_THREADS_PER_MSA=12
export NIM_PARALLEL_MSA_RUNNERS=3


#################################################################################################
# Initiate instance
# April 22, 2025

#################################################################################################

# 0. Download the NGC CLI Tool 
wget --content-disposition https://api.ngc.nvidia.com/v2/resources/nvidia/ngc-apps/ngc_cli/versions/3.41.3/files/ngccli_linux.zip -O ~/ngccli_linux.zip && \
unzip ~/ngccli_linux.zip -d ~/ngc && \
chmod u+x ~/ngc/ngc-cli/ngc && \
echo "export PATH=\"\$PATH:~/ngc/ngc-cli\"" >> ~/.bash_profile && source ~/.bash_profile

# 1. Download container
docker run -it --rm --runtime=nvidia --gpus all ubuntu nvidia-smi
docker pull nvcr.io/nim/deepmind/alphafold2-multimer:2.1

# 2. Run the NIM container with the following command.
sudo chmod -R 777 $LOCAL_NIM_CACHE

export LOCAL_NIM_CACHE=~/.cache/nim # path to the cache directory on the host system
docker run -it --rm --name alphafold2-multimer --runtime=nvidia \
    -e CUDA_VISIBLE_DEVICES=0 \
    -e NGC_CLI_API_KEY \
    -v $LOCAL_NIM_CACHE:/opt/nim/.cache \
    -p 8000:8000 \
    nvcr.io/nim/deepmind/alphafold2-multimer:2.1

export LOCAL_NIM_CACHE=~/.cache/nim
docker run -it --rm --name alphafold2-multimer --runtime=nvidia \
    -e CUDA_VISIBLE_DEVICES=0,1 \
    -e NGC_API_KEY=${NGC_CLI_API_KEY:?Error NGC_CLI_API_KEY not set} \
    -e NIM_CACHE_PATH=/home/nvs/.cache/nim/models \
    -e NIM_DISABLE_MODEL_DOWNLOAD=False \
    -e NIM_PARALLEL_THREADS_PER_MSA=12 \
    -e NIM_PARALLEL_MSA_RUNNERS=3 \
    -v $LOCAL_NIM_CACHE:/home/nvs/.cache/nim \
    -p 8080:8000 \
    nvcr.io/nim/deepmind/alphafold2-multimer:2.1

# 3. In a new terminal, check docker container status
curl -X 'GET' \
    'http://localhost:8000/v1/health/ready' \
    -H 'accept: application/json'
docker ps
nvidia-smi --query-gpu=index,name,utilization.gpu,utilization.memory,memory.used,memory.total,temperature.gpu --format=csv

# 4. Run inference to get a predicted protein structure for an amino acid sequence using the following command.

curl -X 'POST' \
    'http://localhost:8000/protein-structure/alphafold2/multimer/predict-structure-from-sequences' \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{"sequences": ["MNVIDIAIAMAI", "IAMNVIDIAAI"]}' > output.json

# Cycle 1A
cycle="1"
diff=50
temp=0.5

curl -X 'POST' \
    -i \
    "http://localhost:8000/protein-structure/alphafold2/multimer/predict-structure-from-sequences"  \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{"sequences": ["KSAKEREARRREELRRESLE", "MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN"], \
        "databases": ["uniref90", "mgnify", "small_bfd"]}' > 4_cycle1A_${diff}diff_${temp}temp.json

# Cycle 1B
curl -X 'POST' \
    -i \
    "http://localhost:8000/protein-structure/alphafold2/multimer/predict-structure-from-sequences"  \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{"sequences": ["TDGQARRQEQQARQQAQEAG", "MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN"], \
        "databases": ["uniref90", "mgnify", "small_bfd"]}' > 4_cycle1B_${diff}diff_${temp}temp.json

# Cycle 1C
curl -X 'POST' \
    -i \
    "http://localhost:8000/protein-structure/alphafold2/multimer/predict-structure-from-sequences"  \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{"sequences": ["EETLRQIEQELAQQRQLREY", "MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN"], \
        "databases": ["uniref90", "mgnify", "small_bfd"]}' > 4_cycle1C_${diff}diff_${temp}temp.json

# Cycle 1D
curl -X 'POST' \
    -i \
    "http://localhost:8000/protein-structure/alphafold2/multimer/predict-structure-from-sequences"  \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{"sequences": ["PLQAAQEQEQIATQIAQQQA", "MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN"], \
        "databases": ["uniref90", "mgnify", "small_bfd"]}' > 4_cycle1D_${diff}diff_${temp}temp.json


import requests
import json
import os

# Define the AlphaFold2-Multimer endpoint
url = "http://localhost:8000/protein-structure/alphafold2/multimer/predict-structure-from-sequences"

# Define binder-target pairs
binder_target_pairs = [
    ["MSDJDJSBSDMAI", "IAMNVIDIAAI"],  # Binder 1 and target
    ["DGFHJKLAI", "IAMNVIDIAAI"],      # Binder 2 and target
    ["HJGKLNBMNVBM", "IAMNVIDIAAI"]    # Binder 3 and target
]

# Headers for the API request
headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

# Create output directory for PDB files
output_dir = "pdb_outputs"
os.makedirs(output_dir, exist_ok=True)

# Variables for tracking results
multimer_results = []  # Stores PDB strings
multimer_response_codes = []  # Stores response codes (e.g., 200 for success)

# Iterate over binder-target pairs
for idx, binder_target_pair in enumerate(binder_target_pairs):
    print(f"Processing pair {idx + 1} of {len(binder_target_pairs)}: {binder_target_pair}")
    
    # Prepare the API payload
    data = {
        "sequences": binder_target_pair,
        "databases": ["uniref90", "mgnify", "small_bfd"]
    }
    
    # Make the POST request
    response = requests.post(url, headers=headers, data=json.dumps(data))
    
    # Check the response
    if response.ok:
        print(f"Request succeeded for pair {idx + 1}")
        multimer_response_codes.append(response.status_code)
        
        # Save the PDB result to a file
        pdb_string = response.text
        pdb_filename = os.path.join(output_dir, f"structure_pair_{idx + 1}.pdb")
        with open(pdb_filename, "w") as pdb_file:
            pdb_file.write(pdb_string)
        
        # Append the result to multimer_results
        multimer_results.append(pdb_string)
    else:
        print(f"Request failed for pair {idx + 1}: {response.status_code}, {response.text}")
        multimer_response_codes.append(response.status_code)
        multimer_results.append(None)

# Print summary
print(f"\nProcessed {len(binder_target_pairs)} binder-target pairs.")
print(f"Response codes: {multimer_response_codes}")

# Save the results to a JSON file for later use
results_file = os.path.join(output_dir, "multimer_results.json")
with open(results_file, "w") as json_file:
    json.dump({
        "binder_target_pairs": binder_target_pairs,
        "multimer_results": multimer_results,
        "response_codes": multimer_response_codes
    }, json_file)

print(f"Results saved to {results_file}")