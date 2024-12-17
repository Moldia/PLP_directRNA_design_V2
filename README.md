# Installation  
## Step 1: Install Docker
Before using the provided Dockerfile, ensure Docker is installed on your system. Follow the official instructions based on your operating system:

[Docker Installation Guide](https://docs.docker.com/engine/install/)  
Ensure Docker is properly installed by running the following command in your terminal:  
```
docker --version
```
## Step 2: Clone the GitHub Repository
Clone the repository containing the Dockerfile to your local machine. Replace <repository-url> with the actual URL of your repository:
```
git clone https://github.com/Moldia/PLP_directRNA_design_V2.git
cd PLP_directRNA_design_V2
```
## Step 3: Build the Docker Image
Build the Docker image using the docker build command. Replace <image-name> with a name of your choice for the image:
```
docker build -t plp_probe_design_v2 .
```
# Running the docker image   
## Step 1: Run the Docker Container
After successfully building the image, create and run a container using the following command:
```
docker run -it --rm -p 2222:2222 plp_probe_design_v2
```
## Step 2: Checking installation of required tools (`blat` & `cutadapt`)
You can test whether `blat` is installed just running the tool without any parameters which print the help:  
```
blat 
```
You can test whether `cutadapt` is installed:  
```
cutadapt --version
```
## Step 3: Start Jupyter Notebook from Within the Container
After entering the container's shell, start the Jupyter Notebook with the following command:
```
jupyter notebook --ip=0.0.0.0 --port=2222 --no-browser --allow-root
```

# Mounting data to Docker Container
You should mount the data and relevant files (e.g. jupyter notebooks)
```
docker run -it -p 2222:2222 \
  -v /home/nima/Lee_2023/nima_dataset/:/app \
  -v /home/nima/PLP_directRNA_design_V2/codes/Benchmarking.ipynb:/app/Benchmarking.ipynb  \
  plp_probe_design_v2
```
Any changes you make to the mounted files or directories inside the container will be reflected in the original files on your host system. 

Then in your browser open this address:
`http://localhost:2222` 
