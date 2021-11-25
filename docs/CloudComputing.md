# Running pipeline on google cloud platform

## Authentication step (only need once configuration)
1. Install google API from https://anaconda.org/conda-forge/google-cloud-sdk.
   ```angular2html
   conda install -c conda-forge google-cloud-sdk
   ```

2. Authenticate.
   ```angular2html
   gcloud auth login --no-launch-browser
   gcloud auth application-default login
   ```
3. Add Project Name to the config file ~/.config/gcloud/application_default_credentials.json.
   ```angular2html
   "project_id": "[PROJECT-ID]",
   ```

4. Set project.
   ```
   gcloud config set project [PROJECT_ID]
   ```

5. Export environment variables (can be put in ~/.bashrc).
   ```
   export NXF_VER="20.10.0"
   export NXF_MODE=google
   export NXF_DEBUG=3
   export PROJECT="[PROJECT_ID]"
   export GOOGLE_APPLICATION_CREDENTIALS=~/.config/gcloud/application_default_credentials.json
   ```

6. Run GCP nextflow: 
   1. Make sure the below storage bucket exists;
   2. Make sure the service account (Compute Engine default service account) used by nextflow can write to the bucket at `[Google-storage-bucket]`;
   3. Replace PROJECT_ID in google profile with your Project ID.
   
```angular2html
nextflow run TheJacksonLaboratory/nanome\
     -profile test,docker,google\
     -w [Google-storage-bucket]/nanome-work-test\
     --outdir [Google-storage-bucket]/nanome-outputs\
     --googleProjectName  [PROJECT_ID]
```

## Build and submit to container registry of google cloud computing
1. Configure Docker with the following command:
    
    ```angular2html
    gcloud auth configure-docker
    ```
1. Submit to Private Container Registry in the project. You can build the private docker image on Cloud using Dockerfile:
    ```angular2html
    cd nanome
    gcloud builds submit --tag us.gcr.io/[PROJECT_ID]/nanome:latest --timeout=2000s
    ```
Check the Container Regestry link like https://console.cloud.google.com/gcr/images/[PROJECT_ID] for above pushed internal Docker container.

## Running pipeline

Note that our project id is `jax-nanopore-01`, used for `[PROJECT_ID]`, **Data Bucket** `[Google-storage-bucket]` name used in our project is `gs://jax-nanopore-01-project-data`.

```angular2html
nextflow run TheJacksonLaboratory/nanome\
    -profile test,docker,google\
    -w gs://jax-nanopore-01-project-data/nanome-work\
    --outdir gs://jax-nanopore-01-project-data/nanome-outputs\
    --googleProjectName  jax-nanopore-01
```


## Troubleshooting
Make sure the network and subnet is 'default' with 'auto' mode.  
Enable 'Private Google access' for the network/subnet.


## References
* Google cloud for Nextflow: https://cloud.google.com/life-sciences/docs/tutorials/nextflow  
* Nextflow on GCP: https://www.nextflow.io/docs/latest/google.html
* Sandeep sample codes: https://github.com/snamburi3/nextflow-starter-cloud
