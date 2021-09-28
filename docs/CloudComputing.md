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
   1. Make sure the below staging bucket exists;
   2. Make sure the service account (Compute Engine default service account) used by nextflow can write to the bucket at `[Google-storage-bucket]`;
   3. Replace PROJECT_ID in google profile with your Project ID.
   
```angular2html
./nextflow run main.nf\
     -profile ci,docker,google\
     -w [Google-storage-bucket]/nanome-work-test\
     --outputDir [Google-storage-bucket]/nanome-outputs
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
   
    Output will be below.
    ```
    Creating temporary tarball archive of 2 file(s) totalling 4.4 KiB before compression.
    Uploading tarball of [.] to [gs://jax-nanopore-01_cloudbuild/source/1621511968.811599-97a1dc1027de49adba173975e5f98cb3.tgz]
    Created [https://cloudbuild.googleapis.com/v1/projects/jax-nanopore-01/locations/global/builds/13f8d02d-b806-4b8a-b87e-e423508a7373].
    Logs are available at [https://console.cloud.google.com/cloud-build/builds/13f8d02d-b806-4b8a-b87e-e423508a7373?project=619036824958].
    --------------------------------------- REMOTE BUILD OUTPUT ---------------------------------------
    starting build "13f8d02d-b806-4b8a-b87e-e423508a7373"
    ......
    ......
    0ccf7770905c: Pushed
    v1.0: digest: sha256:0a67ba934ad8788bdc2d05ebf4ae30493731f5c5963cd20c259f136bf882dc5f size: 3681
    DONE
    ---------------------------------------------------------------------------------------------------
    
    ID                                    CREATE_TIME                DURATION  SOURCE                                                                                         IMAGES                                 STATUS
    13f8d02d-b806-4b8a-b87e-e423508a7373  2021-05-20T11:59:29+00:00  19M15S    gs://jax-nanopore-01_cloudbuild/source/1621511968.811599-97a1dc1027de49adba173975e5f98cb3.tgz  us.gcr.io/jax-nanopore-01/nanome:v1.0  SUCCESS
    ```
    
Check the Container Regestry link like https://console.cloud.google.com/gcr/images/[PROJECT_ID] for above pushed Docker container.

## Running pipeline

Note that our project id is `jax-nanopore-01`, used for `[PROJECT_ID]`, **Data Bucket** `[Google-storage-bucket]` name used in our project is `gs://jax-nanopore-01-project-data`.

```angular2html
cd nanome
nextflow run main.nf -profile ci,docker,google\
    -w gs://jax-nanopore-01-project-data/nanome-work\
    --outputDir gs://jax-nanopore-01-project-data/nanome-outputs
```


## Troubleshooting
Make sure the network and subnet is 'default' with 'auto' mode
Enable 'Private Google access' for the network/subnet.


## References
* Google cloud for Nextflow: https://cloud.google.com/life-sciences/docs/tutorials/nextflow  
* Nextflow on GCP: https://www.nextflow.io/docs/latest/google.html
* Sandeep sample codes: https://github.com/snamburi3/nextflow-starter-cloud
