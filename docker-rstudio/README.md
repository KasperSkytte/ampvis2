# RStudio docker container with ampvis2 pre-installed
Use the docker container to run RStudio with a specific version of R and ampvis2 with all required dependencies pre-installed. 

Download and run the `run.sh` script to pull a prebuilt image directly from Docker hub (always inspect scripts before running them):
```
curl -fsSL https://raw.githubusercontent.com/MadsAlbertsen/ampvis2/master/docker-rstudio/run.sh | bash
```

The R version used for the images available on Docker hub will only change with every major R release (i.e. once per year). Only the ampvis2 release will change regularly. If you need a specific R version you can also build the image from scratch if needed by supplying `"build"`, fx:
```
curl -fsSL https://raw.githubusercontent.com/MadsAlbertsen/ampvis2/master/docker-rstudio/run.sh | bash -s "build"
```

A few environment variables can be set before running the script to adjust things: `r_ver`, `ampvis2_rel`, `image_name`, `password`, and `port`. Otherwise default values will be used.

You can then access RStudio from a browser at http://127.0.0.1:${port}.

The `run.sh` script will by default mount the current user's home directory at `/home/rstudio`, and also mount a host [renv](https://rstudio.github.io/renv/index.html) library cache to avoid having to install additional R packages with every container launch. Set `renv::settings$use.cache(FALSE)` in R to disable the cache if problems should occur due to fx mismatching shared libraries between the host and container, or if you just want pure isolation with a project specific R library.
