# RStudio docker container with ampvis2 pre-installed
Use the docker container to run RStudio with a specific version of R and ampvis2 with all required dependencies pre-installed. 

Download and run the `buildandrun.sh` script to build from scratch and start RStudio. This can take a while, but you can also just pull a prebuilt image directly from docker hub instead by supplying `"pull"`, fx:

```
curl -fsSL https://raw.githubusercontent.com/MadsAlbertsen/ampvis2/master/docker-rstudio/buildandrun.sh | bash -s "pull"
```

A few environment variables can be set before running the script to adjust things: `r_ver`, `ampvis2_rel`, `image_name`, `password`, and `port`. Otherwise default values will be used.

You can then access RStudio from a browser at http://127.0.0.1:${port}.

The `buildandrun.sh` script will by default also mount a host [renv](https://rstudio.github.io/renv/index.html) library cache to avoid having to install additional R packages with every container launch. Set `renv::settings$use.cache(FALSE)` in R to disable the cache if problems should occur due to fx mismatching shared libraries between the host and container, or if you just want pure isolation with a project specific R library.
