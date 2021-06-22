# RStudio docker container with ampvis2 pre-installed
Use the docker container to run RStudio with a specific version of R as well as ampvis2 and all required dependencies pre-installed. 

Run the `buildandrun.sh` script to build from scratch and start RStudio. This can take a while, but you can also just pull a prebuilt image directly from docker hub and run like below instead:

```
r_ver="4.1.0"
ampvis2_rel="2.7.6"
image_name="kasperskytte/rstudio_r${r_ver}_ampvis2:${ampvis2_rel}"
password="supersafepassword"
port="8787"

#pull image from docker hub
docker pull ${image_name}

#launch container
docker run -d --rm \
  -e "PASSWORD=${password}" \
  -v "$HOME":/home/rstudio \
  -p "$port":8787 \
  ${image_name}
```

Then access RStudio from a browser at http://127.0.0.1:${port}. To just use the latest build, set `ampvis2_rel="latest"`.

The `buildandrun.sh` script will by default also mount a host [renv](https://rstudio.github.io/renv/index.html) library cache to avoid having to install additional R packages with every container launch. Set `renv::settings$use.cache(FALSE)` in R to disable the cache if problems should occur due to fx mismatching shared libraries, or if you just want pure isolation with a project specific R library.
