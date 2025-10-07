## **Terminal (Debian/Ubuntu):**

```
sudo apt-get install -y build-essential libcurl4-openssl-dev libssl-dev libxml2-dev pkg-config
sudo apt-get install -y   build-essential pkg-config   libfreetype6-dev libfontconfig1-dev   libharfbuzz-dev libfribidi-dev   libpng-dev libjpeg-dev libtiff5-dev   libwebp-dev
pkg-config --modversion freetype2 harfbuzz fribidi libpng libjpeg libtiff-4 libwebp # To check if packages install successfully?
sudo apt-get install -y libudunits2-dev libgdal-dev libgeos-dev libproj-dev
sudo apt-get install -y cmake libabsl-dev libprotobuf-dev protobuf-compiler
```

## **R:**
You can run R in terminal by typing `sudo R` or `R`. After finish, you can type `quit()` and select `n` to end the session.
```
install.packages(c("curl","openssl","credentials","urlchecker","rversions"))
install.packages(c("systemfonts","textshaping","ragg"))
install.packages("pkgdown")
install.packages("devtools") # To check if something successfully installed or not, use: packageVersion("ragg"), or whatever packages you want to check
install.packages(c("units", "s2", "sf"))
devtools::install_github("outbreak-info/R-outbreak-info")
```

## **Python**
we also use Python and Snakemake to run the analysis. For the environment setup, we use the envs/workshop environment which was create in the first hand-on tutorial in class. (Link: [here](https://github.com/rki-mf1/2025-SC2-Data-Science/blob/main/day-sc2-seq-and-assembly/hands-on.md))


