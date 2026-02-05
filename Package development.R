install.packages(c("usethis", "devtools"))
usethis::create_package("spatialPipelineR")
usethis::use_git()
usethis::use_readme_md()

usethis::use_git_config(
  user.name = "Chi Zhang",
  user.email = "zhanc183@mcmaster.ca"
)
usethis::use_git()
library(spatialPipelineR)
devtools::install()
ls("package:spatialPipelineR")
usethis::use_package("Seurat")
devtools::document()
devtools::install()
library(spatialPipelineR)

liver2 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/ExciseProject/Sam28",
  slice = "condition2a"
)

devtools::install_github("ChiZhangMcMaster/liverST")
library(spatialPipelineR)


remove.packages("spatialPipelineR")
library(spatialPipelineR)
find.package("spatialPipelineR")

devtools::document()
# Remove old package if it exists
remove.packages("spatialPipelineR")

# Install fresh from local source
devtools::install()
library(liverST)
ls("package:liverST")  # Should show all your functions
usethis::use_git()
# Commit message example:
# "Rename package to liverST"
usethis::use_github()  # if not already linked

