#logo_source.R
library(magick)
source('invert_meth_ge_functions.R')

dat_names

logo_list = list()
print('loading logo images to logo_list...')
for (n in dat_names){
  print(n)
  logo_file = paste('metadata/species_logos/', n, '.png', sep='')
  logo_list[[n]] = magick::image_read(logo_file)
}

