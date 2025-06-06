## ----setup------------------------------------------------------------------------------------------------

# Exporta el cuaderno como script plano .R
#knitr::purl("01_perfiles_zonales_SoilGrids_SWAT.qmd")

# Para cargar librerias se verifica pacman
if ("pacman" %in% installed.packages() == FALSE) install.packages("pacman")

# Se cargan las librerias
pacman::p_load(char = c(
  "here", #manejo de rutas
  "sf", #manipulación de dats espaciales
  "dplyr", #procesamiento de data frames
  "tidyr", #procesamient ode data frames
  "remotes", #prepara formato SWAT
  "ggplot2"  #graficación
  )
)

remotes::install_github("biopsichas/SWATprepR")
library(SWATprepR)

# Ajusta tamaño de letra para las gráficas que genere el script
theme(base_size = 14)



## ---------------------------------------------------------------------------------------------------------

#Descarga mapa de ucs armonizadas para Colombia, escala 1:100k desde Zenodo
source(here("Perfiles_SoilGrids", "Scripts", "00_funcion_descarga_ucs_armonizadas_gpkg.R"), encoding = "UTF-8")

#ucs_sf se crea automáticamente al correr la función. Selecciona columnas y filas relevantes
poli_ucs <- ucs_sf  |>
  dplyr::select(UCSuelo, SUBGRUPO, PERFILES, PORCENTAJE, AREA_HA) |>
  sf::st_make_valid()

# Descarga polígono de área de estudio
poli_estudio <- st_read(here("Perfiles_SoilGrids", "Data", "poligono_tocaima.geojson"))

# Transformar ambos al CRS nativo de SoilGrids (IGH: EPSG:54052)
wkt_igh <- '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
poli_ucs_igh <- st_transform(poli_ucs, wkt_igh)
poli_estudio_igh <- st_transform(poli_estudio, st_crs(poli_ucs_igh))


# Clip
estudio_ucs_clip <- st_intersection(poli_ucs_igh, poli_estudio_igh)
estudio_ucs_clip <- estudio_ucs_clip[!st_is_empty(estudio_ucs_clip), ]

#Verifica visualmente
plot(estudio_ucs_clip)



## ----carga------------------------------------------------------------------------------------------------

# Crea un objeto tipo función al ejecutar un  script externo
source("00_funcion_descarga_soilgrids.R")

# Se llama la función con los argumentos adaptados al proyecto
stack_suelo <- descargar_soilgrids_stack(
  vars = c("bdod", "sand", "silt", "clay", "soc"),
  depths = c("0-5cm", "5-15cm", "15-30cm", "60-100cm"),
  stats = c("mean"),
  resolucion = c(250, 250),
  #define ruta de descarga y verifica si ya existen los archivos
  ruta_vrt = here::here("Perfiles_SoilGrids", "Data", "OUTPUT_SoilGrids_vrt") 
)



## ---------------------------------------------------------------------------------------------------------

# Validación rápida: nombres y visual
print(names(stack_suelo))

# Selecciona algunas capas por número de índice
stack_sub <- stack_suelo[[c(1, 8, 15, 22)]]

# Grafica solo esas capas (en un mismo panel multi-cuadro)
plot(stack_sub)



## ---------------------------------------------------------------------------------------------------------

# Definir la ruta de salida
out_raster <- here("Perfiles_SoilGrids", "Data", "OUT_stack_soilgrids.tif")

# Verificar si el archivo ya existe antes de escribir
if (!file.exists(out_raster)) {
  writeRaster(stack_suelo, filename = out_raster, overwrite = TRUE)
} else {
  message("El archivo ya existe, no se sobrescribirá.")
}

# LEE el .tif guardado —esto ya es solo un archivo físico pequeño
stack_suelo_tif <- rast(out_raster)


## ---------------------------------------------------------------------------------------------------------
# (Si es necesario) Proyecta el stack a CRS de los UCS
crs_estudio_ucs <- terra::crs(estudio_ucs_clip)

if (terra::crs(stack_suelo) != crs_estudio_ucs) {
  stack_suelo_tif <- terra::project(stack_suelo_tif, crs_estudio_ucs, method = "bilinear")
} else {
  stack_suelo_tif <- stack_suelo_tif
}

# Transforma UCS a SpatVector, el formato vectorial de Terra
estudio_ucs_vect <- vect(estudio_ucs_clip)

# Recorta al área de interés
ext <- terra::ext(estudio_ucs_vect)
stack_suelo_crop <- terra::crop(stack_suelo_tif, ext)

# Verificación
plot(stack_suelo_crop[[1]])
plot(estudio_ucs_vect, add = TRUE, border = "white", lwd = 2)



## ---------------------------------------------------------------------------------------------------------

# Extracción
tabla_zonal <- terra::extract(stack_suelo_crop, estudio_ucs_vect, fun=mean, na.rm=TRUE)

# tabla_zonal: una fila por polígono, columnas = capas + ID
head(tabla_zonal)




## ---------------------------------------------------------------------------------------------------------

# Crea soil_id único (con sufijo, si hay duplicados)
#tabla_zonal_id <- tabla_zonal |>
 # mutate(soil_id = make.unique(soil_id))

# Pivot a formato largo
tabla_long <- tabla_zonal %>%
  pivot_longer(
    cols = -c(ID),  # Deja fuera ID (solo pivot propiedades)
    names_to = c("property", "depth", "stat"),
    names_pattern = "([a-z]+)_([0-9]+-[0-9]+)cm_(mean|Q0.5|Q0.05|Q0.95|uncertainty)"
  ) %>%
  mutate(
    top = as.numeric(sub("([0-9]+)-([0-9]+)", "\\1", depth)),
    bottom = as.numeric(sub("([0-9]+)-([0-9]+)", "\\2", depth))
  ) %>%
  select(ID, property, stat, value, top, bottom)



## ---------------------------------------------------------------------------------------------------------

# Tabla de conversión: de unidad de SoilGrids a unidad SWAT
# Solo propiedades que interesan a SWAT

tabla_conversion_unidades <- tibble::tribble(
  ~property,     ~unidad_soilgrids, ~unidad_swat, ~factor_conversion, ~nota,
  "sand",        "g/kg",            "%",           0.1,               "Dividir entre 10",
  "silt",        "g/kg",            "%",           0.1,               "Dividir entre 10",
  "clay",        "g/kg",            "%",           0.1,               "Dividir entre 10",
  "soc",         "dg/kg",           "g/kg",        10,                "Multiplicar por 10",
  "bdod",        "cg/cm³",          "g/cm³",       0.01,              "Dividir entre 100",
  "phh2o",       "pH x10",          "pH",          0.1,               "Dividir entre 10",
  "cec",         "mmol(c)/kg",      "cmol(+)/kg",  0.1,               "Dividir entre 10",
  "nitrogen",    "cg/kg",           "g/kg",        0.01,              "Dividir entre 100"
)

tabla_conversion_unidades



## ---------------------------------------------------------------------------------------------------------

tabla_long_swat <- tabla_long|>
  left_join(tabla_conversion_unidades, by = "property") %>%
  mutate(
    value_swat = if_else(!is.na(factor_conversion), value * factor_conversion, value)
  ) |>
  select(ID, property, stat, value_swat, top, bottom)



## ---------------------------------------------------------------------------------------------------------

# Filtra solo las filas que necesitas: stat == "mean"
tabla_wide_swat <- tabla_long_swat %>%
  filter(stat == "mean") %>%
  select(ID, top, bottom, property, value_swat) %>%
  pivot_wider(names_from = property, values_from = value_swat)



## ---------------------------------------------------------------------------------------------------------
# Valida y estructura para SWAT
soil_swat <- get_usersoil_table(
  soil_data    = tabla_wide,
  soil_id      = "ID",
  layer_top    = "top",
  layer_bottom = "bottom",
  sand         = "sand",
  silt         = "silt",
  clay         = "clay",
  soc          = "soc"        
  # Se puede agregar phh2o, bdod, etc., si están disponibles
)

write_usersoil_table(soil_swat, file = "usersoil.txt")


