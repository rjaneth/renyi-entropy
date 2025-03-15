source("../../../Code/R/Programs/read_ENVI_images.R")

# Cargar la imagen SAR
#x <- myread.ENVI(file='../../../Data/SAR/MUNICH_80/Intensity_VV.img',
#                 headerfile='../../../Data/SAR/MUNICH_80/Intensity_VV.hdr')
#x <- myread.ENVI(file='../../../Data/SAR/envi_sf_water_200/Intensity_HH.img', 
#                 headerfile='../../../Data/SAR/envi_sf_water_200/Intensity_HH.hdr')

x <- myread.ENVI(file='../../../Data/SAR/sinaloa_512/Intensity_VV.img', 
                 headerfile='../../../Data/SAR/sinaloa_512/Intensity_VV.hdr')


window_size <- 50

rows <- nrow(x)
cols <- ncol(x)

# Matriz para almacenar los valores de ENL
enl_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Recorrer la imagen con ventanas deslizantes
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Extraer la ventana local
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    
    # Calcular la media y desviación estándar
    mean_value <- mean(as.vector(window_data))
    std_value <- sd(as.vector(window_data))  # Desviación estándar
    
    # Evitar división por cero
    if (mean_value > 0) {
      enl_value <- (mean_value^2) / (std_value^2)  # ENL correcto
      enl_values[i, j] <- enl_value
    }
  }
}

# Calcular el promedio de ENL en toda la imagen
mean_enl <- mean(enl_values, na.rm = TRUE)

# Imprimir el promedio de ENL
print(mean_enl)
