#!/usr/bin/env python
# -*- encoding: utf8

#  Autores: Yanina Bellini Saibene

# Genera las imagenes con la clasificación de regresion logística
#TODO: tomar los archivos a procesar por parámetro

from identify import Identify

try:
	from osgeo import gdal
	from osgeo import osr 
except ImportError:
	raise ImportError,"Se requiere el modulo osgeo.  Se puede descargar de http://www.lfd.uci.edu/~gohlke/pythonlibs/"

try:
	import numpy 
except ImportError:
	raise ImportError,"Se requiere el modulo numpy. http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy"

import pyodbc
from numpy import e
from math import *

print "Inicio..."

def gepModelHail(d):

	G2C3 = 6.64830583849605
	G3C2 = 9.20590838343455
	G4C5 = -7.728054750206
	G4C9 = -2.29385333017975
	G4C7 = 5.17282326731162

	MxdBZ = 0
	TotdBZ = 2
	AvdBZ = 3
	MnRho = 5
	TotRho = 6
	AvRho = 7
	MxZDR = 8
	MnZDR = 9
	AvZDR = 11

	y = 0.0

	y = (d[MxdBZ]/d[AvZDR])
	y = y + ((((gepAND1(d[TotdBZ],d[MxZDR])*d[TotRho])-(d[AvZDR]*G2C3))-d[AvZDR])+(d[TotRho]-d[MnZDR]))
	y = y + gepLT2E((pow((d[MxdBZ]-(d[MxZDR]-G3C2)),2.0)+(((d[MxZDR]-d[AvRho])+d[AvdBZ])/2.0)),d[MxdBZ])
	y = y + gepLT2E(gepLT2G(d[MxZDR],((((G4C5-G4C9)*(G4C7*d[MnRho]))+gepLT2A(d[TotdBZ],d[MnZDR]))/2.0)),d[MxdBZ])

	SLOPE = 2.46515518080859E-05
	INTERCEPT = -2.14079328432395

	probabilityOne = 1.0 / (1.0 + exp(-(SLOPE * y + INTERCEPT)))
	return probabilityOne

def gepAND1(x, y):
	if ((x < 0.0) and (y < 0.0)):
		return 1.0
	else:
		return 0.0

def gepLT2A(x, y):
	if (x < y):
		return x
	else:
		return y

def gepLT2E(x, y):
	if (x < y):
		return (x+y)
	else:
		return (x*y)

def gepLT2G(x, y):
	if (x < y):
		return (x+y)
	else:
		return atan(x*y)


print "Inicio..."
#Abro los archivos de plantilla
imagenGranizoRL1 = gdal.Open("20121210MAXdBZ.vol_1..tif")

#Obtengo la primera banda
bandaimagenGranizoRL1 = imagenGranizoRL1.GetRasterBand(1)


#Los paso a una matriz
# Esto es la matriz NxM de datos de la imagen, acá están TODOS los 
# valores de pixeles de la imagen
datosimagenGranizoRL1 = bandaimagenGranizoRL1.ReadAsArray(0, 0, imagenGranizoRL1.RasterXSize, imagenGranizoRL1.RasterYSize)


d=[0,0,0,0,0,0,0,0,0,0,0,0,0]
#Imgs= ma.zeros(12,505,487)
Imgs=[0,0,0,0,0,0,0,0,0,0,0,0,0]
#Armo la matriz con las imagenes de totales necesarias
img_tif_name="20121210MAXdBZ.vol_1..tif"
Imgs[0] = Identify(img_tif_name)
img_tif_name="20121210TOTdBZ.vol_1..tif"
Imgs[2] = Identify(img_tif_name)
img_tif_name="20121210AVGdBZ.vol_1..tif"
Imgs[3] = Identify(img_tif_name)
img_tif_name="20121210MINRhoHV.vol_1..tif"
Imgs[5] = Identify(img_tif_name)
img_tif_name="20121210TOTRhoHV.vol_1..tif"
Imgs[6] = Identify(img_tif_name)
img_tif_name="20121210AVGRhoHV.vol_1..tif"
Imgs[7] = Identify(img_tif_name)
img_tif_name="20121210MAXZDR.vol_1..tif"
Imgs[8] = Identify(img_tif_name)
img_tif_name="20121210MINZDR.vol_1..tif"
Imgs[9] = Identify(img_tif_name)
img_tif_name="20121210AVGZDR.vol_1..tif"
Imgs[11] = Identify(img_tif_name)

print 'Iniciando obtención de datos para la fecha determinada.'
for y in xrange(504): #Todas las imágenes tienen el mismo tamaño
	for x in xrange(486): #Todas las imágenes tienen el mismo tamaño
		d[0]=Imgs[0].get_pixel_data(y,x)
		d[2]=Imgs[2].get_pixel_data(y,x)
		d[3]=Imgs[3].get_pixel_data(y,x)
		d[5]= Imgs[5].get_pixel_data(y,x)
		d[6]=Imgs[6].get_pixel_data(y,x)
		d[7]=Imgs[7].get_pixel_data(y,x)
		d[8]=Imgs[8].get_pixel_data(y,x)
		d[9]=Imgs[9].get_pixel_data(y,x)
		d[11]=Imgs[11].get_pixel_data(y,x)
	
		print "d:", d
		granizo=gepModelHail(d)
		print "granizo:", granizo
		datosimagenGranizoRL1[y][x]=granizo
	
	

#Guardo las  imagenes, cuando sale del for.
#Variables para definir la proyección
new_ds_ref = osr.SpatialReference() 
new_ds_ref.ImportFromEPSG(4326)  
driver = gdal.GetDriverByName("GTiff")

print 'Grabo Granizo RL1'
ds = driver.Create("GranizoSiNo1EleCon99ImpImg.tif", imagenGranizoRL1.RasterXSize, imagenGranizoRL1.RasterYSize, 1, gdal.GDT_Float64)
ds.SetGeoTransform(imagenGranizoRL1.GetGeoTransform()) 
ds.SetProjection(new_ds_ref.ExportToWkt())
outband=ds.GetRasterBand(1)
outband.SetNoDataValue(numpy.nan)
outband.WriteArray(datosimagenGranizoRL1)
ds = None 


print "Almacenamiento finalizado"
#con.close()

