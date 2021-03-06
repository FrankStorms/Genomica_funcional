setwd("6to semestre/Genomica funcional/Archivo/")

library("sleuth") #Cargamos la libreria 

#Esta funci�n permite mapear, a partir de la base de datos de
tx2gene <- function(){
#Con la funci�n de abajo metida al objeto mart, indicamos a donde nos queremos conectar, y a
  # que base de datos queremos ingresar, para nuestro caso sera a la base de datos de homo sapiens
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#ahora debemos indicar que es lo queremos extraer de esa base de datos, en este caso nos quedamos
 # con el identificador del transcripto, el identificador del gen y el nombre externo 
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
# ahora renombramos los conceptos, de tal forma que cuando se genere el objeto, venga los nombres
  # de las columnas m�s cortos
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g) # para que despues de que se ejecute la funci�n, nos regrese la tabla que se genero
}

#De esta manera, lo que hacemos es que la funci�n tenga una sintaxis a modo de objeto
t2g <- tx2gene()


#De esta manera nos metemos a la carpeta en donde tenemos nuestros archivos
base_dir<-"~/6to semestre/Genomica funcional/Archivo/"

#hacemos un objeto, en donde se incluyen los nombres de las carpetar que vamos a utilizar
samples <- paste0("sample", c("1","2","3",
                                   "4"))

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id)) #esto no lo entendi, aunque 
#entiendo que tiene que ver con la manera en que se leeran los archivos que tengo en mis carpetas

#De esta forma selecciona dos grupos experimentales, dos controles y dos L1, que corresponden 
# al orden en que se pusieron las carpetas anteriores, los samples.
s2c <- data.frame(path=kal_dirs, sample=samples, muestras = c("CTRL","CTRL","L1", "L1"
                                                              ), stringsAsFactors=FALSE)

#aqui es donde se realiza el analisis final como tal, donde se toman las muestras como las marcamos
#se pide el mapeo, con base en la funci�n que se coloco arriba, se asociara con base a la 
#distribuci�n de genes. Un detective es un grupo de kallistos. Con sleuth_prep se pueden
#usar los resultados obtenidos en kallisto, y operarlos, considerando las covariables
# la profundidad de secuenciaci�n, la variaci�n t�cnica y biol�gica. tarjet_mapping 
#ayuda a que laws columnas que tienen letras se lean como caracteres. y con extra_boostrap_
#summary se pueden resumir estadisticas para los recuentos estimados.
so <- sleuth_prep(s2c, ~muestras, target_mapping = t2g,extra_bootstrap_summary = TRUE)

# Con sleuth_fit se estima la varianza t�cnica a partir de los boostraps, 
#la varianza biol�gica y se estima la contracci�n.
so <- sleuth_fit(so)
#sleuth_wt calcula la prueba de Wald en un coeficiente 'beta' espec�fico en cada 
#transcripci�n (NO S� QUE ES ESTO).
so <- sleuth_wt(so, which_beta="muestrasL1")   
#Para poder visualizarlo de una manera agradable. Salimos oprimiendo ESC en R
sleuth_live(so) 

#descargo la tabla que se obtiene en al abrir sleuth_live
setwd("~/6to semestre/Genomica funcional/Archivo/")

resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE) # hacemos una tabla en donde se agregaran los resultados
significativos<-which(resultados$qval<0.1) #se agregaran los que se consideren como significativos, donde el p valor sea menos a .1
significativos<-resultados[significativos,]
upregulated<-which(significativos$b>0) # De los cuales se "marcaran" con que tienen una baja en la expresi�n, si el logforchange es menor a cero
upregulated<-significativos[upregulated,]
downregulated<-which(significativos$b<0) # y si el logforchange es mayor a cero se guardara como los que tienen una alza en su expresi�n
downregulated<-significativos[downregulated,]

#Se separa en dos tablas por aparte, una donde esten los genes significativos con una menor expresi�n y otra en donde estan 
#los genes significativos con una mayor expresi�n, y la tabla se guarda en donde indico que se me guarde
write.table(upregulated,file="~/6to semestre/Genomica funcional/Archivo/Upregulated_N2vsCyg-25.txt",sep="\t")
write.table(downregulated,file="~/6to semestre/Genomica funcional/Archivo/Downregulated_N2vsCyg-25.txt",sep="\t")

upregulated$ext_gene # de esta manera me saldra el nombre externo que esta en la tabla en donde
# estan los genes con una baja expresi�n y que es significativa. Que si nos guiamos por el 
# volcano, serian los que estan en el cuadrante izquierdo superior.
downregulated$eext_gene # En el primero me sale cero caracteres y en este segundo me arroja NULL

