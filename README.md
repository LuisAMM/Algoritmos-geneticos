# Algoritmos-geneticos
Este proyecto consistió en dimensionar un sistema celular de acuerdo con las características de un mapa genérico proporcionado por el usuario, para que, por medio de algoritmos genéticos, sea posible poder proveer la solución para obtener la mayor cobertura posible en dicha región, con distintos tipos de antenas, las cuales pueden variar según ciertas características, minimizando los costos del sistema.


El programa se compila con 


	g++ proyecto.cpp -o proyecto -lm


se necesitan dos parametros para correrlo, numero de vertices que contenga el mapa y el nombre del archivo en el cual se encuentran dichos vertices, para correrlo se necesitan dichos parametros en el orden correcto


	proyecto 25 mapa.txt


El mapa.txt tiene que estar en lo posible centrado en 0


se incluyo una funcion en la libreria clase_ag la cual consiste en inicializar una poblacion en las coordenadas 0,0 y siempre presentes. 
