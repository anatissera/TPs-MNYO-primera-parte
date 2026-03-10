## Repositorio de la primera parte de la materia Métodos Numéricos y Optimización

### Trabajo Práctico 1: Análisis de esquemas de interpolación y solución de sistemas no lineales
La interpolación es un método de estimación a partir del cual se construyen funciones dado un conjunto de datos discretos conocidos o DataFrame para aproximar aquellos que no. En otras palabras, hallar alguna aproximación de la función que podría estar generando los datos.
Este trabajo analiza el rendimiento de diversos esquemas de interpolación al igual que la reconstrucción trayectorias. En primera instancia, se encuentra que el método de Splines Cúbicos minimiza el error con respecto de la aproximación de funciones de  forma lineal o polinomial. Además, proponemos aprovechar el teorema de Chebyshev para mejorar el funcionamiento de los modelos polinomiales. 

Luego, teniendo en cuenta estas conclusiones, utilizamos este método para reconstruir las trayectorias dadas y escribimos un algoritmo de Newton-Raphson para hallar soluciones de sistemas no lineales.

Desarrollo y conclusiones en [MNYO-TP1.pdf](TP1/MNYO-TP1.pdf)

### Trabajo Práctico 2: Modelos de dinámica poblacional

Este trabajo aborda el estudio de diversos modelos que describen la dinámica poblacional de una o varias especies dentro de un sistema bajo distintas limitaciones.
Se exploran tanto el modelo de crecimiento exponencial como el logístico, junto con el modelo de competencia interespecífica de Lotka-Volterra, así como el modelo Predador-Presa y su versión extendida. 
Se analizan las características de cada modelo, sus restricciones y su capacidad para predecir la evolución de las poblaciones en diferentes escenarios.
Para abordar estos modelos, se llevan a cabo experimentos numéricos con el fin de aproximar las ecuaciones diferenciales que los gobiernan. Se comparan las soluciones numéricas con aquellas exactas, cuando están disponibles, y se estima el error en cada instancia. En particular, se determina que el método de Runge-kutta de orden 4 es el más eficiente en términos de precisión.
Adicionalmente, se utilizan diagramas de isoclinas para determinar las distintas dinámicas de los diferentes modelos, lo que permite plantear diversos casos según su comportamiento. En cada uno se analizan asimismo los puntos de equilibrio y las trayectorias de las poblaciones en relación con diversas condiciones iniciales y parámetros del modelo para poder determinar la estabilidad de cada uno.
De esta manera, el informe ofrece una perspectiva detallada sobre varios modelos de crecimiento y competencia poblacional, explorando su comportamiento mediante experimentos numéricos y evaluando su capacidad predictiva en distintos contextos.

Desarrollo en [MNYO-TP2.pdf](TP2/MNYO-TP2.pdf)
