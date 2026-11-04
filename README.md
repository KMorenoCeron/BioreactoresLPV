# BioreactoresLPV

Control LPV con ubicación de polos en región D para maximizar la productividad de acetato en un proceso de fermentación oscura. Se incluyen las carpetas de **SeDuMi** y **ROLMIP** utilizadas para el diseño del controlador LPV.

> **Nota:** La estrategia de control incluida corresponde a una propuesta preliminar que aún no ha sido validada completamente para el sistema. Es necesario verificar su consistencia con el modelo y las condiciones de operación del proceso.

El diseño del controlador es independiente del diseño del filtro basado en $H_\infty$. Para su formulación, solo se consideran las relaciones entre las ecuaciones del modelo y los puntos de operación.

## Descripción de los códigos

- **bloque1_observador_PI** — Diseña el observador PI de la tasa de reacción y lo valida de forma aislada para los polos elegidos. Muestra la convergencia de los estados estimados desde un error inicial.

- **bloque2_optimizador_w** — Desarrolla las relaciones que permiten no depender de $S_{in}$, sino solo de la combinación observable $w = A_{in} + k\,S_{in} = A + k\,S$. Valida el cálculo del punto de operación a partir de $w$.

- **bloque3_diseno_LPV** — Diseño del controlador LPV por LMI con región D. Entrega las ganancias afines $K_0$, $K_1$, $K_2$ y las guarda para los demás códigos.

- **ValidacionCompleta** — Ejecuta todo el esquema en un solo código: observador, optimizador por $w$, scheduling y controlador LPV en lazo cerrado sobre el proceso.

- **ComparacionPertu** — Comprueba cómo responde el esquema ante variaciones de $A_{in}$ y $S_{in}$. Compara el caso con solo $S_{in}$ variable frente al caso con ambas perturbaciones.

- **graficas_comparacion** — Extrae las gráficas del paper de referencia y las compara con el diseño propuesto.

## Requisitos

- MATLAB con YALMIP, ROLMIP y SeDuMi en el path, para los códigos de diseño (bloque3).
- Los códigos de validación usan las ganancias guardadas.
