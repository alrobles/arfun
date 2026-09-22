# Lean en arfun: manual de la línea experimental

Fecha de apertura: 18 de septiembre de 2026.

Este subproyecto verifica partes de la **especificación matemática** de `arfun`. Es una línea paralela: no cambia los ajustes de R/Stan, no necesita GBIF ni un clúster, y no es una dependencia para instalar el paquete R.

## 1. Empezar con un solo comando

Desde la raíz de este repositorio:

```bash
python3 formal/check.py
```

En la primera ejecución de una copia nueva, o si faltan archivos precompilados de mathlib:

```bash
python3 formal/check.py --bootstrap
```

`--bootstrap` descarga la caché oficial correspondiente a los imports de la revisión fijada de mathlib; no selecciona automáticamente una versión más nueva. Necesita conexión a internet y espacio para las dependencias. Las ejecuciones posteriores normalmente reutilizan esa caché.

En esta instalación también puedes ejecutar, desde cualquier directorio:

```bash
python3 /home/alrobles/GitHub/arfun-devel/formal/check.py
```

La primera integración fue comprobada con:

- Lean **4.33.1**.
- mathlib **v4.33.1**, fijada al commit `0df444a360eaa60ab8c11dca51a86af692955474`.
- Python **3.12**; el verificador requiere **Python 3.11 o posterior** por `tomllib`.
- Linux x86-64.
- **26 lemas/teoremas**, **5 módulos**, **23 pruebas del verificador** y **3 controles negativos**.

La señal de éxito es la salida final del comando, no una captura antigua del editor:

```text
VERIFICADO: 26 teoremas; 5 módulos; 3 controles negativos rechazados.
Alcance: especificación matemática; no valida datos ecológicos ni todo el ejecutable Stan.
```

El número de trabajos que muestra Lake incluye miles de módulos de dependencias: no es el número de teoremas propios de arfun.

El informe de la última ejecución está en `.lake/verification-report.json`, relativo a esta carpeta. Incluye estado, versiones, revisiones de dependencias, axiomas, módulos revisados, controles negativos y hashes SHA-256 de las fuentes Lean. Es un archivo generado e ignorado por Git; no sustituye volver a ejecutar la verificación después de editar.

## 2. Qué es Lean, explicado sin asumir experiencia

En R escribimos un cálculo y lo ejecutamos sobre datos. En Lean también podemos escribir una **afirmación matemática** y construir una prueba que un núcleo pequeño, el *kernel*, comprueba.

Una comprobación numérica puede verificar una identidad para ciertos valores. Una prueba formal puede establecerla para **todos los valores que satisfacen sus hipótesis**. Son garantías diferentes y complementarias.

Los componentes del entorno son:

| Componente | Función |
|---|---|
| **Lean** | Lenguaje y comprobador de pruebas. |
| **Elan** | Selecciona e instala versiones de Lean. |
| **Lake** | Gestiona dependencias y compila el proyecto. |
| **mathlib** | Biblioteca matemática: reales, matrices, logaritmos, desigualdades y tácticas. |
| **VS Code + extensión Lean 4** | Editor opcional para ver objetivos y errores mientras escribes. |
| **`check.py`** | Comando de este proyecto que une compilación, auditoría, revalidación y controles. |

Una táctica como `ring` no es una suposición. Construye una prueba de una identidad polinómica, que después comprueba Lean.

**La garantía siempre depende del enunciado.** Una prueba correcta de una afirmación equivocada o demasiado débil no responde a nuestra pregunta científica. Por eso el registro documenta tanto la afirmación como su alcance.

## 3. Qué está demostrado en este primer hito

El catálogo completo está en [theorems.json](theorems.json). Cada entrada identifica el nombre exacto de la declaración, su significado, sus limitaciones y las etiquetas de ecuaciones del [suplemento matemático](../EVALUATIONS/fundamental-niche-evolution-supplement.tex).

| Módulo | Resultados incluidos | Límite importante |
|---|---|---|
| [Geometry.lean](ArfunFormal/Geometry.lean) | Determinante; forma cuadrática positiva; `Matrix.PosDef`; `L L' = Sigma`; validez de `exp(eta)` y `tanh(zeta)`. | Bivariado y sobre números reales exactos, no sobre redondeo numérico. |
| [Likelihood.lean](ArfunFormal/Likelihood.lean) | Cancelación del normalizador; positividad de masas finitas; constante `-n log(area/m)`; cancelación de constantes comunes; identidad de traslación de log-sum-exp. | No demuestra que el fondo represente M ni que la cuadratura aproxime bien una integral. |
| [Quadrature.lean](ArfunFormal/Quadrature.lean) | Cota del error logarítmico y del log-score bajo un error relativo de masa acotado. | La cota del error relativo es una hipótesis que debe justificarse por separado. |
| [Evolution.lean](ArfunFormal/Evolution.lean) | Semidefinitud de `B diag(longitudes) B'`; identidad matricial no centrada; varianza del contraste BM; cambio de unidad temporal; simetría, diagonal y no negatividad de entradas OU; identidad entre raíces fija y estacionaria. | No es todavía la formalización completa de las leyes gaussianas, de las SDE o de los límites del proceso. |

En particular:

- `quadratic_completed_square` se refiere a la forma `x' Sigma x`, usada para demostrar positividad. **No es** la forma de Mahalanobis `x' Sigma^-1 x`.
- `branchCovariance_posSemidef` prueba semidefinitud. La positividad **estricta** requiere condiciones adicionales, como rango completo y longitudes apropiadas.
- `ouCovariance_nonneg` prueba el signo de una entrada. Entradas no negativas **no equivalen** a una matriz definida positiva.
- `stationary_root_identity` es una identidad entre expresiones. La interpretación del proceso exige `alpha > 0` y tiempos válidos, aunque algunas identidades algebraicas también tengan un valor totalizado en Lean cuando un denominador es cero.

### Ejemplo aplicado: por qué importa la cuadratura

Para masas positivas `Z` y `Zhat`, suponemos:

$$
0 \leq \varepsilon < 1,\qquad
\left|\frac{\widehat Z}{Z}-1\right|\leq\varepsilon.
$$

Manteniendo fijo el numerador del likelihood, está demostrado que:

$$
|\widehat\ell-\ell|\leq n[-\log(1-\varepsilon)].
$$

El resultado formal se llama:

```text
ArfunFormal.Quadrature.logScore_error_bound
```

Con muchos registros, un error pequeño de integración puede afectar el log-score. Es una cota para una evaluación de la log-verosimilitud con el numerador fijado, no una cota ya demostrada para ELPD o estimadores posteriores. Lean **no ha demostrado** que nuestro fondo empírico tenga, por ejemplo, un error relativo menor que 1%. La fórmula es una implicación condicionada, no una certificación de los datos.

## 4. Qué no significa un resultado verde

Esta integración no demuestra:

- que los colibríes carezcan de señal filogenética;
- que el nicho fundamental sea exactamente una elipse;
- que se hayan eliminado las interacciones bióticas o el sesgo de observación;
- que la M elegida sea correcta;
- que los datos sean independientes espacialmente;
- que HMC haya convergido o PSIS-LOO sea fiable;
- que el posterior empírico sea propio;
- que el compilador de Stan o todo su ejecutable sean formalmente correctos;
- que `exp`, `tanh`, Cholesky o log-sum-exp tengan error despreciable en punto flotante.

La correspondencia con el código Stan se documenta en `source_map` dentro del registro. **Es una correspondencia revisada, no una traducción formalmente certificada del ejecutable.** Un cambio posterior de Stan exige revisar esa correspondencia aunque las pruebas matemáticas sigan pasando.

Mantenemos cuatro niveles de evidencia:

1. **Lean:** la conclusión se sigue de las hipótesis escritas.
2. **Pruebas numéricas:** implementaciones evaluadas en casos específicos concuerdan dentro de tolerancias.
3. **Simulaciones:** recuperación, cobertura e identificabilidad bajo escenarios controlados.
4. **Validación empírica:** transferencia espacial y por especie, adecuación ecológica y comparación de modelos.

## 5. Organización de archivos

```text
formal/
  lean-toolchain                 versión exacta de Lean
  lakefile.toml                  biblioteca e importación fijada de mathlib
  lake-manifest.json             versiones exactas de todas las dependencias
  ArfunFormal.lean               importa los cuatro módulos de investigación
  ArfunFormal/
    Geometry.lean
    Likelihood.lean
    Quadrature.lean
    Evolution.lean
  theorems.json                  catálogo de afirmaciones y alcance
  check.py                       comando único de verificación
  tests/test_check.py            pruebas del propio verificador
  MANUAL.es.md                   este manual
  .lake/                        dependencias, compilados e informe generados
```

No edites `.lake/packages/mathlib` para hacer pasar una prueba: es una dependencia fijada, no código propio. No subas `.lake/` a Git.

La carpeta `formal/` y el workflow de GitHub están excluidos de la distribución R por `.Rbuildignore`. Usar `arfun` desde R no ejecuta Lean y no exige instalarlo.

## 6. Instalación y editor

En la máquina donde se abrió esta línea ya estaban instalados Lean, Lake y Elan. No se cambió la versión global de Lean: `lean-toolchain` selecciona la versión del subproyecto.

En otra máquina:

1. Instala Git, Python 3.11 o posterior y el entorno Lean/Elan siguiendo la [guía oficial](https://lean-lang.org/install/).
2. Conserva `lean-toolchain`, `lakefile.toml` y `lake-manifest.json` del repositorio. No sustituyas las versiones por `latest` o `stable`.
3. Desde la raíz del repositorio, ejecuta `python3 formal/check.py --bootstrap`.
4. Comprueba el resultado final y el informe.

La receta fue ejecutada localmente en Linux. El workflow está preparado para Ubuntu 24.04; otros sistemas requieren verificar sus prerrequisitos. Para Windows, WSL permite seguir el flujo Linux.

Para aprender con ayuda visual:

1. Instala VS Code y la [extensión oficial Lean 4](https://marketplace.visualstudio.com/items?itemName=leanprover.lean4).
2. Abre **la carpeta `formal`** como carpeta de trabajo.
3. Abre `ArfunFormal/Geometry.lean`.
4. Coloca el cursor dentro de una prueba y consulta el panel de objetivos, *Infoview*.
5. Usa la navegación a definición para consultar los lemas de mathlib.

No es necesario instalar la extensión para ejecutar el verificador desde la terminal.

## 7. Cómo leer una prueba real del proyecto

En `Evolution.lean` aparece:

```lean
theorem bm_contrast_variance (v ti tj shared : ℝ) :
    v * ti + v * tj - 2 * (v * shared) = v * (ti + tj - 2 * shared) := by
  ring
```

Lectura:

- `theorem`: vamos a demostrar una afirmación.
- `bm_contrast_variance`: nombre para reutilizarla.
- `(v ti tj shared : ℝ)`: para cualesquiera cuatro números reales.
- Lo situado después de `:` es la conclusión.
- `:= by` inicia la construcción de la prueba.
- `ring` comprueba la identidad polinómica y produce una prueba.

Esta prueba no establece por sí sola que los rasgos sigan BM. Establece el paso algebraico que se usa **después** de especificar las varianzas y covarianzas BM.

Un resultado con hipótesis, en `Geometry.lean`, tiene esta estructura:

```lean
theorem covariance_posDef {s₁ s₂ ρ : ℝ}
    (h₁ : 0 < s₁) (h₂ : 0 < s₂) (hρ : |ρ| < 1) :
    (covariance s₁ s₂ ρ).PosDef := by
```

La prueba completa está en el archivo. `h₁`, `h₂` y `hρ` son evidencias de las condiciones del teorema, no estimaciones obtenidas de GBIF. La conclusión solo se aplica bajo esas condiciones.

### Vocabulario mínimo de tácticas

| Expresión | Lectura práctica |
|---|---|
| `exact h` | La evidencia `h` resuelve el objetivo actual. |
| `apply resultado` | Reduce el objetivo a las hipótesis de ese resultado. |
| `have h : P := ...` | Demuestra y nombra un paso intermedio. |
| `rw [identidad]` | Reescribe usando una igualdad demostrada. |
| `simp [...]` | Simplifica con reglas justificadas. |
| `ring` | Resuelve una identidad polinómica. |
| `linarith` / `nlinarith` | Usa consecuencias de igualdades/desigualdades aritméticas. |
| `constructor` | Divide una conclusión compuesta en sus partes. |
| `fin_cases` | Revisa todos los casos de un tipo finito, aquí índices de una matriz 2×2. |

No todas las afirmaciones verdaderas se resuelven automáticamente con estas tácticas. Un error de prueba puede significar que faltan pasos, imports o hipótesis; no es automáticamente un contraejemplo al teorema.

## 8. Primer ejercicio, separado de la biblioteca auditada

Crea con tu editor un archivo `Exercises.lean` dentro de `formal/`:

```lean
module
public import ArfunFormal

example (v t : ℝ) : v * t + v * t = 2 * v * t := by
  ring
```

Desde `formal/`, ejecútalo:

```bash
lake env lean Exercises.lean
```

La ausencia de errores y salida cero indican que el ejemplo fue comprobado. Este ejemplo se verificó al preparar la integración.

Como ejercicio, cambia deliberadamente el `2` por `3`. La identidad general ya no es cierta y la prueba debe fallar. No cambies los teoremas de investigación ni sus controles para realizar el ejercicio.

`Exercises.lean` **no pertenece automáticamente a la biblioteca auditada**. El comando principal comprueba los módulos `ArfunFormal` y el registro, no cualquier archivo de ejercicios. Para promover una prueba al proyecto, sigue el siguiente apartado.

## 9. Añadir un resultado de investigación

1. Formula primero el enunciado en el suplemento o en tus notas, distinguiendo hipótesis y conclusión.
2. Busca si ya existe un lema equivalente en la revisión de mathlib instalada.
3. Añade el resultado al módulo apropiado de `ArfunFormal/`.
4. Usa el único `namespace` del módulo. Las declaraciones del piloto siguen el formato `theorem nombre` o `lemma nombre` al inicio de línea; esto permite comprobar su cobertura en el registro.
5. Si creas un módulo, añádelo explícitamente a `ArfunFormal.lean`.
6. Añade una entrada a `theorems.json` con nombre completo, `claim_es`, `scope_es` y las etiquetas pertinentes del suplemento.
7. Ejecuta `python3 formal/check.py` desde la raíz del repositorio.
8. Revisa el enunciado con criterio matemático y biológico. No debilites sus hipótesis o cambies su significado solo para conseguir que compile.

Por ejemplo, el nombre completo del teorema central de geometría es:

```text
ArfunFormal.Geometry.covariance_posDef
```

`autoImplicit=false` evita crear inadvertidamente ciertos parámetros implícitos no declarados. Las advertencias se consideran errores. El registro también debe actualizarse cuando se añadan lemas auxiliares públicos.

Las definiciones `noncomputable` expresan objetos matemáticos exactos, como funciones reales, que no son necesariamente algoritmos ejecutables de cálculo numérico. No significan que una prueba esté incompleta.

## 10. Qué hace exactamente el verificador

`check.py` realiza estas comprobaciones:

1. Versiones fijadas coherentes entre configuración, registro y lockfile; dependencias resueltas a commits completos.
2. Correspondencia entre declaraciones de las fuentes propias y el catálogo; imports explícitos de los módulos.
3. Ausencia en esas fuentes de construcciones prohibidas por la política del piloto, como `sorry`, `admit`, axiomas ad hoc o `native_decide`.
4. Pruebas unitarias del parser de auditoría y de los controles de cobertura/versiones.
5. Compilación con `lake --wfail build`.
6. Consulta de `#print axioms` para **cada teorema registrado**, incluyendo sus dependencias transitivas.
7. Revalidación de los cinco módulos propios mediante `leanchecker`.
8. Tres controles negativos reales enviados a Lean por entrada estándar, sin incorporarlos a la biblioteca:
   - una prueba falsa con `rfl` debe ser rechazada por Lean;
   - una prueba incompleta expone `sorryAx` y debe ser rechazada por la auditoría;
   - una prueba que depende de un axioma inventado debe ser rechazada por la auditoría.
9. Escritura del informe con estado `passed` solamente al completar todo el recorrido.

El conjunto permitido de axiomas fundacionales es:

```text
propext
Classical.choice
Quot.sound
```

Son fundamentos estándar usados por Lean/mathlib. No son supuestos biológicos sobre el nicho. Decir «sin `sorry` ni axiomas ad hoc» **no significa** «sin fundamentos lógicos».

Los controles negativos no prueban nuevas afirmaciones científicas: prueban que el mecanismo de verificación detecta fallos. Su rechazo es el resultado esperado.

`leanchecker` vuelve a comprobar declaraciones usando el kernel de Lean. Es una defensa adicional frente a problemas del entorno; **no es un verificador externo independiente**. Las bibliotecas importadas, el kernel y su entorno de ejecución siguen formando parte de la base de confianza. El escaneo textual de fuentes tampoco sustituye una auditoría de código no confiable.

## 11. Comandos útiles

Todos los comandos siguientes se ejecutan desde `formal/`, salvo que se indique otra cosa:

```bash
lean --version
lake --version
lake --wfail build
lake env leanchecker -v ArfunFormal
python3 check.py
```

Para revisar un módulo mientras trabajas:

```bash
lake lean ArfunFormal/Geometry.lean
```

Para ejecutar solo las pruebas Python del verificador:

```bash
python3 -m unittest discover -s tests -v
```

Para la comprobación formal completa, usa siempre `check.py`. Compilar un archivo aislado no reemplaza la auditoría del conjunto.

## 12. Problemas comunes

| Mensaje o situación | Qué revisar |
|---|---|
| `lake` no encontrado | Instalación de Elan y su directorio de ejecutables en `PATH`. |
| `No default toolchain` | Ejecuta los comandos Lake desde `formal/`; el script principal cambia a esa carpeta automáticamente. |
| Python no encuentra `tomllib` | Usa Python 3.11 o posterior. |
| Falta un archivo `.olean` | Prueba `python3 check.py --bootstrap`; verifica la conexión y los imports. |
| `Unknown identifier` | Import correcto y nombre en la revisión fijada de mathlib, no necesariamente en la documentación de la versión más reciente. |
| `unsolved goals` | Abre Infoview y revisa qué paso sigue sin demostrar. |
| Argumento de `simp` no utilizado | Limpia la prueba; no desactives las advertencias para ocultarlo. |
| `Registro y fuentes difieren` | Falta registrar una declaración, sobra una entrada o no coincide el nombre completo. |
| `sorryAx` o axioma no permitido | La prueba o una dependencia está incompleta o usa una suposición ajena a la política. No amplíes la lista permitida para hacerla pasar. |
| El informe anterior dice `passed` | Vuelve a ejecutar el verificador: ese informe corresponde a unas fuentes y versiones concretas. |

No uses una versión aleatoria de mathlib, no edites las dependencias descargadas, no desactives el kernel y no borres resultados científicos para solucionar un problema de compilación.

## 13. GitHub Actions y distribución R

El workflow [lean.yml](../.github/workflows/lean.yml) está preparado para cambios relevantes en `formal/`, Stan o el propio workflow, y para ejecución manual desde GitHub Actions.

- El checkout está fijado a un commit de `actions/checkout` v4.2.2.
- El instalador Elan usado por CI está fijado a **v4.2.4** y se verifica con SHA-256.
- Lean y mathlib se seleccionan mediante los archivos del subproyecto.
- El trabajo ejecuta el mismo `python3 formal/check.py --bootstrap` que usamos localmente.
- El token del workflow tiene permisos de lectura; no publica ni modifica resultados científicos.

Crear el archivo no equivale a haber ejecutado un trabajo remoto. La primera integración se validó localmente; el workflow solo podrá ejecutarse en GitHub cuando los cambios se incorporen al repositorio remoto. No se hizo commit ni push al abrir esta línea.

El filtro de Stan sirve para pedir una nueva revisión cuando cambie el modelo. Un workflow verde **no certifica automáticamente** que una fórmula modificada en Stan siga correspondiendo a la especificación formal.

## 14. Cómo continúa la línea experimental

El alcance implementado corresponde a un primer hito de álgebra, geometría y cotas. El trabajo siguiente debe mantenerse explícitamente separado de los resultados ya verificados:

1. **F2: estructura filogenética.** Construir la matriz de caminos desde un árbol finito, demostrar positividad estricta bajo condiciones precisas y formalizar la transformación de vectores gaussianos.
2. **F3: OU.** Derivar la covarianza desde transiciones por ramas, demostrar el límite BM y distinguir formalmente los dos límites de selección intensa.
3. **F4: reconstrucción e identificación.** Condicionamiento gaussiano ancestral, equivalencias entre filtros y respuesta ambiental, y condiciones de normalización del posterior.
4. **F5: cálculo estocástico.** Evaluar las bibliotecas disponibles antes de formalizar Itô y el puente microevolutivo completo; no reconstruir todo desde cero como requisito inicial.
5. **Interfaz espacial.** Formalizar invariantes de particiones y separación, sin confundir folds disjuntos con independencia estadística. El k-fold geográfico sigue siendo un programa empírico separado y todavía no está implementado por este subproyecto.

Los contraejemplos y las hipótesis faltantes son resultados valiosos de esta línea. El objetivo no es demostrar que nuestras conclusiones anteriores siempre eran correctas, sino establecer con precisión cuáles se sostienen y bajo qué condiciones.

## 15. Recursos para aprender

- [Instalación oficial de Lean](https://lean-lang.org/install/).
- [Referencia del lenguaje: validar una prueba](https://lean-lang.org/doc/reference/latest/ValidatingProofs/).
- [Lake y gestión de proyectos](https://lean-lang.org/doc/reference/latest/Build-Tools-and-Distribution/Lake).
- [Axiomas y computación en Lean](https://lean-lang.org/theorem_proving_in_lean4/Axioms-and-Computation/).
- [Matrices positivas en mathlib](https://leanprover-community.github.io/mathlib4_docs/Mathlib/LinearAlgebra/Matrix/PosDef.html).
- [Gaussianas multivariadas en mathlib](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Probability/Distributions/Gaussian/Multivariate.html).

La documentación web puede avanzar más rápido que nuestra versión fijada. Para dudas sobre un nombre o una hipótesis exacta, consulta también la fuente instalada dentro de `.lake/packages/mathlib/` sin modificarla.
